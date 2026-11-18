"""Postulat camkii — port Python fidèle de omegaModel.js (banc_essai_v2).

Unités : concentrations en µM, temps en s, kon en µM⁻¹·s⁻¹,
koff / kcat / taux effectifs en s⁻¹.

Statut épistémique :
- Le modèle est un modèle cinétique déterministe (masse conservée,
  complexes explicites, sous-pas internes).
- La lecture « PP1 / CaMKII / I-1 / calcineurine » reste une analogie
  structurelle avec la physiologie, pas une validation.

Note de portage :
- `bind()` conserve l'approximation de relaxation linéarisée du JS
  d'origine (taux forward calculé sur les concentrations libres
  *courantes*, pas à l'équilibre). C'est stable et borné, mais la
  convergence n'est pas exactement quadratique ; c'est documenté ici
  et non corrigé, pour rester fidèle à la référence.
"""

from __future__ import annotations

import math

# ── Métadonnées du postulat ────────────────────────────────────────────────

POSTULAT = "camkii"
STATUT = (
    "modèle numérique ; lecture physiologique (PP1 / CaMKII / I-1 / "
    "calcineurine) = analogie structurelle, jamais équivalence ontologique"
)
CRITERE = (
    "Invariants de masse (cam, camkii, pp1, i1) conservés à 1e-3 près "
    "sur toute la simulation ; aucune espèce négative ; protocole "
    "effectivement déclenché (Ca max > 0). "
    "Échec si dérive > 1e-3, si une espèce devient négative, "
    "ou si aucune activité n'est observée."
)

# ── Constantes ─────────────────────────────────────────────────────────────

K_EXCHANGE = 1.0
K_I1_PHOS = 2.0
K_I1_DEPHOS = 0.15
K_RECOVERY = 0.02

TOTALS = {"cam": 10.0, "camkii": 50.0, "pp1": 2.0, "i1": 1.0}

DEFAULT_PARAMS: dict = {
    "calmodulin": {
        "kon1": 5.0,
        "koff1": 50.0,
        "kon2": 100.0,
        "koff2": 2000.0,
    },
    "camkii": {
        "konBind": 1.27e-3,
        "koffBind": 0.07,
        "koffBoundPhosphorylated": 0.01,
        "kPhosphorylation": 0.5,
    },
    "dephosphorylation": {"kon": 1.0, "koff": 0.001, "kcat": 0.41},
    "pp1i1": {"kon": 1.0, "koff": 0.001},
    "pkaActivity": 0.05,
    "calcineurinCa50": 5.0,
    "calcineurinHill": 2.0,
}

INITIAL_STATE: dict = {
    "M0": 10.0, "M2N": 0.0, "M2C": 0.0, "M4": 0.0,
    "K0": 50.0, "KC": 0.0, "KpC": 0.0, "Kp": 0.0, "D": 0.0,
    "PP1": 2.0, "I1": 1.0, "I1p": 0.0, "PP1I1p": 0.0, "PP1KpC": 0.0,
}

DEFAULT_PROTOCOL: dict = {
    "mode": "train",
    "amplitude": 20.0,
    "basal": 0.1,
    "width": 2.0,
    "period": 20.0,
    "delay": 5.0,
}

EPS = 1e-12
MAX_INTERNAL_DT = 0.005  # s


# ── Entrée calcique ────────────────────────────────────────────────────────

def calcineurin_activity(ca: float, p: dict) -> float:
    """Activité calcineurine dépendante du calcium (Hill simple, bornée)."""
    ca50 = max(p["calcineurinCa50"], EPS)
    hill = max(p["calcineurinHill"], EPS)
    num = max(ca, 0.0) ** hill
    den = ca50 ** hill + num
    return num / den if den > EPS else 0.0


def calcium_at(t: float, p: dict) -> float:
    """Entrée calcique. Modes : constant | impulse | train."""
    mode = p.get("mode", "train")
    amp = max(p.get("amplitude", 0.0), 0.0)
    basal = max(p.get("basal", 0.0), 0.0)
    delay = p.get("delay", 0.0)
    width = p.get("width", 0.0)

    if mode == "constant":
        return amp
    if t < delay:
        return basal
    if mode == "impulse":
        return amp if t < delay + width else basal

    cycle = max(p.get("period", 20.0), width + 0.01)
    phase = (t - delay) % cycle
    return amp if phase < width else basal


# ── Briques cinétiques ─────────────────────────────────────────────────────

def _clamp_nonneg(x: float) -> float:
    return x if x > 0.0 else 0.0


def relax(s: dict, a: str, b: str, ka: float, kb: float, dt: float) -> None:
    """A <-> B — solution exacte de dA/dt = -ka·A + kb·B. Conserve A+B."""
    total = s[a] + s[b]
    if total <= EPS:
        s[a] = 0.0
        s[b] = 0.0
        return
    k = ka + kb
    if k <= EPS:
        return
    eq_b = (ka / k) * total
    nxt_b = eq_b + (s[b] - eq_b) * math.exp(-k * dt)
    s[b] = _clamp_nonneg(nxt_b)
    s[a] = _clamp_nonneg(total - s[b])


def bind(
    s: dict, a: str, b: str, c: str,
    kf: float, kr: float, dt: float,
) -> None:
    """A + B <-> C — équilibre quadratique + relaxation linéarisée.

    NOTE de portage : le taux forward est calculé sur les concentrations
    libres *courantes* (a0, b0), pas à l'équilibre. Identique au JS.
    """
    a0, b0, c0 = s[a], s[b], s[c]
    total = a0 + b0 + c0
    if total <= EPS or (kf <= EPS and kr <= EPS):
        return

    ta, tb = a0 + c0, b0 + c0
    c_eq = c0

    if kf > EPS:
        qa = kf
        qb = -(kf * (ta + tb) + kr)
        qc = kf * ta * tb
        disc = max(qb * qb - 4.0 * qa * qc, 0.0)
        root = math.sqrt(disc)
        r1 = (-qb + root) / (2.0 * qa)
        r2 = (-qb - root) / (2.0 * qa)
        upper = min(ta, tb)

        if -EPS <= r1 <= upper + EPS:
            c_eq = r1
        elif -EPS <= r2 <= upper + EPS:
            c_eq = r2
        else:
            c_eq = min(max(r2, 0.0), upper)
    elif kr > EPS:
        c_eq = 0.0

    fa, fb = max(a0, 0.0), max(b0, 0.0)
    k = kf * (fa + fb) + kr
    if k <= EPS:
        return

    nxt_c = c_eq + (c0 - c_eq) * math.exp(-k * dt)
    shift = _clamp_nonneg(nxt_c) - c0

    s[c] = _clamp_nonneg(c0 + shift)
    s[a] = _clamp_nonneg(a0 - shift)
    s[b] = _clamp_nonneg(b0 - shift)


def transfer(s: dict, src: str, dst: str, k: float, dt: float) -> None:
    """X -> Y — décroissance exponentielle exacte. Conserve X+Y."""
    if k <= EPS or s[src] <= EPS:
        return
    moved = s[src] * (1.0 - math.exp(-k * dt))
    s[src] = _clamp_nonneg(s[src] - moved)
    s[dst] = _clamp_nonneg(s[dst] + moved)


# ── Pas d'intégration ──────────────────────────────────────────────────────

def step(s: dict, p: dict, ca: float, dt: float) -> None:
    """Un pas, avec sous-pas si dt > MAX_INTERNAL_DT."""
    if not math.isfinite(dt) or dt <= 0.0:
        return
    n = max(1, math.ceil(dt / MAX_INTERNAL_DT))
    h = dt / n
    for _ in range(n):
        step_once(s, p, ca, h)


def step_once(s: dict, p: dict, ca: float, dt: float) -> None:
    c2 = ca * ca
    cm = p["calmodulin"]

    # Calmoduline : M0 <-> M2C <-> M4 <-> M2N <-> M0
    relax(s, "M0", "M2C", cm["kon2"] * c2, cm["koff2"], dt)
    relax(s, "M2C", "M4", cm["kon1"] * c2, cm["koff1"], dt)
    relax(s, "M4", "M2N", cm["koff2"], cm["kon2"] * c2, dt)
    relax(s, "M2N", "M0", cm["koff1"], cm["kon1"] * c2, dt)
    relax(s, "M2C", "M2N", K_EXCHANGE, K_EXCHANGE, dt)

    # I-1 : I1 <-> I1p ; I1p + PP1 <-> PP1I1p
    i1_dephos = K_I1_DEPHOS * calcineurin_activity(ca, p)
    relax(s, "I1", "I1p", K_I1_PHOS * p["pkaActivity"], i1_dephos, dt)
    bind(s, "I1p", "PP1", "PP1I1p",
         p["pp1i1"]["kon"], p["pp1i1"]["koff"], dt)

    # CaMKII : K0 + CaM_active <-> KC -> KpC -> Kp
    ck = p["camkii"]
    cam_active = s["M2C"] + s["M4"]
    relax(s, "K0", "KC", ck["konBind"] * cam_active, ck["koffBind"], dt)
    transfer(s, "KC", "KpC", ck["kPhosphorylation"], dt)
    transfer(s, "KpC", "Kp", ck["koffBoundPhosphorylated"], dt)

    # Déphosphorylation : Kp + PP1 <-> PP1KpC -> D + PP1
    dp = p["dephosphorylation"]
    bind(s, "Kp", "PP1", "PP1KpC", dp["kon"], dp["koff"], dt)
    if s["PP1KpC"] > EPS:
        done = s["PP1KpC"] * (1.0 - math.exp(-dp["kcat"] * dt))
        s["PP1KpC"] = _clamp_nonneg(s["PP1KpC"] - done)
        s["PP1"] = _clamp_nonneg(s["PP1"] + done)
        s["D"] = _clamp_nonneg(s["D"] + done)

    # Recovery : D -> K0
    transfer(s, "D", "K0", K_RECOVERY, dt)


# ── Invariants ─────────────────────────────────────────────────────────────

def invariants_of(s: dict) -> dict:
    return {
        "cam": s["M0"] + s["M2N"] + s["M2C"] + s["M4"],
        "camkii": (
            s["K0"] + s["KC"] + s["D"]
            + s["KpC"] + s["Kp"] + s["PP1KpC"]
        ),
        "pp1": s["PP1"] + s["PP1I1p"] + s["PP1KpC"],
        "i1": s["I1"] + s["I1p"] + s["PP1I1p"],
    }


# ── Boucle de simulation ───────────────────────────────────────────────────

def simulate(
    params: dict | None = None,
    protocol: dict | None = None,
    duration: float = 100.0,
    dt: float = 0.01,
    max_points: int = 520,
    initial_state: dict | None = None,
    check_invariants: bool = False,
) -> dict:
    """Simule le système ; retourne séries + diagnostic."""
    if not math.isfinite(duration) or duration <= 0.0:
        raise ValueError("duration doit être un nombre positif.")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt doit être un nombre positif.")
    if not isinstance(max_points, int) or max_points <= 0:
        raise ValueError("max_points doit être un entier positif.")

    p = _deep_update(DEFAULT_PARAMS, params or {})
    proto = {**DEFAULT_PROTOCOL, **(protocol or {})}
    s = {**INITIAL_STATE, **(initial_state or {})}

    expected = invariants_of(s)
    rows: list[dict] = []
    max_drift = 0.0
    min_value = float("inf")

    def push(time: float) -> None:
        rows.append({
            "t": round(time, 3),
            "ca": calcium_at(time, proto),
            **s,
        })

    push(0.0)

    approx_steps = math.ceil(duration / dt)
    stride = max(1, math.ceil(approx_steps / max_points))
    t = 0.0
    steps_taken = 0

    while t < duration - 1e-9:
        h = min(dt, duration - t)
        ca = calcium_at(t, proto)

        step(s, p, ca, h)
        t += h
        steps_taken += 1

        if check_invariants:
            cur = invariants_of(s)
            for k in expected:
                max_drift = max(max_drift, abs(cur[k] - expected[k]))

        for v in s.values():
            if v < min_value:
                min_value = v

        if steps_taken % stride == 0 or t >= duration - 1e-9:
            push(t)

    return {
        "rows": rows,
        "duration": duration,
        "steps_taken": steps_taken,
        "invariants": invariants_of(s),
        "expected_invariants": expected,
        "max_drift": max_drift,
        "min_value": min_value,
        "params": p,
        "protocol": proto,
    }


# ── Entrée standard du moteur ──────────────────────────────────────────────

def executer(
    duration: float = 100.0,
    dt: float = 0.01,
    max_points: int = 520,
    check_invariants: bool = True,
    protocol: dict | None = None,
    params: dict | None = None,
    initial_state: dict | None = None,
    **extra,
) -> dict:
    """Exécute le postulat ; retourne un dict contenant `critere_ok`."""
    if extra:
        # tolère des overrides plats (ex. kcat=0.5) en plus de params={...}
        params = _deep_update(params or {}, extra)

    res = simulate(
        params=params,
        protocol=protocol,
        duration=duration,
        dt=dt,
        max_points=max_points,
        initial_state=initial_state,
        check_invariants=check_invariants,
    )

    ca_series = [r["ca"] for r in res["rows"]]
    kp_series = [r["Kp"] for r in res["rows"]]
    kpc_series = [r["KpC"] for r in res["rows"]]
    i1p_series = [r["I1p"] for r in res["rows"]]
    pp1_series = [r["PP1"] for r in res["rows"]]
    d_series = [r["D"] for r in res["rows"]]

    ca_max = max(ca_series) if ca_series else 0.0
    activite = {
        "ca_max": ca_max,
        "Kp_max": max(kp_series) if kp_series else 0.0,
        "KpC_max": max(kpc_series) if kpc_series else 0.0,
        "I1p_max": max(i1p_series) if i1p_series else 0.0,
        "PP1_min": min(pp1_series) if pp1_series else 0.0,
        "D_max": max(d_series) if d_series else 0.0,
    }

    critere_ok = (
        res["max_drift"] < 1e-3
        and res["min_value"] >= -1e-9
        and ca_max > 0.0
        and res["steps_taken"] > 0
    )

    return {
        "postulat": POSTULAT,
        "statut": STATUT,
        "critere": CRITERE,
        "critere_ok": critere_ok,
        "diagnostic": {
            "max_drift": res["max_drift"],
            "min_value": res["min_value"],
            "steps_taken": res["steps_taken"],
            "invariants": res["invariants"],
            "expected_invariants": res["expected_invariants"],
        },
        "activite": activite,
        "protocole": res["protocol"],
        "parametres": res["params"],
        "series": res["rows"],
    }


# ── Utilitaires ────────────────────────────────────────────────────────────

def _deep_update(base: dict, override: dict) -> dict:
    out = {k: (dict(v) if isinstance(v, dict) else v) for k, v in base.items()}
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_update(out[k], v)
        else:
            out[k] = v
    return out


# ── Test rapide ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import json

    out = executer()
    print(json.dumps(
        {k: out[k] for k in
         ("postulat", "critere_ok", "diagnostic", "activite")},
        indent=2, ensure_ascii=False,
    ))
