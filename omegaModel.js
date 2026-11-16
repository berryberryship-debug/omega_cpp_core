// omegaModel.js — moteur Omega V1 du simulateur
//
// Unités :
// - Concentrations : µM
// - Temps : s
// - kon : µM⁻¹·s⁻¹
// - koff, kcat, taux effectifs : s⁻¹
//
// Principes :
// - Toutes les réactions conservent la masse.
// - Les complexes sont explicites : PP1I1p et PP1KpC.
// - Seule la PP1 libre est catalytiquement active.
// - La déphosphorylation de Kp passe obligatoirement par PP1KpC.

export const K_EXCHANGE = 1.0; // échange Ca entre lobes C et N (s⁻¹)
export const K_I1_PHOS = 2.0; // phosphorylation de I-1 par PKA (µM⁻¹·s⁻¹)
export const K_I1_DEPHOS = 0.15; // déphosphorylation basale de I-1p (s⁻¹)
export const K_RECOVERY = 0.02; // D -> K0 (s⁻¹)

export const TOTALS = {
  cam: 10,
  camkii: 50,
  pp1: 2,
  i1: 1,
};

export const DEFAULT_PARAMS = {
  calmodulin: {
    kon1: 5.0, // M2C + Ca²⁺ -> M4
    koff1: 50.0, // M4 -> M2C + Ca²⁺
    kon2: 100.0, // M0 + 2 Ca²⁺ -> M2C
    koff2: 2000.0, // M2C -> M0 + 2 Ca²⁺
  },
  camkii: {
    konBind: 1.27e-3, // K0 + CaM active -> KC
    koffBind: 0.07, // KC -> K0 ; KpC -> Kp
    kPhosphorylation: 0.5, // KC -> KpC
  },
  dephosphorylation: {
    kon: 1.0, // Kp + PP1 -> PP1KpC
    koff: 0.001, // PP1KpC -> Kp + PP1
    kcat: 0.41, // PP1KpC -> D + PP1
  },
  pp1i1: {
    kon: 1.0, // I1p + PP1 -> PP1I1p
    koff: 0.001, // PP1I1p -> I1p + PP1
  },
  pkaActivity: 0.05, // µM
};

export const INITIAL_STATE = {
  M0: 10.0,
  M2N: 0.0,
  M2C: 0.0,
  M4: 0.0,

  K0: 50.0,
  KC: 0.0,
  KpC: 0.0,
  Kp: 0.0,
  D: 0.0,

  PP1: 2.0,
  I1: 1.0,
  I1p: 0.0,
  PP1I1p: 0.0,
  PP1KpC: 0.0,
};

export const DEFAULT_PROTOCOL = {
  mode: "train",
  amplitude: 20,
  basal: 0.1,
  width: 2,
  period: 20,
  delay: 5,
};

const EPS = 1e-12;
const MAX_INTERNAL_DT = 0.005; // s

const clampNonNegative = (x) => (x > 0 ? x : 0);

/**
 * Entrée calcique.
 * Modes : constant | impulse | train
 */
export function calciumAt(t, p) {
  if (p.mode === "constant") {
    return Math.max(p.amplitude, 0);
  }

  if (t < p.delay) {
    return Math.max(p.basal, 0);
  }

  if (p.mode === "impulse") {
    return t < p.delay + p.width ? Math.max(p.amplitude, 0) : Math.max(p.basal, 0);
  }

  const cycle = Math.max(p.period, p.width + 0.01);
  const phase = (t - p.delay) % cycle;

  return phase < p.width ? Math.max(p.amplitude, 0) : Math.max(p.basal, 0);
}

/**
 * A <-> B.
 * Solution exacte de dA/dt = -ka A + kb B.
 * Conserve A + B.
 */
function relax(s, a, b, ka, kb, dt) {
  const total = s[a] + s[b];

  if (total <= EPS) {
    s[a] = 0;
    s[b] = 0;
    return;
  }

  const k = ka + kb;

  if (k <= EPS) {
    return;
  }

  const equilibriumB = (ka / k) * total;
  const nextB = equilibriumB + (s[b] - equilibriumB) * Math.exp(-k * dt);

  s[b] = clampNonNegative(nextB);
  s[a] = clampNonNegative(total - s[b]);
}

/**
 * A + B <-> C.
 * Solution exacte à l'équilibre, puis relaxation vers cet équilibre.
 * Conserve A + B + C.
 */
function bind(s, a, b, c, kf, kr, dt) {
  const a0 = s[a];
  const b0 = s[b];
  const c0 = s[c];
  const total = a0 + b0 + c0;

  if (total <= EPS || (kf <= EPS && kr <= EPS)) {
    return;
  }

  // Équilibre : kf (Afree)(Bfree) = kr C
  // Avec x = C : kf (a0 - x)(b0 - x) = kr (c0 + x)
  let complexAtEquilibrium = c0;

  if (kf > EPS) {
    const qa = kf;
    const qb = -(kf * (a0 + b0) + kr);
    const qc = kf * a0 * b0 - kr * c0;

    const discriminant = Math.max(qb * qb - 4 * qa * qc, 0);
    const root = Math.sqrt(discriminant);

    const lower = -c0;
    const upper = Math.min(a0, b0);

    const r1 = (-qb + root) / (2 * qa);
    const r2 = (-qb - root) / (2 * qa);

    const inRange = (x) => x >= lower - EPS && x <= upper + EPS;

    if (inRange(r1)) {
      complexAtEquilibrium = r1;
    } else if (inRange(r2)) {
      complexAtEquilibrium = r2;
    } else {
      complexAtEquilibrium = Math.min(Math.max(r2, lower), upper);
    }
  }

  // Concentrations libres à l'équilibre
  const freeA = Math.max(a0 - complexAtEquilibrium, 0);
  const freeB = Math.max(b0 - complexAtEquilibrium, 0);

  // Relaxation vers l'équilibre
  const forwardRate = kf * (freeA + freeB);
  const k = forwardRate + kr;

  if (k <= EPS) {
    return;
  }

  const nextComplex =
    complexAtEquilibrium + (c0 - complexAtEquilibrium) * Math.exp(-k * dt);

  const shift = clampNonNegative(nextComplex) - c0;

  s[c] = clampNonNegative(c0 + shift);
  s[a] = clampNonNegative(a0 - shift);
  s[b] = clampNonNegative(b0 - shift);
}

/**
 * X -> Y.
 * Décroissance exponentielle exacte.
 * Conserve X + Y.
 */
function transfer(s, from, to, k, dt) {
  if (k <= EPS || s[from] <= EPS) {
    return;
  }

  const moved = s[from] * (1 - Math.exp(-k * dt));

  s[from] = clampNonNegative(s[from] - moved);
  s[to] = clampNonNegative(s[to] + moved);
}

/**
 * Un pas de simulation, avec sous-pas internes si dt est trop grand.
 */
export function step(s, p, ca, dt) {
  const substeps = Math.max(1, Math.ceil(dt / MAX_INTERNAL_DT));
  const h = dt / substeps;

  for (let i = 0; i < substeps; i += 1) {
    stepOnce(s, p, ca, h);
  }
}

function stepOnce(s, p, ca, dt) {
  const c2 = ca * ca;

  // ---------------------------------------------------------------
  // M1 — Calmoduline
  // M0 <-> M2C <-> M4 <-> M2N <-> M0
  // ---------------------------------------------------------------
  relax(s, "M0", "M2C", p.calmodulin.kon2 * c2, p.calmodulin.koff2, dt);
  relax(s, "M2C", "M4", p.calmodulin.kon1 * c2, p.calmodulin.koff1, dt);
  relax(s, "M4", "M2N", p.calmodulin.koff2, p.calmodulin.kon2 * c2, dt);
  relax(s, "M2N", "M0", p.calmodulin.koff1, p.calmodulin.kon1 * c2, dt);
  relax(s, "M2C", "M2N", K_EXCHANGE, K_EXCHANGE, dt);

  // ---------------------------------------------------------------
  // M3 — Inhibition de PP1 par I-1p
  // I1 <-> I1p
  // I1p + PP1 <-> PP1I1p
  // ---------------------------------------------------------------
  relax(
    s,
    "I1",
    "I1p",
    K_I1_PHOS * p.pkaActivity,
    K_I1_DEPHOS,
    dt,
  );

  bind(
    s,
    "I1p",
    "PP1",
    "PP1I1p",
    p.pp1i1.kon,
    p.pp1i1.koff,
    dt,
  );

  // ---------------------------------------------------------------
  // M2 — CaMKII
  // K0 + CaM active <-> KC
  // KC -> KpC
  // KpC -> Kp
  // ---------------------------------------------------------------
  const camActive = s.M2C + s.M4;

  relax(
    s,
    "K0",
    "KC",
    p.camkii.konBind * camActive,
    p.camkii.koffBind,
    dt,
  );

  transfer(s, "KC", "KpC", p.camkii.kPhosphorylation, dt);
  transfer(s, "KpC", "Kp", p.camkii.koffBind, dt);

  // ---------------------------------------------------------------
  // Déphosphorylation explicite
  // Kp + PP1 <-> PP1KpC
  // PP1KpC -> D + PP1
  // ---------------------------------------------------------------
  bind(
    s,
    "Kp",
    "PP1",
    "PP1KpC",
    p.dephosphorylation.kon,
    p.dephosphorylation.koff,
    dt,
  );

  if (s.PP1KpC > EPS) {
    const dephosphorylated =
      s.PP1KpC * (1 - Math.exp(-p.dephosphorylation.kcat * dt));

    s.PP1KpC = clampNonNegative(s.PP1KpC - dephosphorylated);
    s.PP1 = clampNonNegative(s.PP1 + dephosphorylated);
    s.D = clampNonNegative(s.D + dephosphorylated);
  }

  // ---------------------------------------------------------------
  // Recovery
  // D -> K0
  // ---------------------------------------------------------------
  transfer(s, "D", "K0", K_RECOVERY, dt);
}

export function invariantsOf(s) {
  return {
    cam: s.M0 + s.M2N + s.M2C + s.M4,
    camkii: s.K0 + s.KC + s.D + s.KpC + s.Kp + s.PP1KpC,
    pp1: s.PP1 + s.PP1I1p + s.PP1KpC,
    i1: s.I1 + s.I1p + s.PP1I1p,
  };
}

export function simulate({
  params = DEFAULT_PARAMS,
  protocol = DEFAULT_PROTOCOL,
  duration,
  dt = 0.01,
  maxPoints = 520,
}) {
  if (!Number.isFinite(duration) || duration <= 0) {
    throw new Error("duration doit être un nombre positif.");
  }

  if (!Number.isFinite(dt) || dt <= 0) {
    throw new Error("dt doit être un nombre positif.");
  }

  const s = { ...INITIAL_STATE };
  const steps = Math.max(1, Math.round(duration / dt));
  const stride = Math.max(1, Math.ceil(steps / maxPoints));
  const rows = [];

  const push = (t) => {
    rows.push({
      t: Math.round(t * 1000) / 1000,
      ca: calciumAt(t, protocol),
      ...s,
    });
  };

  push(0);

  for (let i = 1; i <= steps; i += 1) {
    const tPrevious = (i - 1) * dt;
    const ca = calciumAt(tPrevious, protocol);

    step(s, params, ca, dt);

    if (i % stride === 0) {
      push(i * dt);
    }
  }

  const last = rows[rows.length - 1];

  if (last.t < duration - 1e-9) {
    push(duration);
  }

  return {
    rows,
    duration,
    invariants: invariantsOf(s),
  };
}
