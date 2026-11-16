// omegaModel.js — Moteur Omega V1 du simulateur
// Unités : Concentrations (µM), Temps (s), kon (µM⁻¹·s⁻¹), koff/kcat (s⁻¹)

export const K_EXCHANGE = 1.0;
export const K_I1_PHOS = 2.0;
export const K_I1_DEPHOS = 0.15;
export const K_RECOVERY = 0.02;

export const TOTALS = { cam: 10, camkii: 50, pp1: 2, i1: 1 };

export const DEFAULT_PARAMS = {
  calmodulin: { kon1: 5.0, koff1: 50.0, kon2: 100.0, koff2: 2000.0 },
  camkii: { konBind: 1.27e-3, koffBind: 0.07, kPhosphorylation: 0.5 },
  dephosphorylation: { kon: 1.0, koff: 0.001, kcat: 0.41 },
  pp1i1: { kon: 1.0, koff: 0.001 },
  pkaActivity: 0.05,
};

export const INITIAL_STATE = {
  M0: 10.0, M2N: 0.0, M2C: 0.0, M4: 0.0,
  K0: 50.0, KC: 0.0, KpC: 0.0, Kp: 0.0, D: 0.0,
  PP1: 2.0, I1: 1.0, I1p: 0.0, PP1I1p: 0.0, PP1KpC: 0.0,
};

export const DEFAULT_PROTOCOL = {
  mode: "train", amplitude: 20, basal: 0.1, width: 2, period: 20, delay: 5,
};

const EPS = 1e-12;
const MAX_INTERNAL_DT = 0.005;
const clampNonNegative = (x) => (x > 0 ? x : 0);

export function calciumAt(t, p) {
  if (p.mode === "constant") return Math.max(p.amplitude, 0);
  if (t < p.delay) return Math.max(p.basal, 0);
  if (p.mode === "impulse") return t < p.delay + p.width ? Math.max(p.amplitude, 0) : Math.max(p.basal, 0);
  const cycle = Math.max(p.period, p.width + 0.01);
  const phase = (t - p.delay) % cycle;
  return phase < p.width ? Math.max(p.amplitude, 0) : Math.max(p.basal, 0);
}

function relax(s, a, b, ka, kb, dt) {
  const total = s[a] + s[b];
  if (total <= EPS) { s[a] = 0; s[b] = 0; return; }
  const k = ka + kb;
  if (k <= EPS) return;
  const equilibriumB = (ka / k) * total;
  const nextB = equilibriumB + (s[b] - equilibriumB) * Math.exp(-k * dt);
  s[b] = clampNonNegative(nextB);
  s[a] = clampNonNegative(total - s[b]);
}

function bind(s, a, b, c, kf, kr, dt) {
  const a0 = s[a], b0 = s[b], c0 = s[c];
  const total = a0 + b0 + c0;
  if (total <= EPS || (kf <= EPS && kr <= EPS)) return;

  let deltaEq = 0; // Variation du complexe (Delta)

  if (kf > EPS) {
    const qa = kf;
    const qb = -(kf * (a0 + b0) + kr);
    const qc = kf * a0 * b0 - kr * c0;
    const discriminant = Math.max(qb * qb - 4 * qa * qc, 0);
    const root = Math.sqrt(discriminant);
    
    // On cherche la racine physiquement valide (0 <= Delta <= min(a0, b0))
    const r1 = (-qb + root) / (2 * qa);
    const r2 = (-qb - root) / (2 * qa);
    const upper = Math.min(a0, b0);
    
    if (r1 >= -EPS && r1 <= upper + EPS) deltaEq = r1;
    else if (r2 >= -EPS && r2 <= upper + EPS) deltaEq = r2;
    else deltaEq = Math.min(Math.max(r2, 0), upper); // Fallback sécurisé
  } else if (kr > EPS) {
    // kf est nul : le complexe ne se forme pas, il ne peut que se dissocier
    deltaEq = -c0; 
  }

  const cEq = c0 + deltaEq; // Concentration absolue à l'équilibre
  const freeA = Math.max(a0 - deltaEq, 0);
  const freeB = Math.max(b0 - deltaEq, 0);

  const forwardRate = kf * (freeA + freeB); // Approximation de la relaxation
  const k = forwardRate + kr;
  if (k <= EPS) return;

  const nextComplex = cEq + (c0 - cEq) * Math.exp(-k * dt);
  const shift = clampNonNegative(nextComplex) - c0;

  s[c] = clampNonNegative(c0 + shift);
  s[a] = clampNonNegative(a0 - shift);
  s[b] = clampNonNegative(b0 - shift);
}

function transfer(s, from, to, k, dt) {
  if (k <= EPS || s[from] <= EPS) return;
  const moved = s[from] * (1 - Math.exp(-k * dt));
  s[from] = clampNonNegative(s[from] - moved);
  s[to] = clampNonNegative(s[to] + moved);
}

export function step(s, p, ca, dt) {
  if (dt <= 0) return;
  const substeps = Math.max(1, Math.ceil(dt / MAX_INTERNAL_DT));
  const h = dt / substeps;
  for (let i = 0; i < substeps; i += 1) stepOnce(s, p, ca, h);
}

// Approximation par découpage des réactions (Operator Splitting).
// Chaque réaction est exacte, mais leur enchaînement ne l'est pas globalement.
function stepOnce(s, p, ca, dt) {
  const c2 = ca * ca;

  // M1
  relax(s, "M0", "M2C", p.calmodulin.kon2 * c2, p.calmodulin.koff2, dt);
  relax(s, "M2C", "M4", p.calmodulin.kon1 * c2, p.calmodulin.koff1, dt);
  relax(s, "M4", "M2N", p.calmodulin.koff2, p.calmodulin.kon2 * c2, dt);
  relax(s, "M2N", "M0", p.calmodulin.koff1, p.calmodulin.kon1 * c2, dt);
  relax(s, "M2C", "M2N", K_EXCHANGE, K_EXCHANGE, dt);

  // M3
  relax(s, "I1", "I1p", K_I1_PHOS * p.pkaActivity, K_I1_DEPHOS, dt);
  bind(s, "I1p", "PP1", "PP1I1p", p.pp1i1.kon, p.pp1i1.koff, dt);

  // M2
  const camActive = s.M2C + s.M4;
  relax(s, "K0", "KC", p.camkii.konBind * camActive, p.camkii.koffBind, dt);
  transfer(s, "KC", "KpC", p.camkii.kPhosphorylation, dt);
  transfer(s, "KpC", "Kp", p.camkii.koffBind, dt);

  // Déphosphorylation
  bind(s, "Kp", "PP1", "PP1KpC", p.dephosphorylation.kon, p.dephosphorylation.koff, dt);
  if (s.PP1KpC > EPS) {
    const dephosphorylated = s.PP1KpC * (1 - Math.exp(-p.dephosphorylation.kcat * dt));
    s.PP1KpC = clampNonNegative(s.PP1KpC - dephosphorylated);
    s.PP1 = clampNonNegative(s.PP1 + dephosphorylated);
    s.D = clampNonNegative(s.D + dephosphorylated);
  }

  // Recovery
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
  if (!Number.isFinite(duration) || duration <= 0) throw new Error("duration > 0 requise.");
  if (!Number.isFinite(dt) || dt <= 0) throw new Error("dt > 0 requis.");
  if (!Number.isInteger(maxPoints) || maxPoints <= 0) throw new Error("maxPoints doit être un entier > 0.");

  const s = { ...INITIAL_STATE };
  const rows = [];
  let t = 0;
  let stepsTaken = 0;
  
  const push = (time) => {
    rows.push({ t: Math.round(time * 1000) / 1000, ca: calciumAt(time, protocol), ...s });
  };

  push(0);

  const approxSteps = Math.ceil(duration / dt);
  const stride = Math.max(1, Math.ceil(approxSteps / maxPoints));

  while (t < duration - 1e-9) {
    let h = Math.min(dt, duration - t); // Intégration exacte jusqu'à duration
    const ca = calciumAt(t, protocol);
    
    step(s, params, ca, h);
    
    t += h;
    stepsTaken += 1;

    if (stepsTaken % stride === 0 || t >= duration - 1e-9) {
      push(t);
    }
  }

  return { rows, duration, invariants: invariantsOf(s) };
  }
