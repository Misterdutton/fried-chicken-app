/*
 * RC Frame Engine — 2D analysis & design of reinforced-concrete moment frames (<= 8 storeys)
 * Units: kN, m (analysis) / N, mm, MPa (section design)
 *
 * - Direct stiffness method, 3 DOF per node, fixed base
 * - Load cases: D (dead), L (live), E (equivalent static seismic, DPT 1301/1302-61 style)
 * - Strength combos per ACI 318-19 / EIT
 * - Beam flexure & shear, column P-M interaction & shear
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.RCFrame = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const G = 9.81;
  const ES = 200000; // MPa
  const MAX_STORIES = 8;
  const MAX_BAYS = 10;

  const SYSTEMS = {
    OMF: { R: 3, Cd: 2.5, label: 'โครงต้านแรงดัดธรรมดา (OMF)' },
    IMF: { R: 5, Cd: 4.5, label: 'โครงต้านแรงดัดความเหนียวปานกลาง (IMF)' },
    SMF: { R: 8, Cd: 5.5, label: 'โครงต้านแรงดัดความเหนียวพิเศษ (SMF)' },
  };

  // Strength combos (ACI 318-19 / EIT) with vertical seismic Ev = 0.2·SDS·D (DPT 1301/1302-61)
  // Wind combos added when wind is on: 1.2D+1.0L±γw·W, 0.9D±γw·W (γw is a user input)
  function makeCombos(SDS, windFactor = 0) {
    const ev = 0.2 * SDS;
    const r = (v) => Math.round(v * 1000) / 1000;
    const gw = windFactor;
    const wind = gw > 0 ? [
      { name: `1.2D+1.0L+${r(gw)}W`, D: 1.2, L: 1.0, E: 0, W: gw },
      { name: `1.2D+1.0L-${r(gw)}W`, D: 1.2, L: 1.0, E: 0, W: -gw },
      { name: `0.9D+${r(gw)}W`, D: 0.9, L: 0, E: 0, W: gw },
      { name: `0.9D-${r(gw)}W`, D: 0.9, L: 0, E: 0, W: -gw },
    ] : [];
    return [
      { name: '1.4D', D: 1.4, L: 0, E: 0 },
      { name: '1.2D+1.6L', D: 1.2, L: 1.6, E: 0 },
      { name: `${r(1.2 + ev)}D+1.0L+1.0E`, D: 1.2 + ev, L: 1.0, E: 1.0 },
      { name: `${r(1.2 + ev)}D+1.0L-1.0E`, D: 1.2 + ev, L: 1.0, E: -1.0 },
      { name: `${r(0.9 - ev)}D+1.0E`, D: 0.9 - ev, L: 0, E: 1.0 },
      { name: `${r(0.9 - ev)}D-1.0E`, D: 0.9 - ev, L: 0, E: -1.0 },
    ].concat(wind);
  }

  // linear interpolation in a user-supplied spectrum table [[T, Sa], ...]
  function interpSpectrum(table, T) {
    const t = table.filter((r) => isFinite(r[0]) && isFinite(r[1])).sort((a, b) => a[0] - b[0]);
    if (t.length < 2) throw new Error('ตารางสเปกตรัมต้องมีอย่างน้อย 2 แถว (คาบ T, S\u2090)');
    if (T <= t[0][0]) return t[0][1];
    for (let i = 0; i < t.length - 1; i++) {
      if (T <= t[i + 1][0]) return t[i][1] + (t[i + 1][1] - t[i][1]) * (T - t[i][0]) / (t[i + 1][0] - t[i][0]);
    }
    return t[t.length - 1][1];
  }

  // ---------------------------------------------------------------- wind (DPT 1311-50, simplified procedure)
  // p = Iw·q·Ce·Cg·Cp,  q = ½ρV²,  V = TF·V50 (strength) or V50 (serviceability)
  function windLoads(model) {
    const { inp, ns, ys, xs } = model;
    const w = inp.wind;
    if (!w || !w.on) return null;
    ['V50', 'TF', 'IwU', 'IwS', 'factor'].forEach((k) => {
      if (!(w[k] > 0)) throw new Error('ค่าแรงลมต้องมากกว่า 0 (' + k + ')');
    });
    const rho = 1.25; // kg/m3
    const H = ys[ns], D = xs[xs.length - 1];
    const ceA = (z) => Math.max(0.9, Math.pow(z / 10, 0.2));
    const ceB = (z) => Math.max(0.7, 0.7 * Math.pow(z / 12, 0.3));
    const Ce = w.terrain === 'B' ? ceB : ceA;
    const Cg = w.Cg ?? 2.0;
    const CpW = w.CpW ?? 0.8;
    const hd = H / D;
    const CpL = w.CpL ?? (hd < 0.25 ? -0.3 : hd >= 1 ? -0.5 : -0.3 - 0.2 * (hd - 0.25) / 0.75);
    const q = 0.5 * rho * (w.TF * w.V50) ** 2 / 1000; // kPa, strength
    const CeL = Ce(H / 2);
    const pL = w.IwU * q * CeL * Cg * Math.abs(CpL);
    const F = new Array(ns + 1).fill(0), levels = [];
    for (let j = 1; j <= ns; j++) {
      const zTop = j === ns ? ys[j] : (ys[j] + ys[j + 1]) / 2;
      const zBot = (ys[j] + ys[j - 1]) / 2;
      const pW = w.IwU * q * Ce(ys[j]) * Cg * CpW;
      F[j] = (pW + pL) * (zTop - zBot) * inp.trib;
      levels.push({ level: j, z: ys[j], Ce: Ce(ys[j]), pW, pL, F: F[j] });
    }
    const svc = (w.IwS / w.IwU) / (w.TF * w.TF); // serviceability / strength pressure ratio
    return { q, CpW, CpL, Cg, CeL, F, levels, V: F.reduce((a, b) => a + b, 0), svc, H, D, hd, factor: w.factor };
  }
  const COMBOS = makeCombos(0);

  // Seismic design category (ก ข ค ง = A B C D) from SDS, SD1 and risk category
  function designCategory(SDS, SD1, Ie) {
    const IV = Ie >= 1.5;
    const bySDS = SDS < 0.167 ? 0 : SDS < 0.33 ? (IV ? 2 : 1) : SDS < 0.5 ? (IV ? 3 : 2) : 3;
    const bySD1 = SD1 < 0.067 ? 0 : SD1 < 0.133 ? (IV ? 2 : 1) : SD1 < 0.20 ? (IV ? 3 : 2) : 3;
    return Math.max(bySDS, bySD1);
  }
  const SDC_NAME = ['ก', 'ข', 'ค', 'ง'];
  const SYSTEM_MAX_SDC = { OMF: 1, IMF: 2, SMF: 3 }; // highest category each system is permitted in

  const barArea = (db) => Math.PI * db * db / 4; // mm2

  function beta1(fc) {
    if (fc <= 28) return 0.85;
    if (fc >= 55) return 0.65;
    return 0.85 - 0.05 * (fc - 28) / 7;
  }

  // ---------------------------------------------------------------- linear algebra
  function luFactor(A, n) {
    const piv = new Int32Array(n);
    for (let k = 0; k < n; k++) {
      let p = k, max = Math.abs(A[k][k]);
      for (let i = k + 1; i < n; i++) {
        const v = Math.abs(A[i][k]);
        if (v > max) { max = v; p = i; }
      }
      if (max < 1e-12) throw new Error('โครงสร้างไม่เสถียร (stiffness matrix singular)');
      piv[k] = p;
      if (p !== k) { const t = A[k]; A[k] = A[p]; A[p] = t; }
      const rk = A[k], akk = rk[k];
      for (let i = k + 1; i < n; i++) {
        const ri = A[i];
        const f = ri[k] / akk;
        if (f === 0) continue;
        ri[k] = f;
        for (let j = k + 1; j < n; j++) ri[j] -= f * rk[j];
      }
    }
    return piv;
  }

  function luSolve(A, piv, b, n) {
    const x = Float64Array.from(b);
    for (let k = 0; k < n; k++) {
      const p = piv[k];
      if (p !== k) { const t = x[k]; x[k] = x[p]; x[p] = t; }
      for (let i = k + 1; i < n; i++) x[i] -= A[i][k] * x[k];
    }
    for (let i = n - 1; i >= 0; i--) {
      let s = x[i];
      for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j];
      x[i] = s / A[i][i];
    }
    return x;
  }

  // ---------------------------------------------------------------- element helpers
  function localStiffness(E, A, I, L) {
    const a = E * A / L, b = 12 * E * I / L ** 3, c = 6 * E * I / L ** 2, d = 4 * E * I / L, e = 2 * E * I / L;
    return [
      [a, 0, 0, -a, 0, 0],
      [0, b, c, 0, -b, c],
      [0, c, d, 0, -c, e],
      [-a, 0, 0, a, 0, 0],
      [0, -b, -c, 0, b, -c],
      [0, c, e, 0, -c, d],
    ];
  }

  // transform global -> local : uL = T uG
  function toLocal(c, s, v) {
    return [
      c * v[0] + s * v[1], -s * v[0] + c * v[1], v[2],
      c * v[3] + s * v[4], -s * v[3] + c * v[4], v[5],
    ];
  }
  function toGlobal(c, s, v) {
    return [
      c * v[0] - s * v[1], s * v[0] + c * v[1], v[2],
      c * v[3] - s * v[4], s * v[3] + c * v[4], v[5],
    ];
  }

  // fixed-end forces (local) for uniform load w acting in local -y (gravity on a beam)
  function fixedEnd(w, L) {
    return [0, w * L / 2, w * L * L / 12, 0, w * L / 2, -w * L * L / 12];
  }

  // ---------------------------------------------------------------- model
  function buildModel(inp) {
    const spans = inp.spans, st = inp.stories;
    if (!spans.length || spans.length > MAX_BAYS) throw new Error(`จำนวนช่วงคานต้องอยู่ระหว่าง 1–${MAX_BAYS}`);
    if (!st.length || st.length > MAX_STORIES) throw new Error(`จำนวนชั้นต้องอยู่ระหว่าง 1–${MAX_STORIES} ชั้น`);
    spans.forEach((L) => { if (!(L > 0)) throw new Error('ความยาวช่วงคานต้องมากกว่า 0'); });
    st.forEach((s, k) => {
      ['h', 'colB', 'colH', 'beamB', 'beamH'].forEach((f) => {
        if (!(s[f] > 0)) throw new Error(`ข้อมูลชั้น ${k + 1} ไม่ถูกต้อง (${f})`);
      });
    });

    const nx = spans.length + 1, ns = st.length;
    const xs = [0]; spans.forEach((L) => xs.push(xs[xs.length - 1] + L));
    const ys = [0]; st.forEach((s) => ys.push(ys[ys.length - 1] + s.h));
    const E = 4700 * Math.sqrt(inp.fc) * 1000; // kN/m2
    const gc = inp.gammaC ?? 24;
    const crB = inp.crackBeam ?? 0.35, crC = inp.crackCol ?? 0.70;

    const nodes = [];
    for (let j = 0; j <= ns; j++)
      for (let i = 0; i < nx; i++) nodes.push({ id: nodes.length, i, j, x: xs[i], y: ys[j], fixed: j === 0 });
    const nid = (i, j) => j * nx + i;

    const elems = [];
    for (let j = 0; j < ns; j++) {
      const s = st[j];
      const b = s.colB / 1000, h = s.colH / 1000;
      for (let i = 0; i < nx; i++) {
        elems.push({
          id: elems.length, type: 'col', name: `C${j + 1}-${i + 1}`, story: j + 1, line: i,
          n1: nid(i, j), n2: nid(i, j + 1), b, h, A: b * h, I: crC * b * h ** 3 / 12, E,
          L: s.h, c: 0, s: 1, wD: 0, wL: 0, selfW: gc * b * h * s.h,
          bars: { n: s.colBarN, db: s.colBarDb },
        });
      }
    }
    for (let j = 1; j <= ns; j++) {
      const s = st[j - 1], roof = j === ns;
      const b = s.beamB / 1000, h = s.beamH / 1000;
      const dl = roof ? inp.roofDL : inp.floorDL, ll = roof ? inp.roofLL : inp.floorLL;
      const wall = roof ? (inp.parapet ?? 0) : (inp.wall ?? 0);
      for (let i = 0; i < nx - 1; i++) {
        elems.push({
          id: elems.length, type: 'beam', name: `B${j}-${i + 1}`, story: j, line: i,
          n1: nid(i, j), n2: nid(i + 1, j), b, h, A: b * h, I: crB * b * h ** 3 / 12, E,
          L: spans[i], c: 1, s: 0,
          wD: dl * inp.trib + gc * b * h + wall, wL: ll * inp.trib,
          barDb: s.beamBarDb,
        });
      }
    }
    return { inp, nx, ns, xs, ys, nodes, elems, ndof: nodes.length * 3 };
  }

  // ---------------------------------------------------------------- seismic loads
  function seismic(model) {
    const { inp, ns, nx, ys, elems } = model;
    const sys = SYSTEMS[inp.system] || SYSTEMS.IMF;
    const R = inp.R ?? sys.R, Cd = inp.Cd ?? sys.Cd, Ie = inp.Ie ?? 1;
    const liveFrac = inp.liveFrac ?? 0;
    const WD = new Array(ns + 1).fill(0), WL = new Array(ns + 1).fill(0);
    elems.forEach((e) => {
      if (e.type === 'beam') { WD[e.story] += e.wD * e.L; WL[e.story] += e.wL * e.L; }
      else { WD[e.story] += e.selfW / 2; if (e.story - 1 > 0) WD[e.story - 1] += e.selfW / 2; else WD[0] += e.selfW / 2; }
    });
    const W = WD.map((d, k) => d + liveFrac * WL[k]);
    const Wtot = W.slice(1).reduce((a, b) => a + b, 0);
    if (!(Wtot > 0)) throw new Error('น้ำหนักอาคารต้องมากกว่า 0');
    const H = ys[ns];
    const Ta = 0.02 * H; // RC moment frame, DPT 1301/1302-61
    const SDS = inp.SDS, SD1 = inp.SD1;
    const Ts = SD1 / SDS;
    const T = Ta;
    const table = inp.spectrum === 'table';
    const Sa = table ? interpSpectrum(inp.saTable || [], T) : (T <= Ts ? SDS : SD1 / T);
    let Cs = Sa * Ie / R;
    const CsMin = Math.max(0.044 * SDS * Ie, 0.01);
    Cs = Math.max(Cs, CsMin);
    const V = Cs * Wtot;
    const k = T <= 0.5 ? 1 : T >= 2.5 ? 2 : 1 + (T - 0.5) / 2;
    let sum = 0;
    for (let j = 1; j <= ns; j++) sum += W[j] * ys[j] ** k;
    const F = new Array(ns + 1).fill(0);
    for (let j = 1; j <= ns; j++) F[j] = V * W[j] * ys[j] ** k / sum;
    const sdc = designCategory(SDS, SD1, Ie);
    const sysKey = SYSTEMS[inp.system] ? inp.system : 'IMF';
    const systemOk = sdc <= SYSTEM_MAX_SDC[sysKey];
    const need = sdc <= 1 ? 'OMF' : sdc === 2 ? 'IMF' : 'SMF';
    return { R, Cd, Ie, H, Ta, T, Ts, Sa, Cs, CsMin, V, k, W, WD, WL, Wtot, F, system: sys.label, nx,
      sdc, sdcName: SDC_NAME[sdc], systemOk, minSystem: need, spectrum: table ? 'table' : 'map' };
  }

  // ---------------------------------------------------------------- analysis
  function analyze(inp) {
    const model = buildModel(inp);
    const { nodes, elems, ndof, nx, ns } = model;
    const eq = seismic(model);
    const wind = windLoads(model);

    // dof map
    const map = new Int32Array(ndof).fill(-1);
    let nf = 0;
    nodes.forEach((n) => { if (!n.fixed) for (let d = 0; d < 3; d++) map[n.id * 3 + d] = nf++; });

    const K = Array.from({ length: nf }, () => new Float64Array(nf));
    elems.forEach((e) => {
      const kl = localStiffness(e.E, e.A, e.I, e.L);
      e.kl = kl;
      const c = e.c, s = e.s;
      // kg = T^T kl T
      const T = [
        [c, s, 0, 0, 0, 0], [-s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
        [0, 0, 0, c, s, 0], [0, 0, 0, -s, c, 0], [0, 0, 0, 0, 0, 1],
      ];
      const kT = kl.map((row) => T[0].map((_, j) => row.reduce((acc, v, m) => acc + v * T[m][j], 0)));
      const kg = T[0].map((_, i) => T[0].map((__, j) => T.reduce((acc, row, m) => acc + row[i] * kT[m][j], 0)));
      const dofs = [e.n1 * 3, e.n1 * 3 + 1, e.n1 * 3 + 2, e.n2 * 3, e.n2 * 3 + 1, e.n2 * 3 + 2];
      e.dofs = dofs;
      for (let a = 0; a < 6; a++) {
        const ra = map[dofs[a]]; if (ra < 0) continue;
        for (let b = 0; b < 6; b++) {
          const cb = map[dofs[b]]; if (cb < 0) continue;
          K[ra][cb] += kg[a][b];
        }
      }
    });
    const piv = luFactor(K, nf);

    function runCase(name) {
      const P = new Float64Array(ndof);
      const w = new Float64Array(elems.length);
      elems.forEach((e) => {
        if (e.type === 'beam') w[e.id] = name === 'D' ? e.wD : name === 'L' ? e.wL : 0;
        else if (name === 'D') { P[e.n1 * 3 + 1] -= e.selfW / 2; P[e.n2 * 3 + 1] -= e.selfW / 2; }
      });
      if (name === 'E' || name === 'W') {
        const Fl = name === 'E' ? eq.F : wind.F;
        nodes.forEach((n) => { if (n.j > 0) P[n.id * 3] += Fl[n.j] / nx; });
      }
      // equivalent nodal loads from member loads
      elems.forEach((e) => {
        if (!w[e.id]) return;
        const fg = toGlobal(e.c, e.s, fixedEnd(w[e.id], e.L));
        e.dofs.forEach((d, a) => { P[d] -= fg[a]; });
      });
      const rhs = new Float64Array(nf);
      for (let d = 0; d < ndof; d++) if (map[d] >= 0) rhs[map[d]] = P[d];
      const u = luSolve(K, piv, rhs, nf);
      const U = new Float64Array(ndof);
      for (let d = 0; d < ndof; d++) if (map[d] >= 0) U[d] = u[map[d]];
      const F = elems.map((e) => {
        const ul = toLocal(e.c, e.s, e.dofs.map((d) => U[d]));
        const fe = fixedEnd(w[e.id], e.L);
        return e.kl.map((row, i) => row.reduce((acc, v, j) => acc + v * ul[j], 0) + fe[i]);
      });
      // reactions
      const Rx = new Float64Array(ndof);
      elems.forEach((e, k) => {
        const fg = toGlobal(e.c, e.s, F[k]);
        e.dofs.forEach((d, a) => { Rx[d] += fg[a]; });
      });
      for (let d = 0; d < ndof; d++) Rx[d] -= P[d];
      return { name, U, F, w, R: Rx };
    }

    const cases = { D: runCase('D'), L: runCase('L'), E: runCase('E') };
    if (wind) cases.W = runCase('W');

    function combine(cf) {
      const U = new Float64Array(ndof), R = new Float64Array(ndof);
      const w = new Float64Array(elems.length);
      const F = elems.map(() => [0, 0, 0, 0, 0, 0]);
      ['D', 'L', 'E', 'W'].forEach((k) => {
        const f = cf[k]; if (!f) return;
        const c = cases[k];
        for (let d = 0; d < ndof; d++) { U[d] += f * c.U[d]; R[d] += f * c.R[d]; }
        elems.forEach((e, m) => { w[m] += f * c.w[m]; for (let i = 0; i < 6; i++) F[m][i] += f * c.F[m][i]; });
      });
      return { name: cf.name, U, F, w, R };
    }
    const combos = makeCombos(inp.SDS, wind ? wind.factor : 0).map(combine);
    const service = combine({ name: 'D+L', D: 1, L: 1 });

    const stories = storyChecks(model, eq, cases.E);
    const design = designMembers(model, combos);
    if (wind) windDrift(model, wind, cases.W);

    return { model, eq, wind, cases, combos, service, stories, design };
  }

  // internal forces along a member (local), x from node 1
  function memberDiagram(f, w, L, n = 20) {
    const pts = [];
    for (let k = 0; k <= n; k++) {
      const x = L * k / n;
      pts.push({ x, N: -f[0], V: f[1] - w * x, M: -f[2] + f[1] * x - w * x * x / 2 });
    }
    // include point of max moment for beams with load
    if (w > 0) {
      const x0 = f[1] / w;
      if (x0 > 0 && x0 < L) pts.push({ x: x0, N: -f[0], V: 0, M: -f[2] + f[1] * x0 - w * x0 * x0 / 2 });
      pts.sort((a, b) => a.x - b.x);
    }
    return pts;
  }

  // serviceability drift under wind (V50, Iw serviceability): recommended limit h/500
  function windDrift(model, wind, Wcase) {
    const { ns, nx, ys, nodes } = model;
    const d = new Array(ns + 1).fill(0);
    for (let j = 1; j <= ns; j++) {
      let s = 0; nodes.forEach((n) => { if (n.j === j) s += Wcase.U[n.id * 3]; });
      d[j] = s / nx * wind.svc;
    }
    let Vx = 0;
    for (let j = ns; j >= 1; j--) {
      Vx += wind.F[j];
      const L = wind.levels[j - 1];
      L.Vx = Vx; L.dSvc = d[j]; L.drift = d[j] - d[j - 1]; L.ratio = L.drift / (ys[j] - ys[j - 1]);
      L.ok = L.ratio <= 1 / 500;
    }
    wind.topRatio = d[ns] / ys[ns];
  }

  // ---------------------------------------------------------------- story drift
  function storyChecks(model, eq, Ecase) {
    const { ns, nx, ys, nodes } = model;
    const out = [];
    const dxe = new Array(ns + 1).fill(0);
    for (let j = 1; j <= ns; j++) {
      let s = 0, m = 0;
      nodes.forEach((n) => { if (n.j === j) { s += Ecase.U[n.id * 3]; m = Math.max(m, Math.abs(Ecase.U[n.id * 3])); } });
      dxe[j] = s / nx;
    }
    // Rayleigh period (information)
    let num = 0, den = 0;
    for (let j = 1; j <= ns; j++) { num += eq.W[j] * dxe[j] ** 2; den += eq.F[j] * dxe[j]; }
    eq.Tray = 2 * Math.PI * Math.sqrt(num / (G * den));

    const limit = eq.Ie >= 1.5 ? 0.010 : eq.Ie >= 1.25 ? 0.015 : 0.020;
    for (let j = 1; j <= ns; j++) {
      const hs = ys[j] - ys[j - 1];
      const dx = eq.Cd * dxe[j] / eq.Ie;
      const dxPrev = eq.Cd * dxe[j - 1] / eq.Ie;
      const drift = dx - dxPrev;
      let Px = 0, Vx = 0;
      for (let k = j; k <= ns; k++) { Px += eq.WD[k] + eq.WL[k]; Vx += eq.F[k]; }
      const theta = Px * (drift) * eq.Ie / (Vx * hs * eq.Cd);
      const thetaMax = Math.min(0.5 / eq.Cd, 0.25);
      out.push({
        level: j, h: hs, elev: ys[j], W: eq.W[j], F: eq.F[j], Vx, dxe: dxe[j], dx, drift,
        ratio: drift / hs, limit, okDrift: drift / hs <= limit, theta, thetaMax,
        okTheta: theta <= thetaMax, pdelta: theta > 0.10,
      });
    }
    return out;
  }

  // ---------------------------------------------------------------- section design
  function beamFlexure(Mu, b, d, fc, fy) {
    // Mu kN·m (>=0), b,d mm
    const phi = 0.9;
    const Rn = Mu * 1e6 / (phi * b * d * d);
    const q = 2 * Rn / (0.85 * fc);
    const rhoMin = Math.max(0.25 * Math.sqrt(fc) / fy, 1.4 / fy);
    const rhoTC = 0.85 * beta1(fc) * fc / fy * 0.003 / (0.003 + 0.005);
    if (q >= 1) return { ok: false, rho: NaN, As: NaN, rhoMin, rhoMax: rhoTC, msg: 'หน้าตัดเล็กเกินไป' };
    const rho = 0.85 * fc / fy * (1 - Math.sqrt(1 - q));
    const use = Mu > 0 ? Math.max(rho, rhoMin) : rhoMin;
    const ok = rho <= rhoTC;
    return { ok, rho, As: use * b * d, rhoMin, rhoMax: rhoTC, msg: ok ? '' : 'ρ เกินขีดจำกัด (ต้องขยายหน้าตัด/ใช้เหล็กรับแรงอัด)' };
  }

  function shearDesign(Vu, b, d, fc, fyt, ds, legs = 2, Nu = 0, Ag = 0) {
    // Vu kN, Nu kN (compression +)
    const phi = 0.75;
    const axial = Ag > 0 ? Math.max(0, 1 + Nu * 1000 / (14 * Ag)) : 1;
    const Vc = 0.17 * axial * Math.sqrt(fc) * b * d / 1000; // kN
    const VsMax = 0.66 * Math.sqrt(fc) * b * d / 1000;
    const Vs = Math.max(0, Vu / phi - Vc);
    const Av = legs * barArea(ds);
    let sReq = Vs > 0 ? Av * fyt * d / (Vs * 1000) : Infinity;
    const sMinRule = Math.min(Av * fyt / (0.062 * Math.sqrt(fc) * b), Av * fyt / (0.35 * b));
    const sMax = Vs > 0.33 * Math.sqrt(fc) * b * d / 1000 ? Math.min(d / 4, 300) : Math.min(d / 2, 600);
    let s = Math.min(sReq, sMinRule, sMax);
    s = Math.floor(s / 25) * 25;
    const ok = Vs <= VsMax && s >= 50;
    return { Vu, Vc, phiVc: phi * Vc, Vs, VsMax, s, ok, msg: Vs > VsMax ? 'แรงเฉือนเกิน (ขยายหน้าตัด)' : '' };
  }

  function columnLayers(h, n, db, dp) {
    const Ab = barArea(db);
    let nf = Math.round(n / 4) + 1;
    if (2 * nf > n) nf = n / 2;
    const rem = n - 2 * nf;
    const m = Math.round(rem / 2);
    const layers = [{ y: dp, As: nf * Ab }];
    for (let k = 1; k <= m; k++) layers.push({ y: dp + (h - 2 * dp) * k / (m + 1), As: 2 * Ab });
    layers.push({ y: h - dp, As: nf * Ab });
    return layers;
  }

  // P-M interaction for rectangular tied column, bending about axis parallel to b
  function interaction(b, h, n, db, fc, fy, cover, dt) {
    const dp = cover + dt + db / 2;
    const layers = columnLayers(h, n, db, dp);
    const Ast = n * barArea(db), Ag = b * h;
    const b1 = beta1(fc), ety = fy / ES, dtExt = h - dp;
    const P0 = 0.85 * fc * (Ag - Ast) + fy * Ast;
    const phiPnMax = 0.65 * 0.80 * P0 / 1000; // kN
    const PtMax = -0.9 * fy * Ast / 1000;
    const raw = [{ P: PtMax, M: 0, phi: 0.9 }];
    const N = 160;
    for (let k = 0; k <= N; k++) {
      const c = h * Math.pow(10, -1.7 + 2.3 * k / N); // 0.02h ... 4h
      const a = Math.min(b1 * c, h);
      const Cc = 0.85 * fc * a * b;
      let Pn = Cc, Mn = Cc * (h / 2 - a / 2);
      layers.forEach((L) => {
        let fs = Math.max(-fy, Math.min(fy, ES * 0.003 * (c - L.y) / c));
        if (L.y < a) fs -= 0.85 * fc;
        Pn += L.As * fs; Mn += L.As * fs * (h / 2 - L.y);
      });
      const et = 0.003 * (dtExt - c) / c;
      const phi = et <= ety ? 0.65 : et >= ety + 0.003 ? 0.9 : 0.65 + 0.25 * (et - ety) / 0.003;
      raw.push({ P: phi * Pn / 1000, M: Math.max(0, phi * Mn / 1e6), phi, c });
    }
    raw.sort((p, q) => p.P - q.P);
    function Mcap(Pu) {
      if (Pu > phiPnMax || Pu < PtMax) return 0;
      for (let i = 0; i < raw.length - 1; i++) {
        const p = raw[i], q = raw[i + 1];
        if (Pu >= p.P && Pu <= q.P) return q.P === p.P ? Math.max(p.M, q.M) : p.M + (q.M - p.M) * (Pu - p.P) / (q.P - p.P);
      }
      return 0;
    }
    const curve = raw.map((p) => ({ P: Math.min(p.P, phiPnMax), M: p.M }));
    return { curve, Mcap, phiPnMax, PtMax, Ast, Ag, rho: Ast / Ag, layers };
  }

  function designMembers(model, combos) {
    const { inp, elems } = model;
    const fc = inp.fc, fy = inp.fy, fyt = inp.fyt, cover = inp.cover ?? 40, ds = inp.stirrupDb ?? 9;
    const beams = [], cols = [];
    elems.forEach((e, m) => {
      const bmm = e.b * 1000, hmm = e.h * 1000;
      if (e.type === 'beam') {
        const db = e.barDb || 16;
        const d = hmm - cover - ds - db / 2;
        let MnegL = 0, MnegR = 0, Mpos = 0, Vu = 0, gL = '', gR = '', gP = '', gV = '';
        combos.forEach((cb) => {
          const pts = memberDiagram(cb.F[m], cb.w[m], e.L, 20);
          const mL = pts[0].M, mR = pts[pts.length - 1].M;
          if (-mL > MnegL) { MnegL = -mL; gL = cb.name; }
          if (-mR > MnegR) { MnegR = -mR; gR = cb.name; }
          pts.forEach((p) => { if (p.M > Mpos) { Mpos = p.M; gP = cb.name; } });
          const v = Math.max(Math.abs(pts[0].V), Math.abs(pts[pts.length - 1].V));
          if (v > Vu) { Vu = v; gV = cb.name; }
        });
        const Ab = barArea(db);
        const fl = (Mu) => { const r = beamFlexure(Mu, bmm, d, fc, fy); r.Mu = Mu; r.nBars = Math.max(2, Math.ceil(r.As / Ab)); return r; };
        const topL = fl(MnegL), topR = fl(MnegR), bot = fl(Mpos);
        const sh = shearDesign(Vu, bmm, d, fc, fyt, ds);
        // bars must fit in one layer (simple check): clear spacing >= max(25, db)
        const nMax = Math.max(topL.nBars, topR.nBars, bot.nBars);
        const clear = (bmm - 2 * (cover + ds) - nMax * db) / Math.max(1, nMax - 1);
        const fit = clear >= Math.max(25, db);
        const ok = topL.ok && topR.ok && bot.ok && sh.ok;
        beams.push({ elem: e, d, db, topL, topR, bot, shear: sh, gov: { L: gL, R: gR, P: gP, V: gV }, fit, ok });
        e.ratio = Math.max(
          MnegL / maxMoment(bmm, d, fc, fy), MnegR / maxMoment(bmm, d, fc, fy), Mpos / maxMoment(bmm, d, fc, fy),
          sh.Vs / sh.VsMax);
        if (!ok) e.ratio = Math.max(e.ratio, 1.01);
      } else {
        const n = e.bars.n || 8, db = e.bars.db || 16;
        const ir = interaction(bmm, hmm, n, db, fc, fy, cover, ds);
        let worst = { ratio: -1 };
        let Vu = 0, NuV = 0;
        combos.forEach((cb) => {
          const f = cb.F[m];
          [[f[0], Math.abs(f[2])], [-f[3], Math.abs(f[5])]].forEach(([Pu, Mu]) => {
            const cap = ir.Mcap(Pu);
            const ratio = cap > 0 ? Mu / cap : Infinity;
            if (ratio > worst.ratio) worst = { ratio, Pu, Mu, cap, combo: cb.name };
          });
          const v = Math.abs(f[1]);
          if (v > Vu) { Vu = v; NuV = f[0]; }
        });
        const d = hmm - (cover + ds + db / 2);
        const sh = shearDesign(Vu, bmm, d, fc, fyt, ds, 2, NuV, bmm * hmm);
        const sTie = Math.floor(Math.min(16 * db, 48 * ds, Math.min(bmm, hmm), sh.s) / 25) * 25;
        const rhoOk = ir.rho >= 0.01 && ir.rho <= 0.08;
        const ok = worst.ratio <= 1 && rhoOk && sh.ok;
        cols.push({ elem: e, n, db, ir, worst, shear: sh, sTie, rhoOk, ok });
        e.ratio = Math.min(worst.ratio, 9.99);
        if (!rhoOk) e.ratio = Math.max(e.ratio, 1.01);
      }
    });
    return { beams, cols };
  }

  // moment capacity at ρ_tc (used for utilisation colouring)
  function maxMoment(b, d, fc, fy) {
    const rho = 0.85 * beta1(fc) * fc / fy * 0.375;
    const Rn = rho * fy * (1 - 0.59 * rho * fy / fc);
    return 0.9 * Rn * b * d * d / 1e6;
  }

  function envelope(res, elemIdx, n = 20) {
    const e = res.model.elems[elemIdx];
    let mx = null, mn = null;
    res.combos.forEach((cb) => {
      const grid = [];
      for (let k = 0; k <= n; k++) {
        const x = e.L * k / n;
        const f = cb.F[elemIdx], w = cb.w[elemIdx];
        grid.push({ x, N: -f[0], V: f[1] - w * x, M: -f[2] + f[1] * x - w * x * x / 2 });
      }
      if (!mx) { mx = grid.map((p) => ({ ...p })); mn = grid.map((p) => ({ ...p })); return; }
      grid.forEach((p, k) => {
        ['N', 'V', 'M'].forEach((q) => { mx[k][q] = Math.max(mx[k][q], p[q]); mn[k][q] = Math.min(mn[k][q], p[q]); });
      });
    });
    return { max: mx, min: mn };
  }

  function defaultInput() {
    const st = [];
    for (let k = 0; k < 4; k++) {
      st.push({
        h: k === 0 ? 4.0 : 3.5,
        colB: k < 2 ? 400 : 350, colH: k < 2 ? 400 : 350, colBarN: 8, colBarDb: k < 2 ? 20 : 16,
        beamB: 250, beamH: 500, beamBarDb: 16,
      });
    }
    return {
      spans: [5, 5, 5], stories: st,
      fc: 24, fy: 400, fyt: 240, cover: 40, stirrupDb: 9, gammaC: 24,
      trib: 4.0, floorDL: 3.9, floorLL: 2.0, roofDL: 3.5, roofLL: 1.0, wall: 5.0, parapet: 0,
      system: 'IMF', SDS: 0.30, SD1: 0.15, Ie: 1.0, liveFrac: 0, spectrum: 'map', saTable: [],
      wind: { on: true, V50: 25, TF: 1.0, IwU: 1.0, IwS: 0.75, terrain: 'A', factor: 1.6 },
      crackBeam: 0.35, crackCol: 0.70,
    };
  }

  return {
    analyze, buildModel, memberDiagram, envelope, interaction, beamFlexure, shearDesign,
    defaultInput, SYSTEMS, COMBOS, makeCombos, designCategory, interpSpectrum, windLoads, MAX_STORIES, MAX_BAYS, beta1, barArea,
  };
});
