// Run: node rc-frame/test/engine.test.js
const assert = require('assert');
const RC = require('../engine.js');
let pass = 0;
const near = (a, b, tol, msg) => { assert(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg}: ${a} vs ${b}`); pass++; };

function base() {
  const inp = RC.defaultInput();
  inp.stories = [{ h: 3, colB: 300, colH: 300, colBarN: 8, colBarDb: 16, beamB: 250, beamH: 500, beamBarDb: 16 }];
  inp.spans = [6];
  inp.gammaC = 0; inp.wall = 0; inp.parapet = 0;
  return inp;
}

// 1) Portal frame, very stiff beam, lateral load only -> each column acts fixed-fixed: V=12EI d/h^3
{
  const inp = base();
  inp.stories[0].beamH = 5000; inp.stories[0].beamB = 5000; // ~rigid beam
  inp.roofDL = 10; inp.roofLL = 0;
  const r = RC.analyze(inp);
  const F = r.eq.F[1];
  const e = r.model.elems[0];
  const dExp = (F / 2) * 3 ** 3 / (12 * e.E * e.I);
  near(r.cases.E.U[r.model.nodes[2].id * 3], dExp, 0.01, 'portal sway (rigid beam)');
  // base moments = V h / 2
  near(Math.abs(r.cases.E.F[0][2]), F / 2 * 3 / 2, 0.01, 'base moment rigid beam');
  // equilibrium: sum base shear = F
  const Rx = r.cases.E.R[0] + r.cases.E.R[3];
  near(-Rx, F, 1e-6, 'horizontal equilibrium');
}

// 2) Gravity on beam of portal with infinitely stiff columns -> fixed-end moments wL^2/12
{
  const inp = base();
  inp.stories[0].colB = 5000; inp.stories[0].colH = 5000;
  inp.roofLL = 0; inp.roofDL = 10; inp.trib = 1; // w=10
  const r = RC.analyze(inp);
  const bi = r.model.elems.findIndex((e) => e.type === 'beam');
  const pts = RC.memberDiagram(r.cases.D.F[bi], r.cases.D.w[bi], 6, 20);
  near(pts[0].M, -10 * 36 / 12, 0.01, 'fixed-end moment');
  near(pts[10].M, 10 * 36 / 24, 0.01, 'midspan moment');
  // vertical equilibrium
  near(r.cases.D.R[1] + r.cases.D.R[4], 60, 1e-6, 'vertical reactions');
}

// 3) Seismic coefficients
{
  const inp = RC.defaultInput();
  const r = RC.analyze(inp);
  const H = 4 + 3 * 3.5;
  near(r.eq.Ta, 0.02 * H, 1e-9, 'Ta');
  const Ts = inp.SD1 / inp.SDS;
  const Sa = r.eq.Ta <= Ts ? inp.SDS : inp.SD1 / r.eq.Ta;
  near(r.eq.Cs, Math.max(Sa * 1 / 5, 0.044 * inp.SDS), 1e-9, 'Cs');
  near(r.eq.F.reduce((a, b) => a + b, 0), r.eq.V, 1e-9, 'sum Fx = V');
  // global horizontal equilibrium
  let Rx = 0; r.model.nodes.forEach((n) => { if (n.fixed) Rx += r.cases.E.R[n.id * 3]; });
  near(-Rx, r.eq.V, 1e-6, 'base shear reaction');
  // global vertical equilibrium for D
  let Ry = 0; r.model.nodes.forEach((n) => { if (n.fixed) Ry += r.cases.D.R[n.id * 3 + 1]; });
  const Wd = r.eq.WD.reduce((a, b) => a + b, 0);
  near(Ry, Wd, 1e-6, 'dead load reaction');
  assert(r.design.beams.length === 12 && r.design.cols.length === 16); pass++;
  console.log('  T(code)=%s s, T(Rayleigh)=%s s, Cs=%s, V=%s kN', r.eq.Ta.toFixed(3), r.eq.Tray.toFixed(3), r.eq.Cs.toFixed(4), r.eq.V.toFixed(1));
}

// 4) Beam flexure vs hand calc: b=250 d=440 fc=24 fy=400 Mu=100
{
  const f = RC.beamFlexure(100, 250, 440, 24, 400);
  const Rn = 100e6 / (0.9 * 250 * 440 * 440);
  const rho = 0.85 * 24 / 400 * (1 - Math.sqrt(1 - 2 * Rn / (0.85 * 24)));
  near(f.As, rho * 250 * 440, 1e-9, 'beam As');
}

// 5) Column interaction: pure axial and pure tension points
{
  const ir = RC.interaction(400, 400, 8, 20, 24, 400, 40, 9);
  const Ast = 8 * Math.PI * 100;
  near(ir.phiPnMax, 0.65 * 0.8 * (0.85 * 24 * (160000 - Ast) + 400 * Ast) / 1000, 1e-9, 'phiPn,max');
  near(ir.PtMax, -0.9 * 400 * Ast / 1000, 1e-9, 'phiPn tension');
  assert(ir.Mcap(0) > 100 && ir.Mcap(0) < 250, 'Mcap(0) plausible ' + ir.Mcap(0)); pass++;
  assert(ir.Mcap(1000) > ir.Mcap(0), 'M increases with moderate axial'); pass++;
}

// 6) 8 storeys, 10 bays runs
{
  const inp = RC.defaultInput();
  inp.spans = new Array(10).fill(6);
  inp.stories = new Array(8).fill(0).map(() => ({ ...inp.stories[0] }));
  const t = Date.now();
  const r = RC.analyze(inp);
  assert(r.stories.length === 8); pass++;
  console.log('  8 storeys x 10 bays solved in %d ms', Date.now() - t);
  assert.throws(() => { inp.stories.push({ ...inp.stories[0] }); RC.analyze(inp); }); pass++;
}
console.log(`OK — ${pass} checks passed`);

// 7) Seismic design category & vertical seismic in combos
{
  assert.strictEqual(RC.designCategory(0.30, 0.15, 1.0), 2); // ค
  assert.strictEqual(RC.designCategory(0.40, 0.30, 1.0), 3); // ง
  assert.strictEqual(RC.designCategory(0.10, 0.05, 1.0), 0); // ก
  assert.strictEqual(RC.designCategory(0.20, 0.10, 1.5), 2); // IV bumps ข -> ค
  const inp = RC.defaultInput(); inp.SDS = 0.4; inp.SD1 = 0.3; inp.system = 'IMF';
  const r = RC.analyze(inp);
  assert(!r.eq.systemOk && r.eq.minSystem === 'SMF');
  near(r.combos[2].name === '1.28D+1.0L+1.0E' ? 1 : 0, 1, 0, 'Ev factor in combo name');
  pass += 5;
  console.log('OK — SDC/Ev checks passed');
}

// 8) Wind (DPT 1311-50 simplified): hand calc of the top-level force
{
  const inp = RC.defaultInput();
  inp.wind = { on: true, V50: 25, TF: 1.2, IwU: 1.0, IwS: 0.75, terrain: 'A', factor: 1.6 };
  const r = RC.analyze(inp);
  const q = 0.5 * 1.25 * (1.2 * 25) ** 2 / 1000;
  const H = 14.5, D = 15, hd = H / D;
  const CpL = -0.3 - 0.2 * (hd - 0.25) / 0.75;
  const pW = q * Math.max(0.9, (H / 10) ** 0.2) * 2 * 0.8;
  const pL = q * Math.max(0.9, (H / 2 / 10) ** 0.2) * 2 * -CpL;
  near(r.wind.q, q, 1e-9, 'q');
  near(r.wind.F[4], (pW + pL) * (3.5 / 2) * 4, 1e-9, 'roof wind force');
  let Rx = 0; r.model.nodes.forEach((n) => { if (n.fixed) Rx += r.cases.W.R[n.id * 3]; });
  near(-Rx, r.wind.V, 1e-6, 'wind base shear');
  assert(r.combos.length === 10 && r.combos[6].name === '1.2D+1.0L+1.6W'); pass++;
  near(r.wind.svc, 0.75 / 1.44, 1e-12, 'serviceability ratio');
  console.log('  wind V = %s kN, top drift (service) = H/%s', r.wind.V.toFixed(1), (1 / r.wind.topRatio).toFixed(0));
}

// 9) User spectrum table (Bangkok basin zones 7–10)
{
  near(RC.interpSpectrum([[0, 0.2], [1, 0.4]], 0.5), 0.3, 1e-12, 'interp mid');
  near(RC.interpSpectrum([[0, 0.2], [1, 0.4]], 3), 0.4, 1e-12, 'interp clamp');
  const inp = RC.defaultInput();
  inp.spectrum = 'table'; inp.saTable = [[0, 0.1], [0.29, 0.1], [1, 0.3]];
  const r = RC.analyze(inp);
  near(r.eq.Sa, 0.1, 1e-9, 'Sa from table');
  inp.saTable = [];
  assert.throws(() => RC.analyze(inp)); pass++;
}
console.log(`OK — ${pass} checks passed (all)`);
