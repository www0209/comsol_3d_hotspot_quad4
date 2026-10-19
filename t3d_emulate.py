"""
Test of the post-processing chain WITHOUT COMSOL (needs numpy, scipy): a small finite-element solver integrates the
rail equation that t3d_run.m gives to COMSOL - kap, G, M, the Joule bump of the fin model, 60 Chebyshev-graded cells
per segment, zero flux at the rail ends - with the HOTSPOT temperatures of segments.csv, and writes the result in the
file formats of t3d_run.m into out_<tag>/ (default out_emul/).  t3d_figures.py --out out_emul then draws every figure.

    python t3d_emulate.py [--max-rails 0] [--tag emul] [--tend-years 10]

It is NOT a COMSOL result: there is no 3-D thermal solve (T_3d_K = T_hotspot_K in the tables), no .mph, no COMSOL
pictures, and its seconds are not a COMSOL timing.  t3d_meta.json says solver = python-emulator.
Nucleation times agree with the reference solver of StackEM (stackem_result.csv, t_nuc_ref_s) to about 1 %.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

try:
    from scipy.linalg import solve_banded
except ImportError:                                                    # Thomas algorithm (slower, same result)
    def solve_banded(lu, ab, b):
        n = len(b); a = ab[2, :-1].copy(); d = ab[1].copy(); c = ab[0, 1:].copy(); x = b.copy()
        for i in range(1, n):
            w = a[i - 1] / d[i - 1]; d[i] -= w * c[i - 1]; x[i] -= w * x[i - 1]
        x[-1] /= d[-1]
        for i in range(n - 2, -1, -1):
            x[i] = (x[i] - c[i] * x[i + 1]) / d[i]
        return x

Y = 365.25 * 86400.0; KB = 1.380649e-23; QE = 1.602176634e-19
HERE = os.path.dirname(os.path.abspath(__file__))


def coefficients(xloc, L, TL, dT, Tm, lam, j, em):
    """kap and kap (G + M) at local positions xloc of a segment: the expressions of t3d_run.m (var2)."""
    axl = np.abs(xloc - 0.5)
    rch = np.exp((axl - 0.5) * lam) * (1 + np.exp(-2 * axl * lam)) / (1 + np.exp(-lam))
    rsh = np.sign(xloc - 0.5) * np.exp((axl - 0.5) * lam) * (1 - np.exp(-2 * axl * lam)) / (1 + np.exp(-lam))
    Tk = TL + dT * xloc + Tm * (1 - rch)
    dTdx = (dT - Tm * lam * rsh) / L
    kap = em["D0"] * np.exp(-em["Ea_eV"] * QE / (KB * Tk)) * em["B"] * em["Omega"] / (KB * Tk)
    Gd = QE * em["Z_eff"] * em["rho"] * j / em["Omega"]
    Md = em["Q_eV"] * QE / (em["Omega"] * Tk) * dTdx
    return kap, kap * (Gd + Md), Tk


def solve_rail(seg, nper, em, t_out, t_sub, theta=1.0):
    """seg: rows of segments.csv of one rail (L, j, T_L, T_R, T_m, Gamma, sigma_T).  Returns (sigma at t_out on every
    mesh vertex (nt, npts), vertex positions along the rail (npts,), node index of every segment end)."""
    xi = 0.5 * (1 - np.cos(np.pi * np.arange(nper + 1) / nper))           # Chebyshev vertices of one segment
    s = [0.0]; kap = []; gam = []; u0 = []; ends = [0]
    for (L, j, TL, TR, Tm, Gam, sT) in seg:
        xs = s[-1] + xi * L; s += list(xs[1:]); ends.append(len(s) - 1)
        xm = 0.5 * (xi[:-1] + xi[1:])                                       # cell mid points: coefficients per element
        k_, g_, _ = coefficients(xm, L, TL, TR - TL, Tm, L / Gam, j, em); kap += list(k_); gam += list(g_)
        u0 += [sT] * nper
    s = np.array(s); kap = np.array(kap); gam = np.array(gam); n = len(s); h = np.diff(s)
    # linear elements: K u = sum over elements kap/h [[1,-1],[-1,1]]; flux source ga: f = gam * [-1, +1] per element
    u = np.zeros(n); cnt = np.zeros(n)
    for e, v in enumerate(u0):                                              # initial stress: sigma_T of the segment, averaged at junctions
        u[e] += v; u[e + 1] += v; cnt[e] += 1; cnt[e + 1] += 1
    u /= cnt
    mass = np.zeros(n); mass[:-1] += 0.5 * h; mass[1:] += 0.5 * h           # lumped mass
    ke = kap / h; f = np.zeros(n); f[:-1] -= gam; f[1:] += gam
    diag = np.zeros(n); diag[:-1] += ke; diag[1:] += ke; off = -ke
    out = np.zeros((len(t_out), n)); out[0] = u; t = 0.0; io = 1
    targets = np.unique(np.concatenate([t_sub, t_out[1:]])); targets = targets[targets > 0]
    for tn in targets:
        dt = tn - t                                                         # theta scheme (theta = 1: implicit Euler, 0.5: Crank-Nicolson)
        ab = np.zeros((3, n)); ab[1] = mass / dt + theta * diag; ab[0, 1:] = theta * off; ab[2, :-1] = theta * off
        rhs = mass / dt * u + f - (1 - theta) * (diag * u + np.r_[0.0, off * u[:-1]] + np.r_[off * u[1:], 0.0])
        u = solve_banded((1, 1), ab, rhs); t = tn
        while io < len(t_out) and abs(t_out[io] - t) <= 1e-9 * max(t, 1.0):
            out[io] = u; io += 1
    return out, s, np.array(ends)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=None); ap.add_argument("--tag", default="emul"); ap.add_argument("--max-rails", type=int, default=0, help="per die and direction (0 = all)")
    ap.add_argument("--tend-years", type=float, default=10.0); ap.add_argument("--dt-decade", type=float, default=0.05); ap.add_argument("--horizons", default="1,2,5,10")
    ap.add_argument("--nper", type=int, default=60); ap.add_argument("--sub-decade", type=float, default=0.01, help="internal step in log10(t)"); ap.add_argument("--theta", type=float, default=0.5, help="time scheme: 0.5 Crank-Nicolson, 1 implicit Euler"); ap.add_argument("--field-pts", type=int, default=11)
    a = ap.parse_args()
    if a.data is None:
        cases = sorted(x for x in os.listdir(os.path.join(HERE, "data")) if os.path.isdir(os.path.join(HERE, "data", x)))
        if len(cases) != 1: raise SystemExit(f"expected one case folder under data/, found {cases}: pass --data")
        a.data = os.path.join(HERE, "data", cases[0])
    C = json.load(open(os.path.join(a.data, "constants.json"))); em = C["em"]; sc = em["sigma_crit"]; dies = C["dies"]
    R = np.genfromtxt(os.path.join(a.data, "rails.csv"), delimiter=",", names=True); S = np.genfromtxt(os.path.join(a.data, "segments.csv"), delimiter=",", names=True)
    HT = np.genfromtxt(os.path.join(a.data, "hotspot_T.csv"), delimiter=",", names=True)
    out = os.path.join(HERE, f"out_{a.tag}"); os.makedirs(out, exist_ok=True)
    hz = sorted(set([h for h in (float(v) for v in a.horizons.split(",")) if h <= a.tend_years] + [a.tend_years])); nh = len(hz)
    tl = np.unique(np.round(np.concatenate([[0.0], 10 ** np.arange(4, np.log10(a.tend_years * Y) + 1e-9, a.dt_decade), np.array(hz) * Y]), 3)); nt = len(tl)
    ih = [int(np.argmin(np.abs(tl - h * Y))) for h in hz]
    t_sub = 10 ** np.arange(0, np.log10(a.tend_years * Y) + 1e-9, a.sub_decade)
    keep = np.zeros(len(R), bool)
    for d in np.unique(R["die_id"]):
        for kv in (0, 1):
            m = np.where((R["die_id"] == d) & (R["kind_is_V"] == kv))[0]; m = m[np.argsort(R["index"][m], kind="stable")]
            keep[m[:a.max_rails] if a.max_rails > 0 else m] = True
    idx = np.where(keep)[0]; t0 = time.time()
    fr = open(os.path.join(out, "t3d_rail_result.csv"), "w"); fr.write("rail,die_id,kind_is_V,index,n_seg,t_nuc_3d_s" + "".join(f",smax_node_{h:g}yr_Pa" for h in hz) + "\n")
    fs = open(os.path.join(out, "t3d_rail_smax.csv"), "w"); fs.write("rail" + "".join(f",{t:.9g}" for t in tl) + "\n")
    fn = open(os.path.join(out, "t3d_nodes.csv"), "w"); fn.write("rail,node,die_id,x_m,y_m,T_3d_K,T_hotspot_K" + "".join(f",sigma_{h:g}yr_Pa" for h in hz) + "\n")
    ff = {h: open(os.path.join(out, f"t3d_field_{h:g}yr.csv"), "w") for h in hz}
    for f_ in ff.values(): f_.write("rail,seg,xloc,s_m,sigma_Pa\n")
    hist_ids = []; hist = []; tn_all = []
    for q in idx:
        r = int(R["rail"][q]); ns = int(R["n_seg"][q]); o0 = int(R["seg_offset"][q]); isV = int(R["kind_is_V"][q]); die = int(R["die_id"][q])
        rows = S[o0:o0 + ns]
        if not (rows["rail"] == r).all(): raise SystemExit(f"segments of rail {r} are not at offset {o0}")
        seg = [(float(w["L_m"]), float(w["j_A_per_m2"]), float(w["T_L_K"]), float(w["T_R_K"]), float(w["T_m_K"]), float(w["Gamma_m"]), float(w["sigma_T_Pa"])) for w in rows]
        sig, s, ends = solve_rail(seg, a.nper, em, tl, t_sub, a.theta)
        sn = sig[:, ends]; smax = sn.max(axis=1); kc = np.where(smax >= sc)[0]
        if len(kc) == 0: tn = np.inf
        elif kc[0] == 0: tn = tl[0]
        else:
            k = kc[0]; tn = (np.sqrt(tl[k - 1]) + (sc - smax[k - 1]) / (smax[k] - smax[k - 1]) * (np.sqrt(tl[k]) - np.sqrt(tl[k - 1]))) ** 2 if smax[k] > smax[k - 1] else tl[k]
        tn_all.append(tn)
        x0, y0 = float(R["x0_m"][q]), float(R["y0_m"][q])
        X = np.column_stack([np.full(ns + 1, x0), y0 + s[ends]]) if isV else np.column_stack([x0 + s[ends], np.full(ns + 1, y0)])
        Tn = np.concatenate([rows["T_L_K"], [rows["T_R_K"][-1]]])
        fr.write(f"{r},{die},{isV},{int(R['index'][q])},{ns},{tn:.9e}" + "".join(f",{smax[i]:.6e}" for i in ih) + "\n")
        fs.write(f"{r}" + "".join(f",{v:.6e}" for v in smax) + "\n")
        for k in range(ns + 1):
            fn.write(f"{r},{k},{die},{X[k, 0]:.9e},{X[k, 1]:.9e},{Tn[k]:.4f},{Tn[k]:.4f}" + "".join(f",{sn[i, k]:.6e}" for i in ih) + "\n")
        xi = np.linspace(0, 1, a.field_pts)
        for h, i in zip(hz, ih):
            for k in range(ns):
                sl = s[ends[k]] + xi * (s[ends[k + 1]] - s[ends[k]]); val = np.interp(sl, s, sig[i])
                ff[h].write("".join(f"{r},{k},{xx:.4f},{ss:.9e},{vv:.6e}\n" for xx, ss, vv in zip(xi, sl, val)))
        hist_ids += [r * 1000 + k for k in range(ns + 1)]; hist.append(sn)
    for f_ in [fr, fs, fn] + list(ff.values()): f_.close()
    H = np.concatenate(hist, axis=1)
    with open(os.path.join(out, "t3d_nodes_history.csv"), "w") as f:
        f.write("time_s" + "".join(f",{i}" for i in hist_ids) + "\n")
        for it in range(nt): f.write(f"{tl[it]:.9g}" + "".join(f",{v:.6e}" for v in H[it]) + "\n")
    with open(os.path.join(out, "t3d_die_temperature.csv"), "w") as f:                       # no thermal solve: both columns HotSpot
        f.write("die_id,x_m,y_m,T_comsol_K,T_hotspot_K\n")
        for w in HT: f.write(f"{int(w['die_id'])},{w['x_m']:.9e},{w['y_m']:.9e},{w['T_K']:.4f},{w['T_K']:.4f}\n")
    tn_all = np.array(tn_all)
    meta = dict(solver=f"python-emulator (NOT COMSOL): HotSpot temperatures, lumped-mass linear elements, theta scheme theta={a.theta}", case_name=C.get("case_name"),
                n_rails=int(len(idx)), n_segments=int(R["n_seg"][idx].sum()), nper=a.nper, times_s=tl.tolist(), horizons_years=hz, seconds=time.time() - t0,
                earliest_years=float(tn_all.min() / Y) if np.isfinite(tn_all).any() else None, n_within_tend=int(np.sum(tn_all <= a.tend_years * Y)))
    json.dump(meta, open(os.path.join(out, "t3d_meta.json"), "w"), indent=1)
    ref = os.path.join(a.data, "stackem_result.csv")
    if os.path.exists(ref):
        rr = np.genfromtxt(ref, delimiter=",", names=True); tr = rr["t_nuc_ref_s"][R["rail"][idx].astype(int)]; m = np.isfinite(tr) & np.isfinite(tn_all)
        if m.any():
            rel = np.abs(tn_all[m] - tr[m]) / tr[m]
            print(f"nucleation time against the reference solver on {m.sum()} rails: median {np.median(rel):.3%}, max {rel.max():.3%}; "
                  f"nucleating within {a.tend_years:g} yr: emulator {int(np.sum(tn_all <= a.tend_years * Y))}, reference {int(np.sum(tr <= a.tend_years * Y))}")
    print(f"emulator: {len(idx)} rails in {time.time() - t0:.1f} s; earliest {meta['earliest_years']} yr; wrote {out}")


if __name__ == "__main__":
    main()
