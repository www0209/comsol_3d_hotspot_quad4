"""
All result figures of a 3-D run (needs numpy and matplotlib; emfig.py must be next to this file).

    python t3d_figures.py --out out_run

Reads   out_<tag>/t3d_nodes.csv, t3d_rail_result.csv, t3d_rail_smax.csv, t3d_die_temperature.csv   (t3d_run.m)
        data/<case>/constants.json, stack.json, stackem_result.csv   (the one case folder next to this file)
Writes  out_<tag>/figures/ :
    fig_exploded_tnuc              exploded stack: COMSOL die temperature on the slabs, rails coloured by nucleation time
    fig_exploded_stress_<H>yr      the same view with every segment coloured by its stress at horizon H
    fig_pdn_stress_all             stress map of the power-delivery rails: every die x every horizon
    fig_pdn_stress_<die>_<H>yr     one sheet per die and horizon: horizontal layer, vertical layer, both
    fig_pdn_tnuc                   nucleation time of every rail, per die
    fig_pdn_temperature            temperature seen by the rails, per die
    fig_smax_history               largest stress of every rail against time, per die
    fig_tnuc_cdf                   nucleated rails against time, per die, with the reference solver and StackEM
    fig_stress_hist                distribution of the rail stress per die and horizon
    fig_rail_profiles_<die>        stress along the three earliest rails of every die
    fig_tnuc_parity_all            nucleation time against the reference solver and StackEM
and then runs t3d_compare.py (temperature against HotSpot, t3d_compare.json).
"""
import argparse, json, os, subprocess, sys
import numpy as np
import emfig

here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--out", default="out_run"); ap.add_argument("--data", default=None, help="data/<case> (default: the one folder under data/)")
a = ap.parse_args(); d = a.out if os.path.isabs(a.out) else os.path.join(here, a.out)
if a.data is None:
    cases = sorted(x for x in os.listdir(os.path.join(here, "data")) if os.path.isdir(os.path.join(here, "data", x)))
    if len(cases) != 1: raise SystemExit(f"expected one case folder under data/, found {cases}: pass --data")
    a.data = os.path.join(here, "data", cases[0])
C = json.load(open(os.path.join(a.data, "constants.json"))); K = json.load(open(os.path.join(a.data, "stack.json")))
N = np.genfromtxt(os.path.join(d, "t3d_nodes.csv"), delimiter=",", names=True)
RR = np.atleast_1d(np.genfromtxt(os.path.join(d, "t3d_rail_result.csv"), delimiter=",", names=True))
SM = np.atleast_2d(np.genfromtxt(os.path.join(d, "t3d_rail_smax.csv"), delimiter=",", skip_header=1))
times = np.array([float(v) for v in open(os.path.join(d, "t3d_rail_smax.csv")).readline().strip().split(",")[1:]])
T = np.genfromtxt(os.path.join(d, "t3d_die_temperature.csv"), delimiter=",", names=True)
ref = np.genfromtxt(os.path.join(a.data, "stackem_result.csv"), delimiter=",", names=True)
hz = [float(c[6:-5].replace("_", ".")) for c in N.dtype.names if c.startswith("sigma_") and c.endswith("yr_Pa")]; rid = RR["rail"].astype(int); order = np.argsort(N["rail"], kind="stable")
cut = np.r_[0, np.where(np.diff(N["rail"][order]) != 0)[0] + 1, len(order)]; grp = {int(N["rail"][order[c0]]): order[c0:c1] for c0, c1 in zip(cut[:-1], cut[1:])}
smrow = {int(r): k for k, r in enumerate(SM[:, 0])}
Tg = {}
for k in range(len(C["dies"])):
    m = T["die_id"] == k
    if m.any():
        nx = len(np.unique(np.round(T["x_m"][m], 12))); ny = len(np.unique(np.round(T["y_m"][m], 12))); o = np.lexsort((T["x_m"][m], T["y_m"][m]))
        if nx * ny != m.sum(): raise SystemExit(f"die {k}: {m.sum()} temperature cells are not a {ny} x {nx} grid")
        Tg[k] = (T["x_m"][m][o].reshape(ny, nx) * 1e3, T["y_m"][m][o].reshape(ny, nx) * 1e3, T["T_comsol_K"][m][o].reshape(ny, nx))
D = dict(case=C.get("case_name", C.get("case", "case")), label="COMSOL 3-D, HotSpot T", dies=C["dies"], hz=hz, sc=C["em"]["sigma_crit"], rail=rid, die=RR["die_id"].astype(int),
         kind=RR["kind_is_V"].astype(int), tn=RR["t_nuc_3d_s"].astype(float), times=times, smax=np.array([SM[smrow[r], 1:] for r in rid]),
         nodes=[np.column_stack([N["x_m"][grp[r]], N["y_m"][grp[r]]]) for r in rid],
         sig=[np.column_stack([N[f"sigma_{h:g}yr_Pa"][grp[r]] for h in hz]) for r in rid], Tn=[N["T_3d_K"][grp[r]] for r in rid], Tgrid=Tg,
         die_xy=(K["die_size_x"] * 1e3, K["die_size_y"] * 1e3),
         T_source="HotSpot temperature (no thermal solve)",
         ref={"reference solver, HotSpot temperature": ref["t_nuc_ref_s"][rid], "StackEM, HotSpot temperature": ref["t_nuc_stackem_s"][rid]})
made = emfig.all_figures(D, os.path.join(d, "figures"))
print(f"{len(made)} figures written to {os.path.join(d, 'figures')}")
subprocess.call([sys.executable, os.path.join(here, "t3d_compare.py"), "--out", d, "--data", a.data])
