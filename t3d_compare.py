"""
Figures and numbers for the 3-D COMSOL run of one case (needs only numpy and matplotlib).

    python t3d_compare.py --out out_run

Reads   out_<tag>/t3d_die_temperature.csv, t3d_nodes.csv, t3d_rail_result.csv   (t3d_run.m)
        data/<case>/stackem_result.csv  reference solver and StackEM nucleation times with the HotSpot temperatures
Writes  out_<tag>/fig_t3d_temperature.png   every die: COMSOL temperature, and COMSOL minus HotSpot
        out_<tag>/fig_t3d_stress_<H>yr.png  stress on every rail of every die
        out_<tag>/fig_t3d_tnuc.png          nucleation time: 3-D COMSOL chain against the HotSpot-based reference
        out_<tag>/t3d_compare.json
A difference in nucleation time here is mostly a difference between the two THERMAL models: near 345 K one kelvin
changes the diffusivity, and with it the nucleation time, by about 8 %.
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

Y = 365.25 * 86400.0
INK, INK2, MUTED, GRID, SURF, CRIT = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#d03b3b"
SEQ = LinearSegmentedColormap.from_list("seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#ec835a", "#b8431f"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF, "axes.titlesize": 9, "axes.titlecolor": INK, "legend.frameon": False})
here = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--out", default="out_run"); ap.add_argument("--data", default=None, help="data/<case> (default: the one folder under data/)")
a = ap.parse_args(); d = a.out if os.path.isabs(a.out) else os.path.join(here, a.out)
if a.data is None:
    cases = sorted(x for x in os.listdir(os.path.join(here, "data")) if os.path.isdir(os.path.join(here, "data", x)))
    if len(cases) != 1: raise SystemExit(f"expected one case folder under data/, found {cases}: pass --data")
    a.data = os.path.join(here, "data", cases[0])
C = json.load(open(os.path.join(a.data, "constants.json"))); dies = C["dies"]; sc = C["em"]["sigma_crit"]; case = C.get("case_name", C.get("case", "case"))
N0 = np.genfromtxt(os.path.join(d, "t3d_nodes.csv"), delimiter=",", names=True, max_rows=2)
hz = [float(c[6:-5].replace("_", ".")) for c in N0.dtype.names if c.startswith("sigma_") and c.endswith("yr_Pa")]
yr = lambda h: f"{h:g} year" + ("" if h == 1 else "s")
res = {}
# ---------------- temperature ----------------
T = np.genfromtxt(os.path.join(d, "t3d_die_temperature.csv"), delimiter=",", names=True)
fig, axs = plt.subplots(2, len(dies), figsize=(3.3 * len(dies) + 1.0, 6.2), squeeze=False)
lim = float(np.max(np.abs(T["T_comsol_K"] - T["T_hotspot_K"]))); res["temperature"] = {}
for k, name in enumerate(dies):
    m = T["die_id"] == k; nx = len(np.unique(np.round(T["x_m"][m], 12))); ny = len(np.unique(np.round(T["y_m"][m], 12))); o = np.lexsort((T["x_m"][m], T["y_m"][m]))
    x = T["x_m"][m][o].reshape(ny, nx) * 1e3; y = T["y_m"][m][o].reshape(ny, nx) * 1e3
    tc = T["T_comsol_K"][m][o].reshape(ny, nx); dT = tc - T["T_hotspot_K"][m][o].reshape(ny, nx)
    p0 = axs[0][k].pcolormesh(x, y, tc, cmap=SEQ, shading="nearest", vmin=T["T_comsol_K"].min(), vmax=T["T_comsol_K"].max())
    p1 = axs[1][k].pcolormesh(x, y, dT, cmap=DIV, shading="nearest", vmin=-lim, vmax=lim)
    axs[0][k].set_title(f"{name}: COMSOL {tc.min():.1f} to {tc.max():.1f} K", loc="left")
    axs[1][k].set_title(f"{name}: mean {dT.mean():+.2f} K, largest {np.abs(dT).max():.2f} K", loc="left")
    for r in (0, 1): axs[r][k].set_aspect("equal")
    axs[1][k].set_xlabel("x (mm)")
    res["temperature"][name] = dict(comsol_min=float(tc.min()), comsol_max=float(tc.max()), diff_mean=float(dT.mean()), diff_max_abs=float(np.abs(dT).max()))
axs[0][0].set_ylabel("y (mm)"); axs[1][0].set_ylabel("y (mm)")
c0 = fig.colorbar(p0, ax=list(axs[0]), fraction=0.03, pad=0.02); c0.set_label("COMSOL temperature (K)"); c0.outline.set_visible(False)
c1 = fig.colorbar(p1, ax=list(axs[1]), fraction=0.03, pad=0.02); c1.set_label("COMSOL minus HotSpot (K)"); c1.outline.set_visible(False)
fig.suptitle(f"{case}: die temperature from the 3-D COMSOL model (top) and its difference from HotSpot (bottom)", x=0.01, y=0.99, ha="left", color=INK, fontsize=10)
fig.savefig(os.path.join(d, "fig_t3d_temperature.png"), dpi=200, bbox_inches="tight"); plt.close(fig)
# ---------------- stress maps ----------------
N = np.genfromtxt(os.path.join(d, "t3d_nodes.csv"), delimiter=",", names=True)
RR = np.atleast_1d(np.genfromtxt(os.path.join(d, "t3d_rail_result.csv"), delimiter=",", names=True))
for h in hz:
    col = f"sigma_{h:g}yr_Pa"; fig, axs = plt.subplots(1, len(dies), figsize=(3.3 * len(dies) + 0.9, 3.6), squeeze=False)
    allv = []; packs = []
    for k in range(len(dies)):
        segs, val = [], []
        for r in np.unique(N["rail"][N["die_id"] == k]):
            m = N["rail"] == r; P = np.stack([N["x_m"][m], N["y_m"][m]], 1) * 1e3; s_ = N[col][m] / 1e6
            segs.append(np.stack([P[:-1], P[1:]], 1)); val.append(np.maximum(s_[:-1], s_[1:]))
        packs.append((np.concatenate(segs), np.concatenate(val)) if segs else None)
        if segs: allv.append(packs[-1][1])
    allv = np.concatenate(allv); sub = allv[allv < sc / 1e6]; vmin = float(np.percentile(sub, 2)) if len(sub) else 0.0
    for k, ax in enumerate(axs[0]):
        if packs[k] is None: ax.set_axis_off(); continue
        sg, v = packs[k]; o = np.argsort(v)
        lc = LineCollection(sg[o], array=v[o], cmap=SEQ, linewidths=1.1, capstyle="butt"); lc.set_clim(vmin, sc / 1e6); ax.add_collection(lc)
        if (v >= sc / 1e6).any(): ax.add_collection(LineCollection(sg[v >= sc / 1e6], colors=CRIT, linewidths=1.8, capstyle="butt", zorder=3))
        nf = int(np.sum((RR["die_id"] == k) & (RR["t_nuc_3d_s"] <= h * Y)))
        ax.set_title(f"{dies[k]}: {nf} of {int(np.sum(RR['die_id'] == k))} rails nucleated", loc="left"); ax.autoscale(); ax.set_aspect("equal", adjustable="datalim"); ax.set_xlabel("x (mm)")
    axs[0][0].set_ylabel("y (mm)")
    cb = fig.colorbar(lc, ax=list(axs[0]), fraction=0.03, pad=0.02); cb.set_label("larger end-node stress of the segment (MPa)"); cb.outline.set_visible(False)
    fig.suptitle(f"{case}, 3-D COMSOL chain: stress on the rails after {yr(h)}. Red: segments with an end at or above {sc/1e6:.0f} MPa", x=0.01, y=1.04, ha="left", color=INK, fontsize=10)
    fig.savefig(os.path.join(d, f"fig_t3d_stress_{h:g}yr.png"), dpi=200, bbox_inches="tight"); plt.close(fig)
# ---------------- nucleation time against the HotSpot-based reference ----------------
ref = np.genfromtxt(os.path.join(a.data, "stackem_result.csv"), delimiter=",", names=True)
rid = RR["rail"].astype(int); t3 = RR["t_nuc_3d_s"]; tr = ref["t_nuc_ref_s"][rid]; ts = ref["t_nuc_stackem_s"][rid]
m = np.isfinite(t3) & np.isfinite(tr)
res["tnuc"] = dict(n_rails=int(len(rid)), n_nucleating_3d=int(np.isfinite(t3).sum()), n_nucleating_reference=int(np.isfinite(tr).sum()),
                   within_10yr_3d=int(np.sum(t3 <= 10 * Y)), within_10yr_reference=int(np.sum(tr <= 10 * Y)),
                   earliest_years_3d=float(t3.min() / Y) if np.isfinite(t3).any() else None, earliest_years_reference=float(np.nanmin(np.where(np.isfinite(tr), tr, np.nan)) / Y) if np.isfinite(tr).any() else None)
if m.any():
    rel = (t3[m] - tr[m]) / tr[m]; res["tnuc"].update(rel_diff_median=float(np.median(rel)), rel_diff_abs_median=float(np.median(np.abs(rel))), rel_diff_abs_max=float(np.abs(rel).max()))
    fig, ax = plt.subplots(figsize=(4.2, 4.0)); lo = min(t3[m].min(), tr[m].min()) / Y; hi = max(t3[m].max(), tr[m].max()) / Y
    for k, (name, mk, colr) in enumerate(zip(dies, "osD^v", ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"])):
        mm = m & (RR["die_id"] == k)
        if mm.any(): ax.scatter(tr[mm] / Y, t3[mm] / Y, s=16, marker=mk, color=colr, edgecolor=SURF, linewidths=0.5, label=f"{name} ({mm.sum()} rails)", zorder=3)
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8, zorder=1); ax.set_xscale("log"); ax.set_yscale("log")
    tk = [t for t in (0.5, 1, 2, 5, 10, 20, 50, 100) if lo * 0.8 <= t <= hi * 1.25]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(tk)); axis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g")); axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid(True, color=GRID, lw=0.5); ax.set_xlabel("reference solver with HotSpot temperatures (years)"); ax.set_ylabel("3-D COMSOL chain (years)")
    ax.set_title(f"{case}: nucleation time of every rail that nucleates in both", loc="left"); ax.legend(loc="upper left")
    fig.savefig(os.path.join(d, "fig_t3d_tnuc.png"), dpi=200, bbox_inches="tight"); plt.close(fig)
dTn = N["T_3d_K"] - N["T_hotspot_K"]; res["rail_node_temperature_diff_K"] = dict(mean=float(dTn.mean()), max_abs=float(np.abs(dTn).max()))
json.dump(res, open(os.path.join(d, "t3d_compare.json"), "w"), indent=1); print(json.dumps(res, indent=1)); print("figures written to", d)
