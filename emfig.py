"""
emfig.py - result figures shared by the 1-D / 2-D runs (comsol_full/cf_figures.py) and the 3-D run (comsol_3d/t3d_figures.py).
Needs numpy and matplotlib only.

The exploded views use the visual grammar of the paper's stack figure (stackem/viz/stack3d.py: orthographic view,
elev 24, azim -58, die slabs textured with the temperature in "inferno", rails in "viridis_r" on a log scale of the
nucleation time, earliest rail of every die outlined) so that they can be put next to it.

A run is described by one dict D:
    dies      list of die names, bottom to top           hz        horizons in years          sc   critical stress (Pa)
    label     what solved it, e.g. "COMSOL 1-D"          case      case name
    rail      (n,) rail ids       die (n,)  kind (n,) 0 = horizontal, 1 = vertical            tn   (n,) nucleation time (s), inf = none
    nodes     list of n arrays (S+1, 2): x, y in metres  sig       list of n arrays (S+1, nh): node stress (Pa) at the horizons
    Tn        list of n arrays (S+1,): node temperature (K) or None
    times     (m,) output times (s)                      smax      (n, m) max node stress of every rail at every output time
    ref       {name: (n,) nucleation time (s)} of other methods on the same rails (may be empty)
    Tgrid     {die index: (x (r,c) mm, y (r,c) mm, T (r,c) K)} for the slab texture; built from Tn when missing
    die_xy    (size_x, size_y) of the dies in mm (optional; the rail extent otherwise)
    field     {rail id: (s (p,) m, sigma (p, nh) Pa)} stress along the rail (optional)
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

Y = 365.25 * 86400.0
INK, INK2, MUTED, GRID, SURF, CRIT = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#d03b3b"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
SEQ = LinearSegmentedColormap.from_list("stress", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
CMAP_T, CMAP_TNUC, SILICON, BOND, VERMILION = "inferno", "viridis_r", "#8E9BA6", "#D9C9A3", "#D55E00"


def style():
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
                         "savefig.facecolor": SURF, "axes.titlesize": 9, "axes.titlecolor": INK, "legend.frameon": False})


def save(fig, d, name, made):
    fig.savefig(os.path.join(d, name + ".png"), dpi=220, bbox_inches="tight"); fig.savefig(os.path.join(d, name + ".pdf"), bbox_inches="tight")
    plt.close(fig); made.append(name + ".png")


def yr(h):
    return f"{h:g} year" + ("" if h == 1 else "s")


def seg_values(D, i, q):
    """Segments of rail i in mm and the larger end-node stress of every segment (MPa) at horizon q."""
    P = D["nodes"][i] * 1e3; s = D["sig"][i][:, q] / 1e6
    return np.stack([P[:-1], P[1:]], 1), np.maximum(s[:-1], s[1:])


def stress_limits(D, q):
    v = np.concatenate([seg_values(D, i, q)[1] for i in range(len(D["rail"]))]); sub = v[v < D["sc"] / 1e6]
    return (float(np.percentile(sub, 2)) if len(sub) else 0.0), D["sc"] / 1e6


def temperature_grids(D, n=64):
    """Die temperature on an n x n grid for the slab texture: given, or interpolated from the rail-node temperatures."""
    if D.get("Tgrid"):
        return D["Tgrid"]
    out = {}
    for k in range(len(D["dies"])):
        ii = [i for i in range(len(D["rail"])) if D["die"][i] == k and D["Tn"][i] is not None]
        if not ii:
            continue
        P = np.concatenate([D["nodes"][i] for i in ii]) * 1e3; T = np.concatenate([D["Tn"][i] for i in ii])
        key = np.round(P, 6); _, u = np.unique(key, axis=0, return_index=True); P, T = P[u], T[u]
        if len(P) < 4 or np.ptp(P[:, 0]) == 0 or np.ptp(P[:, 1]) == 0:
            continue
        gx = np.linspace(P[:, 0].min(), P[:, 0].max(), n); gy = np.linspace(P[:, 1].min(), P[:, 1].max(), n); X, Yg = np.meshgrid(gx, gy)
        Z = mtri.LinearTriInterpolator(mtri.Triangulation(P[:, 0], P[:, 1]), T)(X, Yg)
        out[k] = (X, Yg, np.ma.filled(Z, np.nanmean(T)))
    return out


# =====================================================================================================================
# exploded views
# =====================================================================================================================
def _slab(ax, S, z, t, zo, Sy=None):
    Sy = S if Sy is None else Sy; x = [0, S]; y = [0, Sy]; f = []
    for a in (0, 1):
        f.append([(x[a], 0, z), (x[a], Sy, z), (x[a], Sy, z + t), (x[a], 0, z + t)]); f.append([(0, y[a], z), (S, y[a], z), (S, y[a], z + t), (0, y[a], z + t)])
    f.append([(0, 0, z), (S, 0, z), (S, Sy, z), (0, Sy, z)])
    ax.add_collection3d(Poly3DCollection(f, facecolors=SILICON, edgecolors="#333333", linewidths=0.3, alpha=0.95, zorder=zo()))


def exploded(D, d, made, mode="tnuc", q=None, name=None, elev=24.0, azim=-58.0):
    """mode 'tnuc': rails coloured by nucleation time (the paper's figure); mode 'stress': segments coloured by stress at horizon q."""
    style(); nd = len(D["dies"]); n = len(D["rail"]); G = temperature_grids(D)
    P_ = np.concatenate(D["nodes"]) * 1e3
    S, Sy = (float(D["die_xy"][0]), float(D["die_xy"][1])) if D.get("die_xy") else (max(float(P_[:, 0].max()), 1e-9), max(float(P_[:, 1].max()), 1e-9))
    gap = (0.26 + 0.05 * max(nd - 3, 0)) * max(S, Sy); thick = 0.03 * max(S, Sy)
    fig = plt.figure(figsize=(7.16, 3.9 + 0.5 * max(nd - 3, 0))); ax = fig.add_axes([0.02, 0.0, 0.80, 1.0], projection="3d")
    ax.set_proj_type("ortho"); ax.computed_zorder = False; cnt = [0]

    def zo():
        cnt[0] += 1; return cnt[0]
    Tall = np.concatenate([g[2].ravel() for g in G.values()]) if G else np.array([0.0, 1.0])
    tnorm = Normalize(float(Tall.min()), float(Tall.max())); tn_norm = LogNorm(0.03, 100.0); cm_tn = plt.get_cmap(CMAP_TNUC)
    if mode == "stress":
        vmin, vmax = stress_limits(D, q); snorm = Normalize(vmin, vmax)
    labels = []; z = 0.0
    for k, dname in enumerate(D["dies"]):
        _slab(ax, S, z, thick, zo, Sy); zt = z + thick
        if k in G:
            X, Yg, T = G[k]
            ax.plot_surface(X, Yg, np.full_like(X, zt), facecolors=plt.get_cmap(CMAP_T)(tnorm(T)), rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False, zorder=zo())
        ii = [i for i in range(n) if D["die"][i] == k]; zr = zt + 0.02 * max(S, Sy); tn = np.array([D["tn"][i] for i in ii])
        if mode == "tnuc":
            imm = [[(*D["nodes"][i][0] * 1e3, zr), (*D["nodes"][i][-1] * 1e3, zr)] for i in ii if not np.isfinite(D["tn"][i])]
            mor = [[(*D["nodes"][i][0] * 1e3, zr), (*D["nodes"][i][-1] * 1e3, zr)] for i in ii if np.isfinite(D["tn"][i])]
            col = [cm_tn(tn_norm(min(max(D["tn"][i] / Y, 0.03), 100.0))) for i in ii if np.isfinite(D["tn"][i])]
            if imm: ax.add_collection3d(Line3DCollection(imm, colors=[(1, 1, 1, 0.3)], linewidths=0.35, zorder=zo()))
            if mor:
                ax.add_collection3d(Line3DCollection(mor, colors=[(1, 1, 1, 0.9)], linewidths=2.2, zorder=zo()))
                ax.add_collection3d(Line3DCollection(mor, colors=col, linewidths=1.2, zorder=zo()))
        else:
            sg, v = [], []
            for i in ii:
                a, b = seg_values(D, i, q); sg.append(a); v.append(b)
            if sg:
                sg = np.concatenate(sg); v = np.concatenate(v); o = np.argsort(v); s3 = np.concatenate([sg, np.full(sg.shape[:2] + (1,), zr)], 2)
                ax.add_collection3d(Line3DCollection(s3[o], colors=SEQ(snorm(v[o])), linewidths=1.0, zorder=zo()))
                if (v >= D["sc"] / 1e6).any():
                    ax.add_collection3d(Line3DCollection(s3[v >= D["sc"] / 1e6], colors=[CRIT], linewidths=1.7, zorder=zo()))
        if len(tn) and np.isfinite(tn).any():
            i0 = ii[int(np.argmin(tn))]; e = [[(*D["nodes"][i0][0] * 1e3, zt + 0.04 * max(S, Sy)), (*D["nodes"][i0][-1] * 1e3, zt + 0.04 * max(S, Sy))]]
            ax.add_collection3d(Line3DCollection(e, colors=[VERMILION], linewidths=2.8, capstyle="round", zorder=zo()))
            ax.add_collection3d(Line3DCollection(e, colors=["white"], linewidths=0.8, capstyle="round", zorder=zo()))
            nv = int(np.sum(tn <= 10 * Y)); body = f"earliest {tn.min() / Y:.1f} yr, {nv} rail{'s' if nv != 1 else ''} < 10 yr"
        else:
            body = "no rail nucleates"
        labels.append((dname + f"  ({len(ii)} rails)", body, zt))
        if k < nd - 1:
            zb = z + 0.55 * gap
            ax.add_collection3d(Poly3DCollection([[(0, 0, zb), (S, 0, zb), (S, Sy, zb), (0, Sy, zb)]], facecolors=BOND, edgecolors="#B8A67A", linewidths=0.4, alpha=0.18, zorder=zo()))
        z += gap
    zs = z - 0.6 * gap
    ax.add_collection3d(Line3DCollection([[(0, 0, zs), (S, 0, zs), (S, Sy, zs), (0, Sy, zs), (0, 0, zs)]], colors=["#7F8B95"], linewidths=0.9, zorder=zo()))
    ax.set_xlim(0, S); ax.set_ylim(0, Sy); ax.set_zlim(0, z - 0.35 * gap); ax.set_box_aspect((1, Sy / S, 0.55 + 0.12 * nd)); ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x (mm)", labelpad=-9); ax.set_ylabel("y (mm)", labelpad=-9); ax.set_zlabel("")
    ax.set_xticks([0, S]); ax.set_xticklabels(["0", f"{S:g}"]); ax.set_yticks([0, Sy]); ax.set_yticklabels(["0", f"{Sy:g}"]); ax.set_zticks([]); ax.grid(False); ax.tick_params(pad=-4, labelsize=6.5)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False; axis.pane.set_edgecolor((1, 1, 1, 0))
    ax.zaxis.line.set_color((1, 1, 1, 0))
    for head, body, zt in labels:
        x2, y2, _ = proj3d.proj_transform(0.0, 0.0, zt, ax.get_proj())
        ax.annotate(head + "\n" + body, (x2, y2), xycoords="data", xytext=(-6, 0), textcoords="offset points", fontsize=6.3, color="#222222", ha="right", va="center", linespacing=1.15, zorder=1000)
    x2, y2, _ = proj3d.proj_transform(S, Sy, zs, ax.get_proj())
    ax.annotate("  TIM + heat spreader / sink", (x2, y2), xycoords="data", fontsize=5.5, color="#555555", ha="left", va="center", zorder=1000)
    sm = plt.cm.ScalarMappable(norm=tnorm, cmap=CMAP_T); sm.set_array([])
    cb = fig.colorbar(sm, cax=fig.add_axes([0.845, 0.56, 0.018, 0.34])); cb.set_label("die temperature (K)", fontsize=7); cb.ax.tick_params(labelsize=6.5)
    cax = fig.add_axes([0.845, 0.10, 0.018, 0.34])
    if mode == "tnuc":
        sm = plt.cm.ScalarMappable(norm=tn_norm, cmap=CMAP_TNUC); sm.set_array([]); cb = fig.colorbar(sm, cax=cax)
        cb.set_label("rail nucleation time (years)", fontsize=7); cb.set_ticks([0.1, 1, 10, 100]); cb.set_ticklabels(["0.1", "1", "10", "100"])
        cb.ax.axhline(10, color="white", lw=1.2); cb.ax.axhline(10, color=VERMILION, lw=0.7)
        ttl = f"{D['case']}, {D['label']}: die temperature and rail nucleation time"
    else:
        sm = plt.cm.ScalarMappable(norm=snorm, cmap=SEQ); sm.set_array([]); cb = fig.colorbar(sm, cax=cax)
        cb.set_label("segment stress (MPa); red = critical", fontsize=7)
        ttl = f"{D['case']}, {D['label']}: die temperature and rail stress after {yr(D['hz'][q])}"
    cb.ax.tick_params(labelsize=6.5); fig.text(0.02, 0.97, ttl, fontsize=8.5, color=INK, ha="left", va="top")
    save(fig, d, name or ("fig_exploded_tnuc" if mode == "tnuc" else f"fig_exploded_stress_{D['hz'][q]:g}yr"), made)


# =====================================================================================================================
# flat maps of the power-delivery network of every die
# =====================================================================================================================
def _die_axes(nrow, ncol, w=3.4, h=3.5):
    fig, axs = plt.subplots(nrow, ncol, figsize=(w * ncol + 0.9, h * nrow), squeeze=False)
    return fig, axs


def _finish_map(ax, title):
    ax.set_title(title, loc="left"); ax.autoscale(); ax.set_aspect("equal", adjustable="datalim")


def pdn_stress(D, d, made):
    """(a) all dies x all horizons on one sheet; (b) one large sheet per die: horizontal layer, vertical layer, both."""
    style(); nd = len(D["dies"]); nh = len(D["hz"]); n = len(D["rail"])
    lo = min(stress_limits(D, q)[0] for q in range(nh)); lims = [(lo, D["sc"] / 1e6)] * nh       # one colour scale for all horizons
    fig, axs = _die_axes(nd, nh); last = None
    for k in range(nd):
        for q in range(nh):
            ax = axs[nd - 1 - k][q]; lc = _draw_stress(ax, D, [i for i in range(n) if D["die"][i] == k], q, lims[q])
            nf = int(np.sum((D["die"] == k) & (D["tn"] <= D["hz"][q] * Y)))
            _finish_map(ax, f"{D['dies'][k]} after {yr(D['hz'][q])}: {nf} of {int(np.sum(D['die'] == k))} rails nucleated")
            if k == 0: ax.set_xlabel("x (mm)")
            if q == 0: ax.set_ylabel("y (mm)")
            last = lc if lc is not None else last
    if last is not None:
        cb = fig.colorbar(last, ax=list(axs.ravel()), fraction=0.02, pad=0.02); cb.set_label("larger end-node stress of the segment (MPa)"); cb.outline.set_visible(False)
    fig.suptitle(f"{D['case']}, {D['label']}: stress on the power-delivery rails, every die (rows, top die first) at every horizon (columns)\n"
                 f"Red: an end of the segment is at or above the critical {D['sc']/1e6:.0f} MPa", x=0.01, y=0.93, ha="left", va="bottom", color=INK, fontsize=10)
    save(fig, d, "fig_pdn_stress_all", made)
    for k in range(nd):
        for q in range(nh):
            fig, axs = _die_axes(1, 3, 4.2, 4.3)
            sets = [("horizontal rails", [i for i in range(n) if D["die"][i] == k and D["kind"][i] == 0]),
                    ("vertical rails", [i for i in range(n) if D["die"][i] == k and D["kind"][i] == 1]),
                    ("both layers", [i for i in range(n) if D["die"][i] == k])]
            lc = None
            for ax, (ttl, ii) in zip(axs[0], sets):
                l2 = _draw_stress(ax, D, ii, q, lims[q], lw=1.5); lc = l2 if l2 is not None else lc
                _finish_map(ax, f"{ttl} ({len(ii)})"); ax.set_xlabel("x (mm)")
            axs[0][0].set_ylabel("y (mm)")
            if lc is not None:
                cb = fig.colorbar(lc, ax=list(axs[0]), fraction=0.03, pad=0.02); cb.set_label("larger end-node stress of the segment (MPa)"); cb.outline.set_visible(False)
            fig.suptitle(f"{D['case']}, {D['label']}: {D['dies'][k]} rail stress after {yr(D['hz'][q])}. Red: at or above {D['sc']/1e6:.0f} MPa",
                         x=0.01, y=1.03, ha="left", color=INK, fontsize=10)
            save(fig, d, f"fig_pdn_stress_{D['dies'][k]}_{D['hz'][q]:g}yr", made)


def _draw_stress(ax, D, ii, q, lim, lw=1.1):
    sg, v = [], []
    for i in ii:
        a, b = seg_values(D, i, q); sg.append(a); v.append(b)
    if not sg:
        ax.set_axis_off(); return None
    sg = np.concatenate(sg); v = np.concatenate(v); o = np.argsort(v)
    lc = LineCollection(sg[o], array=v[o], cmap=SEQ, linewidths=lw, capstyle="butt"); lc.set_clim(*lim); ax.add_collection(lc)
    if (v >= D["sc"] / 1e6).any():
        ax.add_collection(LineCollection(sg[v >= D["sc"] / 1e6], colors=CRIT, linewidths=lw + 0.7, capstyle="butt", zorder=3))
    return lc


def pdn_tnuc(D, d, made):
    style(); nd = len(D["dies"]); n = len(D["rail"]); fig, axs = _die_axes(1, nd); norm = LogNorm(0.03, 100.0); cm = plt.get_cmap(CMAP_TNUC)
    for k, ax in enumerate(axs[0]):
        ii = [i for i in range(n) if D["die"][i] == k]
        imm = [D["nodes"][i][[0, -1]] * 1e3 for i in ii if not np.isfinite(D["tn"][i])]
        mor = sorted([i for i in ii if np.isfinite(D["tn"][i])], key=lambda i: -D["tn"][i])
        if imm: ax.add_collection(LineCollection(imm, colors="#c3c2b7", linewidths=0.5))
        if mor:
            ax.add_collection(LineCollection([D["nodes"][i][[0, -1]] * 1e3 for i in mor], colors=[cm(norm(min(max(D["tn"][i] / Y, 0.03), 100))) for i in mor], linewidths=1.6))
            i0 = mor[-1]; P = D["nodes"][i0][[0, -1]] * 1e3; ax.plot(P[:, 0], P[:, 1], color=VERMILION, lw=3.0, solid_capstyle="round", zorder=4); ax.plot(P[:, 0], P[:, 1], color="white", lw=0.9, zorder=5)
        tn = np.array([D["tn"][i] for i in ii]); e = f"earliest {tn.min()/Y:.2f} yr" if len(tn) and np.isfinite(tn).any() else "none nucleates"
        _finish_map(ax, f"{D['dies'][k]}: {e}\n{int(np.sum(tn <= 10 * Y))} rails within 10 yr, {len(mor)} of {len(ii)} within 100 yr"); ax.set_xlabel("x (mm)")
    axs[0][0].set_ylabel("y (mm)"); sm = plt.cm.ScalarMappable(norm=norm, cmap=CMAP_TNUC); sm.set_array([])
    cb = fig.colorbar(sm, ax=list(axs[0]), fraction=0.03, pad=0.02); cb.set_label("rail nucleation time (years)"); cb.set_ticks([0.1, 1, 10, 100]); cb.set_ticklabels(["0.1", "1", "10", "100"]); cb.outline.set_visible(False)
    fig.suptitle(f"{D['case']}, {D['label']}: nucleation time of every rail (grey: never within 100 years; outlined: the earliest rail of the die)", x=0.01, y=1.08, ha="left", color=INK, fontsize=10)
    save(fig, d, "fig_pdn_tnuc", made)


def pdn_temperature(D, d, made):
    if all(t is None for t in D["Tn"]):
        return
    style(); nd = len(D["dies"]); n = len(D["rail"]); fig, axs = _die_axes(1, nd)
    allT = np.concatenate([t for t in D["Tn"] if t is not None]); norm = Normalize(float(allT.min()), float(allT.max())); lc = None
    for k, ax in enumerate(axs[0]):
        sg, v = [], []
        for i in range(n):
            if D["die"][i] == k and D["Tn"][i] is not None:
                P = D["nodes"][i] * 1e3; sg.append(np.stack([P[:-1], P[1:]], 1)); v.append(0.5 * (D["Tn"][i][:-1] + D["Tn"][i][1:]))
        if not sg:
            ax.set_axis_off(); continue
        v = np.concatenate(v); lc = LineCollection(np.concatenate(sg), array=v, cmap=CMAP_T, norm=norm, linewidths=1.2); ax.add_collection(lc)
        _finish_map(ax, f"{D['dies'][k]}: {v.min():.1f} to {v.max():.1f} K on the rails"); ax.set_xlabel("x (mm)")
    axs[0][0].set_ylabel("y (mm)")
    if lc is not None:
        cb = fig.colorbar(lc, ax=list(axs[0]), fraction=0.03, pad=0.02); cb.set_label("rail temperature at the nodes (K)"); cb.outline.set_visible(False)
    fig.suptitle(f"{D['case']}, {D['label']}: temperature seen by the rails ({D.get('T_source', 'input of the stress solve')})", x=0.01, y=1.04, ha="left", color=INK, fontsize=10)
    save(fig, d, "fig_pdn_temperature", made)


# =====================================================================================================================
# curves
# =====================================================================================================================
def smax_history(D, d, made):
    style(); nd = len(D["dies"]); fig, axs = plt.subplots(1, nd, figsize=(3.5 * nd, 3.3), squeeze=False, sharey=True); t = D["times"] / Y; m = t > 0
    for k, ax in enumerate(axs[0]):
        ii = np.where(D["die"] == k)[0]; fin = np.isfinite(D["tn"][ii])
        for sel, col, lw, al in ((ii[~fin], "#c3c2b7", 0.5, 0.8), (ii[fin], SERIES[0], 0.6, 0.55)):
            if len(sel): ax.add_collection(LineCollection([np.column_stack([t[m], D["smax"][i][m] / 1e6]) for i in sel], colors=col, linewidths=lw, alpha=al))
        ax.axhline(D["sc"] / 1e6, color=CRIT, lw=0.9, ls=(0, (4, 3))); ax.axvline(10, color=MUTED, lw=0.7, ls=":")
        ax.set_xscale("log"); ax.set_xlim(t[m].min(), t[m].max()); ax.autoscale(axis="y"); ax.grid(True, color=GRID, lw=0.5)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g")); ax.set_xlabel("time (years)")
        ax.set_title(f"{D['dies'][k]}: {int(fin.sum())} of {len(ii)} rails cross the critical stress", loc="left")
    axs[0][0].set_ylabel("largest node stress of the rail (MPa)")
    axs[0][0].annotate("critical stress", (0.02, D["sc"] / 1e6), xycoords=("axes fraction", "data"), xytext=(0, 4), textcoords="offset points", color=INK2, fontsize=8)
    fig.suptitle(f"{D['case']}, {D['label']}: stress history of every rail (blue: nucleates within 100 years; grey: does not; dotted: 10 years)", x=0.01, y=1.04, ha="left", color=INK, fontsize=10)
    save(fig, d, "fig_smax_history", made)


def tnuc_cdf(D, d, made):
    style(); nd = len(D["dies"]); fig, axs = plt.subplots(1, nd, figsize=(3.5 * nd, 3.3), squeeze=False, sharey=False)
    series = [(D["label"], D["tn"])] + [(k, v) for k, v in D.get("ref", {}).items()]
    for k, ax in enumerate(axs[0]):
        mk = D["die"] == k
        for q, (name, tn) in enumerate(series):
            t = np.sort(tn[mk][np.isfinite(tn[mk])]) / Y
            if len(t): ax.step(np.concatenate([[t[0]], t]), np.arange(0, len(t) + 1), where="post", color=SERIES[q], lw=2.0 if q == 0 else 1.3, ls="-" if q == 0 else (0, (4, 2)) if q == 1 else ":", label=name)
        ax.axvline(10, color=MUTED, lw=0.7, ls=":"); ax.set_xscale("log"); ax.set_xlim(0.5, 100); ax.grid(True, color=GRID, lw=0.5)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g")); ax.set_xlabel("time (years)"); ax.set_title(f"{D['dies'][k]} ({int(mk.sum())} rails)", loc="left")
    axs[0][0].set_ylabel("rails that have nucleated"); axs[0][0].legend(loc="best", fontsize=8)
    fig.suptitle(f"{D['case']}: number of nucleated rails against time, per die (dotted: 10 years)", x=0.01, y=1.04, ha="left", color=INK, fontsize=10)
    save(fig, d, "fig_tnuc_cdf", made)


def stress_hist(D, d, made):
    style(); nd = len(D["dies"]); nh = len(D["hz"]); fig, axs = plt.subplots(1, nd, figsize=(3.5 * nd, 3.2), squeeze=False, sharey=False)
    ih = [int(np.argmin(np.abs(D["times"] - h * Y))) for h in D["hz"]]; allv = D["smax"][:, ih] / 1e6; bins = np.linspace(allv.min(), allv.max(), 31)
    for k, ax in enumerate(axs[0]):
        for q in range(nh):
            ax.hist(D["smax"][D["die"] == k][:, ih[q]] / 1e6, bins=bins, histtype="step", lw=1.6, color=SERIES[q], label=f"after {yr(D['hz'][q])}")
        ax.axvline(D["sc"] / 1e6, color=CRIT, lw=0.9, ls=(0, (4, 3))); ax.grid(True, axis="y", color=GRID, lw=0.5); ax.set_xlabel("largest node stress of the rail (MPa)")
        ax.set_title(D["dies"][k], loc="left")
    axs[0][0].set_ylabel("rails"); axs[0][-1].legend(loc="upper right")
    fig.suptitle(f"{D['case']}, {D['label']}: distribution of the rail stress per die (dashed red: critical stress)", x=0.01, y=1.04, ha="left", color=INK, fontsize=10)
    save(fig, d, "fig_stress_hist", made)


def rail_profiles(D, d, made, n_per_die=3):
    """Stress along the n earliest rails of every die at every horizon (nodes; the in-segment field when available)."""
    style(); nh = len(D["hz"])
    for k, dname in enumerate(D["dies"]):
        ii = np.where(D["die"] == k)[0]; order = ii[np.argsort(D["tn"][ii])][:n_per_die]
        order = [i for i in order if np.isfinite(D["tn"][i])] or list(ii[np.argsort(-D["smax"][ii].max(1))][:n_per_die])
        if not order:
            continue
        fig, axs = plt.subplots(len(order), nh, figsize=(3.6 * nh, 2.5 * len(order)), squeeze=False, sharey=True)
        for r, i in enumerate(order):
            P = D["nodes"][i]; s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P, axis=0).T))]) * 1e6
            for q in range(nh):
                ax = axs[r][q]; fld = D.get("field", {}).get(int(D["rail"][i]))
                if fld is not None: ax.plot(fld[0] * 1e6, fld[1][:, q] / 1e6, color=SERIES[1], lw=1.6, zorder=2)
                ax.plot(s, D["sig"][i][:, q] / 1e6, ls="-" if fld is None else "none", lw=0.8, marker="o", ms=3, color=SERIES[0], zorder=3)
                ax.axhline(D["sc"] / 1e6, color=CRIT, lw=0.9, ls=(0, (4, 3))); ax.grid(True, color=GRID, lw=0.5)
                tt = f"{D['tn'][i]/Y:.2f} yr" if np.isfinite(D["tn"][i]) else "no nucleation"
                ax.set_title(f"rail {int(D['rail'][i])} ({'V' if D['kind'][i] else 'H'}), {tt}; after {yr(D['hz'][q])}", loc="left", fontsize=8)
                if r == len(order) - 1: ax.set_xlabel("position along the rail (um)")
            axs[r][0].set_ylabel("stress (MPa)")
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle(f"{D['case']}, {D['label']}: stress along the earliest rails of {dname} (dots: nodes; line: in-segment field where exported; dashed red: critical)",
                     x=0.01, y=0.99, ha="left", color=INK, fontsize=10)
        save(fig, d, f"fig_rail_profiles_{dname}", made)


def tnuc_parity(D, d, made):
    if not D.get("ref"):
        return
    style(); fig, ax = plt.subplots(figsize=(4.3, 4.1)); lo, hi = np.inf, 0
    for q, (name, tn) in enumerate(D["ref"].items()):
        m = np.isfinite(tn) & np.isfinite(D["tn"])
        if m.any():
            ax.scatter(D["tn"][m] / Y, tn[m] / Y, s=14, marker="so"[q % 2], color=SERIES[q + 1], edgecolor=SURF, linewidths=0.5, label=f"{name} ({m.sum()} rails)", zorder=3)
            lo = min(lo, D["tn"][m].min() / Y, tn[m].min() / Y); hi = max(hi, D["tn"][m].max() / Y, tn[m].max() / Y)
    if hi == 0:
        plt.close(fig); return
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8, zorder=1); ax.set_xscale("log"); ax.set_yscale("log")
    tk = [t for t in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100) if lo * 0.8 <= t <= hi * 1.25]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.FixedLocator(tk)); axis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g")); axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid(True, color=GRID, lw=0.5); ax.set_xlabel(f"{D['label']} nucleation time (years)"); ax.set_ylabel("nucleation time (years)")
    ax.set_title(f"{D['case']}: nucleation time of every rail that nucleates in both", loc="left"); ax.legend(loc="upper left")
    save(fig, d, "fig_tnuc_parity_all", made)


def all_figures(D, d):
    """Every figure of a run.  Returns the list of files written."""
    made = []; os.makedirs(d, exist_ok=True)
    exploded(D, d, made, "tnuc")
    for q in range(len(D["hz"])):
        exploded(D, d, made, "stress", q)
    pdn_stress(D, d, made); pdn_tnuc(D, d, made); pdn_temperature(D, d, made)
    smax_history(D, d, made); tnuc_cdf(D, d, made); stress_hist(D, d, made); rail_profiles(D, d, made); tnuc_parity(D, d, made)
    return made


# =====================================================================================================================
# COMSOL | StackEM | difference, node stress at one horizon, one row per die
# =====================================================================================================================
CMAP_SIG, CMAP_DIFF = "RdBu_r", "BrBG"


def _node_segments(D, i, vals):
    """Every segment of rail i split into 8 pieces whose colour runs linearly from one node value to the next (mm, MPa)."""
    P = D["nodes"][i] * 1e3; n = len(P) - 1; f = np.linspace(0, 1, 9)
    a = P[:-1, None, :] + (P[1:] - P[:-1])[:, None, :] * f[None, :, None]            # (n, 9, 2)
    v = vals[:-1, None] + (vals[1:] - vals[:-1])[:, None] * (f[:-1] + 0.0625)[None, :]
    return np.stack([a[:, :-1], a[:, 1:]], 2).reshape(-1, 2, 2), v.reshape(-1) / 1e6


def stress_compare(D, d, made, other, other_label="StackEM", q=None, name=None):
    """other: list of (S+1,) node-stress arrays (Pa) at the horizons, NaN where the other method gives no transient
    (StackEM: rails the Blech screen called immortal).  Columns: D['label'] | other | difference."""
    style(); nd = len(D["dies"]); n = len(D["rail"]); qs = range(len(D["hz"])) if q is None else [q]
    for q in qs:
        fig, axs = _die_axes(nd, 3, 3.4, 3.4); h = D["hz"][q]
        allv = np.concatenate([D["sig"][i][:, q] for i in range(n)]) / 1e6; lim = float(np.nanpercentile(np.abs(allv), 99.5))
        dif = np.concatenate([other[i][:, q] - D["sig"][i][:, q] for i in range(n) if other[i] is not None and not np.isnan(other[i][:, q]).all()]) / 1e6
        dlim = max(float(np.nanpercentile(np.abs(dif), 99.0)), 0.5) if len(dif) else 1.0
        lcs = [None, None, None]
        for k in range(nd):
            ii = [i for i in range(n) if D["die"][i] == k]; row = axs[nd - 1 - k]
            for col, (ax, src) in enumerate(zip(row, ("ref", "other", "diff"))):
                sg, v, grey = [], [], []
                for i in ii:
                    o = other[i][:, q] if other[i] is not None else None; miss = o is None or np.isnan(o).all()
                    if src == "ref": a, b = _node_segments(D, i, D["sig"][i][:, q])
                    elif miss: grey.append(D["nodes"][i][[0, -1]] * 1e3); continue
                    elif src == "other": a, b = _node_segments(D, i, o)
                    else: a, b = _node_segments(D, i, o - D["sig"][i][:, q])
                    sg.append(a); v.append(b)
                if grey: ax.add_collection(LineCollection(grey, colors="#c3c2b7", linewidths=0.6))
                if sg:
                    sg = np.concatenate(sg); v = np.concatenate(v)
                    lc = LineCollection(sg, array=v, cmap=CMAP_DIFF if src == "diff" else CMAP_SIG, linewidths=1.2, capstyle="butt")
                    lc.set_clim(-dlim, dlim) if src == "diff" else lc.set_clim(-lim, lim); ax.add_collection(lc); lcs[col] = lc
                nm = int(sum(1 for i in ii if other[i] is not None and not np.isnan(other[i][:, q]).all()))
                ttl = {"ref": f"{D['dies'][k]}: {D['label']}", "other": f"{D['dies'][k]}: {other_label} ({nm} of {len(ii)} rails closed)",
                       "diff": f"{D['dies'][k]}: {other_label} minus {D['label']}"}[src]
                if src == "diff" and sg is not None and len(v): ttl += f"; largest {np.abs(v).max():.1f} MPa"
                _finish_map(ax, ttl)
                if k == 0: ax.set_xlabel("x (mm)")
            row[0].set_ylabel("y (mm)")
        for col, lab in ((0, "hydrostatic stress at the nodes (MPa)"), (1, "hydrostatic stress at the nodes (MPa)"), (2, f"difference (MPa); grey: rails {other_label} did not close")):
            if lcs[col] is not None:
                cb = fig.colorbar(lcs[col], ax=list(axs[:, col]), fraction=0.025, pad=0.1, location="bottom"); cb.set_label(lab); cb.outline.set_visible(False)
        fig.suptitle(f"{D['case']}: node stress after {yr(h)}, {D['label']} against {other_label}. Each segment is shaded linearly between its two node values; "
                     f"red = tensile, blue = compressive.", x=0.01, y=0.995, ha="left", va="bottom", color=INK, fontsize=10)
        save(fig, d, name or f"fig_stress_compare_{h:g}yr", made)
