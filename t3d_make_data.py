"""
Builds data/<case>/ for the 3-D COMSOL model from a finished StackEM run (Linux / WSL, the StackEM environment).

    python t3d_make_data.py --case quad4 --repo ~/StackEM_clean --root ~/stackem_work/outputs --hotspot ~/stackem_work/HotSpot \
                            --weights ~/stackem_work/outputs/base3/skn_distill/skn_student.npz --ref-workers 16

What it needs (written by the StackEM experiments):
    base3 / quad4 / dense3   <root>/<case>/e2/e2_full_truth_rails.json      (stackem.experiments.e2_rail_accuracy)
                             <root>/<case>/e2/hotspot/                       stack.lcf, *.flp, stack.ptrace, stack.config, stack.grid.steady
    ev6                      <root>/base3/e8/e8_full_rails.json             (stackem.experiments.e8_ev6)
                             <root>/base3/e8/hotspot_full/                   ev6_3D.lcf, *.flp, ev6_3D.ptrace, ev6.config, ev6.grid.steady
What it writes (data/<case>/):
    rails.csv, segments.csv, constants.json, rails.json      every rail of the case      (comsol_full/cf_export.py of the repo)
    stackem_result.csv, stackem_timing.json                  reference solver + StackEM nucleation time of every rail
                                                             (comsol_full/cf_bench_stackem.py, numpy engine)
    stack.json                                               every layer of the HotSpot stack bottom -> top (name, kind, thickness,
                                                             conductivity, footprint), the inserts of a layer (ev6: TSV strips in the
                                                             TIM), spreader, sink, package and the mesh defaults of the 3-D model
    power_blocks.csv                                         die_id, x0, x1, y0, y1, P_W  of every powered floorplan block
    hotspot_T.csv                                            die_id, x, y, T  of every HotSpot grid cell of every die
Everything is read from the files HotSpot actually ran on, nothing is re-derived from the config.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# mesh defaults of the 3-D model per case: lateral element size in the die stack / spreader / sink (m) and the number of
# element layers through [die, bond, TIM, spreader, sink].  'thermal' = comsol_3d (temperature is solved on this mesh),
# 'pictures' = comsol_3d_hotspot (nothing is solved in the solid; the mesh only carries the pictures).
MESH = {
    "base3":  dict(thermal=dict(hdie=50e-6,  hspr=400e-6, hsink=1.5e-3, nz=[4, 2, 2, 5, 8]), pictures=dict(hdie=200e-6, hspr=1e-3,   hsink=3e-3, nz=[2, 1, 1, 2, 2])),
    "dense3": dict(thermal=dict(hdie=50e-6,  hspr=400e-6, hsink=1.5e-3, nz=[4, 2, 2, 5, 8]), pictures=dict(hdie=200e-6, hspr=1e-3,   hsink=3e-3, nz=[2, 1, 1, 2, 2])),
    "quad4":  dict(thermal=dict(hdie=75e-6,  hspr=500e-6, hsink=1.5e-3, nz=[4, 2, 2, 5, 8]), pictures=dict(hdie=300e-6, hspr=1.5e-3, hsink=3e-3, nz=[2, 1, 1, 2, 2])),
    "ev6":    dict(thermal=dict(hdie=100e-6, hspr=600e-6, hsink=2e-3,   nz=[4, 2, 2, 5, 8]), pictures=dict(hdie=400e-6, hspr=2e-3,   hsink=4e-3, nz=[2, 1, 1, 2, 2])),
}
EV6_DIE_NAMES = {"ev6_3D_cache_1.flp": "cache1", "ev6_3D_cache_2.flp": "cache2", "ev6_3D_core_layer.flp": "cores"}   # stackem.experiments.e8_ev6.EV6_DIES


def read_lcf(path):
    lines = [l.strip() for l in open(path) if l.strip() and not l.strip().startswith("#")]
    if len(lines) % 7:
        raise SystemExit(f"{path}: {len(lines)} non-comment lines, expected a multiple of 7")
    out = []
    for k in range(0, len(lines), 7):
        out.append(dict(index=int(lines[k]), lateral=lines[k + 1].upper() == "Y", power=lines[k + 2].upper() == "Y",
                        cp=float(lines[k + 3]), res=float(lines[k + 4]), thick=float(lines[k + 5]), flp=lines[k + 6]))
    out.sort(key=lambda L: L["index"])
    if [L["index"] for L in out] != list(range(len(out))):
        raise SystemExit(f"{path}: layer indices are not 0..{len(out)-1}")
    return out


def read_flp(path):
    """(name, w, h, x, y, cp or None, res or None) of every block."""
    out = []
    for l in open(path):
        s = l.strip()
        if not s or s.startswith("#"):
            continue
        v = s.split()
        if len(v) < 5:
            raise SystemExit(f"{path}: cannot read block line '{s}'")
        out.append((v[0], float(v[1]), float(v[2]), float(v[3]), float(v[4]), float(v[5]) if len(v) > 5 else None, float(v[6]) if len(v) > 6 else None))
    return out


def read_ptrace(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    names = lines[0].split(); vals = np.array([[float(v) for v in l.split()] for l in lines[1:]], float)
    if vals.shape[1] != len(names):
        raise SystemExit(f"{path}: {len(names)} names but {vals.shape[1]} values per row")
    return dict(zip(names, vals.mean(axis=0)))


def read_config(path):
    out = {}
    for l in open(path):
        m = re.match(r"\s*-(\w+)\s+(\S+)", l)
        if m:
            try: out[m.group(1)] = float(m.group(2))
            except ValueError: out[m.group(1)] = m.group(2)
    return out


def read_grid_steady(path, rows, cols):
    layers = {}; cur = None
    for line in open(path):
        m = re.match(r"Layer (\d+):", line)
        if m:
            cur = int(m.group(1)); layers[cur] = np.full(rows * cols, np.nan); continue
        if cur is not None and line.strip():
            k, v = line.split(); layers[cur][int(k)] = float(v)
    for k, v in layers.items():
        if np.isnan(v).any():
            raise SystemExit(f"{path}: layer {k} has {int(np.isnan(v).sum())} missing cells")
    return {k: v.reshape(rows, cols) for k, v in layers.items()}       # row 0 = top edge (y = size_y)


def run(cmd, cwd, env):
    print("$ " + " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=cwd, env=env)
    if r.returncode:
        raise SystemExit(f"command failed ({r.returncode}): {' '.join(cmd)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, choices=["base3", "quad4", "dense3", "ev6"])
    ap.add_argument("--repo", required=True, help="the StackEM_clean checkout (its comsol_full/ scripts and the stackem package are used)")
    ap.add_argument("--root", required=True, help="StackEM output root (paths.root_out)")
    ap.add_argument("--hotspot", required=True, help="HotSpot directory (paths.hotspot_dir)")
    ap.add_argument("--weights", required=True, help="skn_student.npz")
    ap.add_argument("--out", default=None, help="output folder (default <this folder>/data/<case>)")
    ap.add_argument("--ref-workers", type=int, default=max(1, (os.cpu_count() or 1)), help="processes for the reference solver")
    ap.add_argument("--skip-export", action="store_true", help="rails.csv etc. already there")
    ap.add_argument("--skip-bench", action="store_true", help="stackem_result.csv already there")
    a = ap.parse_args()
    repo = os.path.abspath(os.path.expanduser(a.repo)); root = os.path.abspath(os.path.expanduser(a.root))
    out = os.path.abspath(a.out or os.path.join(HERE, "data", a.case)); os.makedirs(out, exist_ok=True)
    env = dict(os.environ); env["PYTHONPATH"] = repo + os.pathsep + os.path.join(repo, "comsol_full") + os.pathsep + env.get("PYTHONPATH", "")
    env.setdefault("MPLBACKEND", "Agg")
    cf = os.path.join(repo, "comsol_full")
    for fn in ("cf_export.py", "cf_bench_stackem.py", "cf_common.py"):
        if not os.path.exists(os.path.join(cf, fn)):
            raise SystemExit(f"{os.path.join(cf, fn)} not found: --repo must be the StackEM_clean checkout that carries comsol_full/")
    # ---- 1. every rail: rails.csv, segments.csv, constants.json (cf_export writes to <out_root>/<case>) ----------------
    cfg = os.path.join(repo, "configs", f"{a.case}.json" if a.case in ("quad4", "dense3") else "base3.json")
    if not a.skip_export:
        run([sys.executable, os.path.join(cf, "cf_export.py"), "--case", a.case, "--config", cfg, "--root", root, "--hotspot", os.path.expanduser(a.hotspot),
             "--out", os.path.dirname(out), "--horizons", "1,2,5,10"], repo, env)
        if os.path.basename(out) != a.case:
            raise SystemExit("--out must end with the case name (cf_export writes <out parent>/<case>)")
    # ---- 2. reference solver + StackEM nucleation time of every rail -------------------------------------------------
    if not a.skip_bench:
        run([sys.executable, os.path.join(cf, "cf_bench_stackem.py"), "--case-dir", out, "--weights", os.path.expanduser(a.weights),
             "--engine", "numpy", "--ref-workers", str(a.ref_workers), "--repeat", "1", "--tag", "numpy_for_3d"], repo, env)
    C = json.load(open(os.path.join(out, "constants.json"))); dies = list(C["dies"])
    # ---- 3. the HotSpot stack as it was run ---------------------------------------------------------------------------
    if a.case == "ev6":
        hs = os.path.join(root, "base3", "e8", "hotspot_full"); lcf = "ev6_3D.lcf"; ptr = "ev6_3D.ptrace"; cfgf = "ev6.config"; tag = "ev6"
    else:
        hs = os.path.join(root, a.case, "e2", "hotspot"); lcf = "stack.lcf"; ptr = "stack.ptrace"; cfgf = "stack.config"; tag = "stack"
    for fn in (lcf, ptr, cfgf, f"{tag}.grid.steady"):
        if not os.path.exists(os.path.join(hs, fn)):
            raise SystemExit(f"{os.path.join(hs, fn)} not found: run the StackEM experiment of this case first")
    layers = read_lcf(os.path.join(hs, lcf)); powers = read_ptrace(os.path.join(hs, ptr)); conf = read_config(os.path.join(hs, cfgf))
    rows, cols = int(conf["grid_rows"]), int(conf["grid_cols"])
    grids = read_grid_steady(os.path.join(hs, f"{tag}.grid.steady"), rows, cols)
    # die layers = powered layers; their names must be the die names of the rails
    die_layers = [L for L in layers if L["power"]]
    names = []
    for L in die_layers:
        nm = EV6_DIE_NAMES.get(L["flp"]) if a.case == "ev6" else os.path.splitext(L["flp"])[0]
        if nm is None:
            raise SystemExit(f"no die name for floorplan {L['flp']}")
        names.append(nm)
    if sorted(names) != sorted(dies):
        raise SystemExit(f"die names of the HotSpot stack {names} differ from the dies of the rails {dies}")
    if names != dies:
        # cf_export numbers the dies by name; the 3-D model numbers them bottom -> top: renumber die_id in rails.csv and
        # the dies list of constants.json (ev6: cf_export gives cores, cache1, cache2; the stack is cache1, cache2, cores)
        remap = {dies.index(n): k for k, n in enumerate(names)}
        lines = open(os.path.join(out, "rails.csv")).read().splitlines(); hdr = lines[0].split(","); jd = hdr.index("die_id")
        with open(os.path.join(out, "rails.csv"), "w") as f:
            f.write(lines[0] + "\n")
            for l in lines[1:]:
                v = l.split(","); v[jd] = str(remap[int(v[jd])]); f.write(",".join(v) + "\n")
        C["dies"] = list(names); C["dies_note"] = "numbered bottom -> top (the order of the HotSpot layer stack)"
        json.dump(C, open(os.path.join(out, "constants.json"), "w"), indent=1); dies = list(names)
        print(f"die ids renumbered bottom -> top: {names}")
    # footprint and conductivity of every layer, inserts (blocks with their own resistivity)
    stack_layers = []; inserts = []; blocks = []
    for L in layers:
        flp = read_flp(os.path.join(hs, L["flp"]))
        sx = round(max(b[3] + b[1] for b in flp), 12); sy = round(max(b[4] + b[2] for b in flp), 12)       # no 1e-18 slivers between layers
        if min(b[3] for b in flp) < -1e-12 or min(b[4] for b in flp) < -1e-12:
            raise SystemExit(f"{L['flp']}: a block starts at a negative coordinate")
        kind = "die" if L["power"] else ("bond" if "bond" in L["flp"].lower() else "tim")
        name = names[die_layers.index(L)] if L["power"] else (f"{kind}{L['index']}" if kind == "bond" else ("TIM" if L is layers[-1] else f"TIM{L['index']}"))
        rec = dict(name=name, kind=kind, lcf_index=L["index"], thickness=L["thick"], k=1.0 / L["res"], cp=L["cp"], size_x=sx, size_y=sy, flp=L["flp"])
        if L["power"]:
            rec["die_index"] = names.index(name)
        stack_layers.append(rec)
        for b in flp:
            if b[6] is not None and abs(1.0 / b[6] - 1.0 / L["res"]) > 1e-9:            # a block with its own material (ev6: TSV strips)
                inserts.append(dict(layer=name, name=b[0], x0=round(b[3], 12), x1=round(b[3] + b[1], 12), y0=round(b[4], 12), y1=round(b[4] + b[2], 12), k=1.0 / b[6], cp=b[5]))
            if L["power"]:
                P = float(powers.get(b[0], 0.0))
                if P > 0:
                    blocks.append((names.index(name), round(b[3], 12), round(b[3] + b[1], 12), round(b[4], 12), round(b[4] + b[2], 12), P))
    missing = [n for n in powers if powers[n] > 0 and not any(abs(bk[5] - powers[n]) < 1e-12 for bk in blocks)]
    if missing:
        raise SystemExit(f"powered blocks of the trace without a floorplan rectangle on a die layer: {missing[:5]}")
    size_x = max(L["size_x"] for L in stack_layers if L["kind"] == "die"); size_y = max(L["size_y"] for L in stack_layers if L["kind"] == "die")
    spreader = dict(name="spreader", kind="spreader", thickness=conf["t_spreader"], k=conf["k_spreader"], size_x=conf["s_spreader"], size_y=conf["s_spreader"])
    sink = dict(name="sink", kind="sink", thickness=conf["t_sink"], k=conf["k_sink"], size_x=conf["s_sink"], size_y=conf["s_sink"])
    if stack_layers[-1]["kind"] == "die":
        raise SystemExit("the last HotSpot layer is a die: expected a TIM between the top die and the spreader")
    for L in stack_layers:                                                          # the rails sit 2 and 4 um above every die, inside the next layer
        if L["kind"] == "die":
            above = stack_layers[stack_layers.index(L) + 1]
            if above["thickness"] < 5e-6:
                raise SystemExit(f"layer {above['name']} above die {L['name']} is only {above['thickness']*1e6:.1f} um thick: the rail levels (2 and 4 um) do not fit")
    K = dict(case=a.case, die_names=names, die_size_x=size_x, die_size_y=size_y, layers=stack_layers + [spreader, sink], inserts=inserts,
             r_convec=conf["r_convec"], ambient=conf["ambient"], T_zero_stress=C["em"]["T_zero_stress"], dalpha=C["em"]["dalpha"],
             hotspot=dict(rows=rows, cols=cols, workdir=hs, lcf=lcf, config=cfgf, grid_model=conf.get("model_type", "grid")),
             mesh=MESH[a.case],
             note="layers bottom -> top as in the HotSpot .lcf (layer 0 farthest from the sink); the spreader and the sink are centred "
                  "on the die stack; convection 1/(r_convec * sink area) on the top face of the sink; inserts = floorplan blocks with "
                  "their own conductivity inside a layer (ev6: the TSV strips of the TIM layers)")
    json.dump(K, open(os.path.join(out, "stack.json"), "w"), indent=1)
    with open(os.path.join(out, "power_blocks.csv"), "w") as f:
        f.write("die_id,x0_m,x1_m,y0_m,y1_m,P_W\n")
        for b in blocks:
            f.write(f"{b[0]},{b[1]:.9e},{b[2]:.9e},{b[3]:.9e},{b[4]:.9e},{b[5]:.9e}\n")
    # ---- 4. the HotSpot temperature of every die on its grid ---------------------------------------------------------
    xs = (np.arange(cols) + 0.5) / cols * size_x; ys = size_y - (np.arange(rows) + 0.5) / rows * size_y          # row 0 = top edge
    with open(os.path.join(out, "hotspot_T.csv"), "w") as f:
        f.write("die_id,x_m,y_m,T_K\n")
        for L in stack_layers:
            if L["kind"] != "die":
                continue
            G = grids[L["lcf_index"]]
            for r in range(rows):
                for c in range(cols):
                    f.write(f"{L['die_index']},{xs[c]:.9e},{ys[r]:.9e},{G[r, c]:.4f}\n")
    # ---- 5. checks against the rails -----------------------------------------------------------------------------------
    S = np.genfromtxt(os.path.join(out, "segments.csv"), delimiter=",", names=True); R = np.genfromtxt(os.path.join(out, "rails.csv"), delimiter=",", names=True)
    P_tot = {k: sum(b[5] for b in blocks if b[0] == k) for k in range(len(names))}
    print(f"{a.case}: {len(names)} dies {names}, die {size_x*1e3:.3f} x {size_y*1e3:.3f} mm, {len(stack_layers)} HotSpot layers + spreader + sink, "
          f"{len(inserts)} inserts, {len(blocks)} powered blocks, P = " + ", ".join(f"{names[k]} {P_tot[k]:.2f} W" for k in P_tot))
    print("layers bottom -> top: " + " | ".join(f"{L['name']} {L['thickness']*1e6:.0f} um k={L['k']:.3g}" for L in K["layers"]))
    worst = 0.0
    for L in stack_layers:
        if L["kind"] != "die":
            continue
        G = grids[L["lcf_index"]]; d = L["die_index"]; m = R["die_id"] == d
        x0 = R["x0_m"][m]; y0 = R["y0_m"][m]
        if x0.min() < -1e-9 or y0.min() < -1e-9 or x0.max() > size_x + 1e-9 or y0.max() > size_y + 1e-9:
            raise SystemExit(f"die {L['name']}: rails outside the die footprint")
        # bilinear HotSpot temperature at the first node of every rail against T_L of its first segment (StackEM's own look-up;
        # ev6 was resampled to a square grid by e8, so a small difference is expected there)
        seg0 = S[(S["k"] == 0)]; TL = {int(r): float(t) for r, t in zip(seg0["rail"], seg0["T_L_K"])}
        fx = np.clip(x0 / size_x * cols - 0.5, 0, cols - 1); fy = np.clip((size_y - y0) / size_y * rows - 0.5, 0, rows - 1)
        c0 = np.floor(fx).astype(int); r0 = np.floor(fy).astype(int); c1 = np.minimum(c0 + 1, cols - 1); r1 = np.minimum(r0 + 1, rows - 1); wx = fx - c0; wy = fy - r0
        Tb = (1 - wx) * (1 - wy) * G[r0, c0] + wx * (1 - wy) * G[r0, c1] + (1 - wx) * wy * G[r1, c0] + wx * wy * G[r1, c1]
        dT = np.array([Tb[i] - TL[int(r)] for i, r in enumerate(R["rail"][m])])
        worst = max(worst, float(np.abs(dT).max()))
        print(f"  {L['name']}: HotSpot {G.min():.2f}..{G.max():.2f} K; rails {int(m.sum())}; grid T at the first rail node vs T_L of the rails: largest |diff| {np.abs(dT).max():.3f} K")
    lim = 2.0 if a.case == "ev6" else 0.05
    if worst > lim:
        raise SystemExit(f"the HotSpot grid does not match the temperatures of the rails ({worst:.3f} K): wrong HotSpot folder?")
    print(f"wrote {out}: " + ", ".join(sorted(os.listdir(out))))


if __name__ == "__main__":
    main()
