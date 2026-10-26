#!/usr/bin/env python3
"""PEG timing benchmark (lean). Runs jobs one after another and measures pegSerial / pegMPI over

    N  truncation index      P  points of the grating profile (blazed, rounded corners, "custom" profile)
    T  threads per process   R  MPI ranks

Everything goes into ONE folder (default benchmark_results/timing_<host>_<date>/):
    peg_results/<job>/   PEG output + progress file of every job
    timing.csv           one row per run                  (columns: stage,prog,N,P,T,R,S,rep,ok,time)
    summary.csv          one row per configuration        (best time, time/point, speed-up, efficiency)
    tables/*.tex         booktabs tables, \\input{} them  |  summary.txt  short interpretation
    timing.pdf/png       plots (matplotlib, optional)     |  surface_NP.pdf with --plot3d

time = pegSerial: wall time of the process; pegMPI: "Run time" printed by rank 0.
--measureTiming is only printed to the terminal, so the script reads it from the program output: every solver call
prints a "Timing Profile" block; the blocks of all calls (and ranks) are summed per section (columns alloc, refr,
layers, integ, matrix, rayleigh, effs, solver = "Total (solver) time"). summary.csv gives them per sweep point.
S = R*rounds sweep points per job; speed-up = throughput / throughput of (pegSerial, T=1) at the same N,P.
Run from the project root:   python3 run/timing_benchmark.py [--quick] [--plot3d]
"""
import argparse, csv, datetime, math, os, re, socket, statistics, subprocess, sys, time


CATS = [("alloc", "Allocate Memory"), ("refr", "Look up refractive index"),
        ("layers", "Compute alpha, beta values and layers"), ("integ", "Numerically integrating trial solutions"),
        ("matrix", "Matrix operations"), ("rayleigh", "Computing Rayleigh coeffients B_n"),
        ("effs", "Compute and package efficiencies"), ("solver", "Total (solver) time")]


def parse_timing(text):
    """Sum of every section of all 'Timing Profile' blocks in the program output, and the number of blocks."""
    labels = {lab: key for key, lab in CATS}
    sums = {key: 0.0 for key, _ in CATS}
    for line in text.splitlines():
        m = re.match(r"\s+(.+?):\s+([0-9.eE+-]+)\s*$", line)
        if m and m.group(1) in labels:
            sums[labels[m.group(1)]] += float(m.group(2))
    return sums, text.count("Timing Profile:")


# ---------- geometry: blazed grating with rounded corners as point list -----------------------------
def rounded_blaze(period, alpha, beta, rounding=0.05, M=8192):
    ta, tb = math.tan(math.radians(alpha)), math.tan(math.radians(beta))
    h, xa = period * ta * tb / (ta + tb), tb / (ta + tb)
    ys = [h * (i / M) / xa if i / M <= xa else h * (1 - i / M) / (1 - xa) for i in range(M)]
    w = max(1, int(rounding * M / 3)) | 1
    for _ in range(3 if rounding > 0 else 0):                   # 3 moving averages -> smooth corners
        ext = ys[-(w // 2):] + ys + ys[:w // 2]
        pre = [0.0]
        for v in ext:
            pre.append(pre[-1] + v)
        ys = [(pre[i + w] - pre[i]) / w for i in range(M)]
    i0 = min(range(M), key=ys.__getitem__)                      # valley at x=0
    ys = ys[i0:] + ys[:i0]
    ys = [v - ys[0] for v in ys] + [0.0]
    return max(ys), ys


def geometry(depth, fine, P):
    """--gratingGeometry string: depth, then (x,y) for P points from (0,0) to (1,0)."""
    M, vals = len(fine) - 1, ["%.10g" % depth]
    for k in range(P):
        t = k / (P - 1) * M
        i = min(int(t), M - 1)
        y = 0.0 if k in (0, P - 1) else ((1 - (t - i)) * fine[i] + (t - i) * fine[i + 1]) / depth
        vals += ["%.10g" % (k / (P - 1)), "%.10g" % y]
    return ",".join(vals)


# ---------- job list --------------------------------------------------------------------------------
def pow2(n):
    k, out = 1, []
    while k <= n:
        out.append(k)
        k *= 2
    return out + ([n] if out[-1] != n else [])


def make_jobs(a):
    jobs, seen = [], set()

    def add(stage, prog, N, P, T, R):
        if (prog, N, P, T, R) not in seen and T * R <= a.cores and (prog == "serial" or not a.no_mpi):
            seen.add((prog, N, P, T, R))
            jobs.append((stage, prog, N, P, T, R))

    for T in a.T_list: add("threads", "serial", a.N0, a.P0, T, 1)
    for N in a.N_list: add("orders", "serial", N, a.P0, 1, 1)
    for P in a.P_list: add("points", "serial", a.N0, P, 1, 1)
    for R in a.R_list: add("ranks", "mpi", a.N0, a.P0, 1, R)
    for R in a.R_list:
        for T in a.T_list: add("hybrid", "mpi", a.N0, a.P0, T, R)
    for N in a.N_list:
        for P in a.P_list: add("grid", "serial", N, P, 1, 1)
    return jobs


# ---------- running ---------------------------------------------------------------------------------
def run_job(a, geo, outdir, job, rep):
    stage, prog, N, P, T, R = job
    S = R * a.rounds
    jid = "%s_N%d_P%d_T%d_R%d" % (prog, N, P, T, R)
    d = os.path.join(outdir, "peg_results", jid)
    os.makedirs(d, exist_ok=True)
    out, prgf = os.path.join(d, "rep%d_output.txt" % rep), os.path.join(d, "rep%d_progress.txt" % rep)
    cmd = ["--mode", "constantIncidence", "--incidenceAngle", "88", "--eV", "--min", "100", "--max", str(100 + S - 1),
           "--increment", "1", "--gratingType", "custom", "--gratingMaterial", "Au", "--gratingPeriod", "1",
           "--gratingGeometry", geo(P), "--N", str(N), "--threads", str(T), "--computeTE", "--measureTiming",
           "--outputFile", out, "--progressFile", prgf]
    if prog == "serial":
        cmd = [a.serial] + cmd
    else:
        cmd = a.mpi_launcher.format(ranks=R).split() + [a.mpi] + cmd
    env = dict(os.environ, OMP_NUM_THREADS=str(T), PEG_MATERIAL_DB=a.material_db)
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=a.max_job_seconds)
        text, rc = r.stdout, r.returncode
    except subprocess.TimeoutExpired:
        text, rc = "", -9
    wall = time.perf_counter() - t0
    m = re.search(r"Run time \(s\):\s*([0-9.eE+-]+)", text)
    t = float(m.group(1)) if (prog == "mpi" and m) else wall
    status = ""
    try:
        status = [l.strip() for l in open(prgf) if l.startswith("status=")][0]
    except (OSError, IndexError):
        pass
    ok = rc == 0 and status == "status=succeeded"
    sums, calls = parse_timing(text)
    if calls != S:                  # e.g. a launcher that does not forward the output of all ranks
        sums = {k: float("nan") for k in sums}
    return dict(stage=stage, prog=prog, N=N, P=P, T=T, R=R, S=S, rep=rep, ok=int(ok), time=t, calls=calls, **sums)


# ---------- analysis --------------------------------------------------------------------------------
def power_exponent(xs, ys):
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
    n = len(pts)
    mx, my = sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n
    return sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx) ** 2 for p in pts)


def analyse(outdir, rows, a):
    best = {}
    for r in rows:
        if r["ok"]:
            k = (r["prog"], r["N"], r["P"], r["T"], r["R"])
            if k not in best or r["time"] < best[k]["time"]:
                best[k] = r                                     # fastest repetition
    summ = []
    for (prog, N, P, T, R), r in sorted(best.items()):
        t, S = r["time"], r["S"]
        base = best.get(("serial", N, P, 1, 1))
        sp = (S / t) / (base["S"] / base["time"]) if base else float("nan")
        summ.append(dict(prog=prog, N=N, P=P, T=T, R=R, S=S, time=t, tpp=t / S, speedup=sp, efficiency=sp / (T * R),
                         **{k + "_pp": r[k] / S for k, _ in CATS}))
    with open(os.path.join(outdir, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0].keys()))
        w.writeheader()
        w.writerows(summ)

    os.makedirs(os.path.join(outdir, "tables"), exist_ok=True)

    def table(name, header, data):
        with open(os.path.join(outdir, "tables", name + ".tex"), "w") as f:
            f.write("\\begin{tabular}{%s}\n\\toprule\n%s \\\\\n\\midrule\n" % ("r" * len(header), " & ".join(header)))
            for row in data:
                f.write(" & ".join("%.3g" % v if isinstance(v, float) else str(v) for v in row) + " \\\\\n")
            f.write("\\bottomrule\n\\end{tabular}\n")

    sel = lambda **kw: [s for s in summ if all(s[k] == v for k, v in kw.items())]
    text = []
    o = sel(prog="serial", P=a.P0, T=1, R=1)
    if len(o) > 1:
        table("orders", ["$N$", "$t_\\mathrm{pt}$ / s", "solver / s", "integration / s", "matrix ops / s"],
              [(s["N"], s["tpp"], s["solver_pp"], s["integ_pp"], s["matrix_pp"]) for s in o])
        text.append("time per point ~ N^%.2f (N=%d..%d, P=%d); solver time ~ N^%.2f, integration ~ N^%.2f, matrix operations ~ N^%.2f" %
                    ((power_exponent([s["N"] for s in o], [s["tpp"] for s in o]), o[0]["N"], o[-1]["N"], a.P0) +
                     tuple(power_exponent([s["N"] for s in o], [s[c + "_pp"] for s in o]) for c in ("solver", "integ", "matrix"))))
        text.append("integration is %.0f %% of the solver time at N=%d" % (100 * o[-1]["integ_pp"] / o[-1]["solver_pp"], o[-1]["N"]))
    p = sel(prog="serial", N=a.N0, T=1, R=1)
    if len(p) > 1:
        table("points", ["$P$", "$t_\\mathrm{pt}$ / s"], [(s["P"], s["tpp"]) for s in p])
        text.append("P=%d -> P=%d changes the time per point by %+.0f %% (N=%d)" % (p[0]["P"], p[-1]["P"], 100 * (p[-1]["tpp"] / p[0]["tpp"] - 1), a.N0))
    t = sel(prog="serial", N=a.N0, P=a.P0, R=1)
    if len(t) > 1:
        table("threads", ["$T$", "$t_\\mathrm{pt}$ / s", "speed-up", "efficiency"], [(s["T"], s["tpp"], s["speedup"], s["efficiency"]) for s in t])
        text.append("threads: speed-up %.2f at T=%d (efficiency %.0f %%)" % (t[-1]["speedup"], t[-1]["T"], 100 * t[-1]["efficiency"]))
    r_ = sel(prog="mpi", N=a.N0, P=a.P0, T=1)
    if len(r_) > 1:
        table("ranks", ["$R$", "$t_\\mathrm{pt}$ / s", "speed-up", "efficiency"], [(s["R"], s["tpp"], s["speedup"], s["efficiency"]) for s in r_])
        text.append("ranks: speed-up %.2f at R=%d (efficiency %.0f %%)" % (r_[-1]["speedup"], r_[-1]["R"], 100 * r_[-1]["efficiency"]))
    pair = [(sel(prog="serial", N=a.N0, P=a.P0, T=1, R=1), sel(prog="mpi", N=a.N0, P=a.P0, T=1, R=1))]
    if pair[0][0] and pair[0][1]:
        text.append("pegMPI(R=1) / pegSerial: run time per point x%.2f, pure solver time per point x%.2f" %
                    (pair[0][1][0]["tpp"] / pair[0][0][0]["tpp"], pair[0][1][0]["solver_pp"] / pair[0][0][0]["solver_pp"]))
    grid = sel(prog="serial", T=1, R=1)
    if any(math.isnan(s["solver_pp"]) for s in summ):
        text.append("WARNING: the timing profile of some runs is incomplete (number of 'Timing Profile' blocks != sweep points); "
                    "check that the MPI launcher forwards the output of all ranks")
    bad = sum(1 for r in rows if not r["ok"])
    if bad:
        text.append("WARNING: %d run(s) failed or timed out and are excluded" % bad)
    open(os.path.join(outdir, "summary.txt"), "w").write("\n".join(text) + "\n")
    print("\n".join(text))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed: no plots (the CSV/TeX files are complete)")
        return
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))
    if len(o) > 1:
        ax[0][0].loglog([s["N"] for s in o], [s["tpp"] for s in o], "o-"); ax[0][0].set_xlabel("N"); ax[0][0].set_ylabel("time per point / s")
    if len(p) > 1:
        ax[0][1].semilogx([s["P"] for s in p], [s["tpp"] for s in p], "o-", base=2); ax[0][1].set_xlabel("P (profile points)")
    for axx, data, lab in ((ax[1][0], t, "threads T"), (ax[1][1], r_, "MPI ranks R")):
        if len(data) > 1:
            x = [s["T"] if lab[0] == "t" else s["R"] for s in data]
            axx.plot(x, [s["speedup"] for s in data], "o-", label="measured"); axx.plot(x, x, "k--", label="ideal")
            axx.set_xlabel(lab); axx.set_ylabel("speed-up"); axx.legend()
    for x in ax.flat: x.grid(True, alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "timing.pdf")); fig.savefig(os.path.join(outdir, "timing.png"), dpi=150)
    if a.plot3d and len({s["N"] for s in grid}) > 1 and len({s["P"] for s in grid}) > 1:
        fig = plt.figure(figsize=(6.5, 5)); ax3 = fig.add_subplot(111, projection="3d")
        xs, ys, zs = [math.log2(s["N"]) for s in grid], [math.log2(s["P"]) for s in grid], [math.log10(s["tpp"]) for s in grid]
        ax3.plot_trisurf(xs, ys, zs, cmap="viridis", alpha=.9); ax3.scatter(xs, ys, zs, c="k", s=8)
        ax3.set_xlabel("log2 N"); ax3.set_ylabel("log2 P"); ax3.set_zlabel("log10 time per point / s")
        fig.savefig(os.path.join(outdir, "surface_NP.pdf")); fig.savefig(os.path.join(outdir, "surface_NP.png"), dpi=150)


# ---------- main ------------------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here) if os.path.basename(here) in ("run", "tools") else os.getcwd()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--serial", default=os.path.join(root, "pegSerial")); ap.add_argument("--mpi", default=os.path.join(root, "pegMPI"))
    ap.add_argument("--no-mpi", action="store_true"); ap.add_argument("--mpi-launcher", default="mpiexec -n {ranks}")
    ap.add_argument("--material-db", default=os.environ.get("PEG_MATERIAL_DB", os.path.join(root, "materialDatabase")))
    ap.add_argument("--out"); ap.add_argument("--quick", action="store_true", help="small lists, 1 repeat")
    ap.add_argument("--cores", type=int, default=len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count())
    ap.add_argument("--N-list"); ap.add_argument("--P-list"); ap.add_argument("--T-list"); ap.add_argument("--R-list")
    ap.add_argument("--rounds", type=int, default=2, help="sweep points per rank and job")
    ap.add_argument("--repeats", type=int, default=None); ap.add_argument("--max-job-seconds", type=float, default=600)
    ap.add_argument("--plot3d", action="store_true")
    a = ap.parse_args()
    ints = lambda s, d: [int(x) for x in s.split(",")] if s else d
    a.N_list = ints(a.N_list, [6, 10, 16, 24] if a.quick else [6, 8, 12, 16, 24, 32])
    a.P_list = ints(a.P_list, [9, 65, 513] if a.quick else [9, 17, 33, 65, 129, 257, 513, 1025])
    a.T_list = ints(a.T_list, pow2(a.cores)); a.R_list = ints(a.R_list, pow2(a.cores))
    a.N0 = min(a.N_list, key=lambda n: abs(n - 12)); a.P0 = min(a.P_list, key=lambda p: abs(p - 65))
    a.repeats = a.repeats or (1 if a.quick else 3)
    if not os.path.exists(a.mpi) and not a.no_mpi:
        print("pegMPI not found -> only pegSerial"); a.no_mpi = True
    outdir = a.out or os.path.join(root, "benchmark_results", "timing_%s_%s" % (socket.gethostname(), datetime.datetime.now().strftime("%Y%m%d-%H%M")))
    os.makedirs(outdir, exist_ok=True)
    depth, fine = rounded_blaze(1.0, 2.5, 30.0)
    geo = lambda P, c={}: c.setdefault(P, geometry(depth, fine, P))

    jobs, rows, dead = make_jobs(a), [], set()
    with open(os.path.join(outdir, "timing.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["stage", "prog", "N", "P", "T", "R", "S", "rep", "ok", "time", "calls"] + [k for k, _ in CATS])
        w.writeheader()
        for i, job in enumerate(jobs, 1):
            if job[0] == "orders" and any(job[2] >= n for n in dead):
                continue                                          # a smaller N already hit the time limit
            for rep in range(1, a.repeats + 1):
                r = run_job(a, geo, outdir, job, rep)
                rows.append(r); w.writerow(r); f.flush()
                print("[%d/%d] %s N=%d P=%d T=%d R=%d rep %d: %s %.2f s" % (i, len(jobs), job[1], job[2], job[3], job[4], job[5], rep, "ok" if r["ok"] else "FAILED", r["time"]))
                if not r["ok"]:
                    dead.add(job[2]); break
    analyse(outdir, rows, a)
    print("results in", outdir)


if __name__ == "__main__":
    main()
