"""
AI-Generated scipt

PEG timing benchmark (interactive): how do pegSerial / pegMPI scale with different parameters?

The script asks for the range (from / to) of every parameter. Enter = default value, '-' = do not scan this parameter.
All other flags of pegSerial / pegMPI stay at their defaults (see run_job).

  N    truncation index                      --N
  P    points of the grating profile         --gratingGeometry   (values 2^k+1)
  T    threads per process                   --threads           (values 2^k)
  R    MPI ranks (= nodes)                   mpiexec -n R        (values 2^k)
  E    integration tolerance                 --integrationTolerance <eps>   (values 10^k, scanned from loose to tight)

SCANS: in a scan ONE parameter runs through its whole range, the other four stay at their fixed value (REF).
    serial:       N-, P-, T-, E-scan
    mpi, per R:   N-, P-, T-, E-scan (T>1 = hybrid)   + one R-scan (T=1)
Every point of a scan is run LOOPS times. S = R*rounds sweep points (energies) per job, so every rank has the same work.

OUTPUT  (default  <script folder>/benchmark_results/timing_<host>_<date>/)
Plain-text .xyz tables: '#' header (context + description of every column), then whitespace separated numbers
(nan = not available). Load with gnuplot, numpy.loadtxt, pandas.read_csv(sep=r"\\s+", comment="#").
    raw.xyz      every single run
    results.xyz  one row per point = simple average over the loops (nothing else is calculated)
    fails.xyz    points with a scatter of the wall time above 10 %, failed / timed out runs, skipped points
"""
import argparse, datetime, math, os, re, shlex, signal, socket, statistics, subprocess, time

# ============================== DEFAULTS (what the prompts offer when you just press Enter) ====================================================
LOOPS = 3                       # repetitions per point
ROUNDS = 2                      # sweep points per rank and job
MAX_REL_STD = 0.10              # relative standard deviation of the wall time above which a point is written to fails.xyz
REF = dict(N=16, P=65, T=1, R=1, E=1e-5)           # fixed values while another parameter is scanned
PARAMS = ["N", "P", "T", "R", "E"]
SPEC = {                        # label, type, default from, default to
    "N": ("trunc. index N", int, 8, 32),
    "P": ("profile points P   (2^k+1)", int, 9, 513),
    "T": ("threads T          (2^k)", int, 1, 16),
    "R": ("MPI ranks R        (2^k)", int, 1, 16),
    "E": ("integration tolerance eps   (10^k)", float, 1e-8, 1e-3),
}
# ==========================================================================================================================================

CATS = [("alloc", "Allocate Memory"), ("refr", "Look up refractive index"),
        ("layers", "Compute alpha, beta values and layers"), ("integ", "Numerically integrating trial solutions"),
        ("matrix", "Matrix operations"), ("rayleigh", "Computing Rayleigh coeffients B_n"),
        ("effs", "Compute and package efficiencies"), ("solver", "Total (solver) time")]
NAN = float("nan")


def parse_timing(text):
    """Sum of every section of all 'Timing Profile' blocks in the program output, and the number of blocks."""
    labels = {lab: key for key, lab in CATS}
    sums = {key: 0.0 for key, _ in CATS}
    for line in text.splitlines():
        m = re.match(r"\s*(.+?):\s+([0-9.eE+-]+)\s*(?:s|sec|seconds)?\s*$", line)
        if m and m.group(1) in labels:
            sums[labels[m.group(1)]] += float(m.group(2))
    return sums, text.count("Timing Profile:")


# ---------- output: plain-text tables with human-readable header -----------------------------------
class Table:
    """'#' header (title, context, one line per column), then whitespace separated numeric columns. Rows are flushed immediately."""
    def __init__(self, path, title, about, cols, meta):
        self.cols = [(n, k, d, max(len(n), 6 if k == "i" else 13)) for n, k, d in cols]
        self.f = open(path, "w")
        w = self.f.write
        w("# %s\n#\n" % title)
        for l in about:
            w("# %s\n" % l)
        w("#\n# --- run context ---\n")
        for l in meta:
            w("# %s\n" % l)
        w("#\n# --- columns (whitespace separated, numeric, nan = not available; times in seconds) ---\n")
        for i, (n, k, d, wd) in enumerate(self.cols, 1):
            w("#  %2d  %-12s %s\n" % (i, n, d))
        w("#\n# " + " ".join("%*s" % (wd, n) for n, k, d, wd in self.cols) + "\n")
        self.f.flush()

    def row(self, d):
        self.f.write("  " + " ".join(("%*d" % (wd, int(d.get(n, 0)))) if k == "i" else ("%*.6e" % (wd, d.get(n, NAN)))
                                     for n, k, _, wd in self.cols) + "\n")
        self.f.flush()

    def close(self):
        self.f.close()


MPI_COL = ("mpi", "i", "program: 0 = pegSerial, 1 = pegMPI")
PAR_COLS = [("N", "i", "truncation index N"), ("P", "i", "number of points of the grating profile P"),
            ("T", "i", "threads per process (--threads / OMP_NUM_THREADS)"), ("R", "i", "number of MPI ranks (1 for pegSerial)"),
            ("eps", "f", "integration tolerance (--integrationTolerance)")]
S_COL = ("S", "i", "sweep points (energies) computed by one job = R*rounds")
RAW_COLS = [MPI_COL] + PAR_COLS + [S_COL,
    ("rep", "i", "repetition (loop) number"), ("ok", "i", "1 = run succeeded, 0 = failed / timed out"),
    ("wall_s", "f", "wall time of the whole process incl. startup [s]"),
    ("run_s", "f", "'Run time' printed by pegMPI rank 0 [s] (nan for pegSerial)"),
    ("calls", "i", "number of 'Timing Profile' blocks found in the output (should equal S)")] + \
    [(k + "_s", "f", "sum over all solver calls (and ranks) of '%s' [s]" % lab) for k, lab in CATS]
RES_COLS = [MPI_COL] + PAR_COLS + [S_COL,
    ("loops", "i", "number of runs done for this point"), ("ok_runs", "i", "number of successful runs (the average uses only these)"),
    ("time_s", "f", "wall time of the job, average over the successful runs [s]"),
    ("std_s", "f", "standard deviation of the wall time (nan for < 2 runs) [s]"),
    ("rel_std", "f", "std_s / time_s"),
    ("run_s", "f", "average 'Run time' printed by pegMPI rank 0 [s] (nan for pegSerial)")] + \
    [(k + "_s", "f", "average of the sum over all solver calls (and ranks) of '%s' [s]" % lab) for k, lab in CATS]
FAIL_COLS = [MPI_COL] + PAR_COLS + [S_COL,
    ("reason", "i", "1 = scatter of the wall time > %d %%, 2 = a run failed / timed out, 3 = skipped (a smaller case failed before)" % round(100 * MAX_REL_STD)),
    ("loops", "i", "number of runs done for this point"), ("ok_runs", "i", "number of successful runs"),
    ("time_s", "f", "average wall time of the successful runs [s]"), ("std_s", "f", "standard deviation of the wall time [s]"),
    ("rel_std", "f", "std_s / time_s"), ("min_s", "f", "fastest successful run [s]"), ("max_s", "f", "slowest successful run [s]")]


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


# ---------- the benchmark definition ----------------------------------------------------------------
def make_values(p, lo, hi, step=4):
    """values of parameter p between lo and hi (both included if they are on the grid), ascending in computational cost."""
    lo, hi = min(lo, hi), max(lo, hi)
    if p == "N":
        return list(range(lo, hi + 1, max(1, step)))
    if p == "P":
        return [2 ** k + 1 for k in range(1, 24) if lo <= 2 ** k + 1 <= hi]
    if p in "TR":
        return [2 ** k for k in range(0, 14) if lo <= 2 ** k <= hi]
    v = [float("1e%d" % k) for k in range(-16, 0) if lo <= float("1e%d" % k) <= hi]
    return sorted(v, reverse=True)                              # loose -> tight tolerance = more work


def benchmark(vals, ref, cores, no_mpi=False):
    """List of scans. In a scan ONE parameter runs through its range, the others stay at ref.
       serial: N-, P-, T-, E-scan (R=1).   mpi: for EVERY R: N-, P-, T-, E-scan (T>1 = hybrid), plus one R-scan (T=1)."""
    scans = []

    def add(prog, vary, **fixed):
        base = dict(ref, **fixed)
        pts = [dict(base, **{vary: v}) for v in vals[vary]]
        pts = [p for p in pts if p["T"] * p["R"] <= cores]
        if len(pts) > 1:
            scans.append(dict(prog=prog, vary=vary, base=base, pts=pts))

    for vary in "NPTE":
        add("serial", vary, R=1)
    if not no_mpi:
        for R in (vals["R"] or [ref["R"]]):
            for vary in "NPTE":
                add("mpi", vary, R=R)
        add("mpi", "R", T=1)
    return scans


def jobs_of(scans):
    """all distinct points of the scans as ((prog, T, R), N, P, eps), ascending along every scan."""
    seen, jobs = set(), []
    for s in scans:
        for p in s["pts"]:
            k = (s["prog"], p["N"], p["P"], p["T"], p["R"], p["E"])
            if k not in seen:
                seen.add(k)
                jobs.append(((s["prog"], p["T"], p["R"]), p["N"], p["P"], p["E"]))
    return jobs


# ---------- interactive input -----------------------------------------------------------------------
def ask(text, default, cast):
    """Enter = default, '-' = None (do not scan)."""
    while True:
        try:
            s = input("      %-6s [%g]: " % (text, default)).strip()
        except EOFError:
            return default
        if s == "":
            return default
        if s == "-":
            return None
        try:
            return cast(s)
        except ValueError:
            print("      please enter a number")


def yes(text, default=True):
    try:
        s = input("  %s [%s]: " % (text, "Y/n" if default else "y/N")).strip().lower()
    except EOFError:
        return default
    return default if not s else s[0] in "yj"


def interactive(have_mpi):
    print("\nPEG timing benchmark     (Enter = default, '-' = do not scan this parameter)\n")
    print("I do a benchmark test with:")
    vals = {}
    for p in PARAMS:
        label, cast, lo, hi = SPEC[p]
        print("  %s" % label)
        a = ask("from", lo, cast)
        b = ask("to", hi, cast) if a is not None else None
        step = (ask("step", 4, int) if p == "N" and b is not None else None) or 4
        vals[p] = make_values(p, a, b, step) if a is not None and b is not None else []
    print("\nFixed values while another parameter is scanned:")
    ref = {}
    for p in PARAMS:
        x = ask(p, REF[p], SPEC[p][1])
        ref[p] = REF[p] if x is None else x
    print()
    loops = ask("loops (repetitions) per point", LOOPS, int) or LOOPS
    no_mpi = not (have_mpi and yes("also run pegMPI (all rank counts)?"))
    return vals, ref, max(1, loops), no_mpi


# ---------- running ---------------------------------------------------------------------------------
def run_job(a, geo, outdir, cfg, N, P, E, rep):
    prog, T, R = cfg
    S = R * a.rounds
    d = os.path.join(outdir, "peg_results", "%s_N%d_P%d_T%d_R%d_E%g" % (prog, N, P, T, R, E))
    os.makedirs(d, exist_ok=True)
    out, prgf = os.path.join(d, "rep%d_output.txt" % rep), os.path.join(d, "rep%d_progress.txt" % rep)
    geo_str = geo(P)
    cmd = ["--mode", "constantIncidence", "--incidenceAngle", "88", "--eV", "--min", "100", "--max", str(100 + S - 1),
           "--increment", "1", "--gratingType", "custom", "--gratingMaterial", "Au", "--gratingPeriod", "1",
           "--gratingGeometry", geo_str, "--N", str(N), "--threads", str(T), "--integrationTolerance", "%g" % E,
           "--TE", "--measureTiming", "--outputFile", out, "--progressFile", prgf]
    cmd = [a.serial] + cmd if prog == "serial" else a.mpi_launcher.format(ranks=R).split() + [a.mpi] + cmd

    job = os.path.join(d, "rep%d_job.sh" % rep)           # job file: can be restarted by hand
    with open(job, "w") as f:
        f.write("#!/usr/bin/env bash\nexport OMP_NUM_THREADS=%d\nexport PEG_MATERIAL_DB=%s\nGEO='%s'\n"
                % (T, shlex.quote(a.material_db), geo_str))
        f.write(" ".join('"$GEO"' if c == geo_str else shlex.quote(c) for c in cmd) + "\n")

    env = dict(os.environ, OMP_NUM_THREADS=str(T), PEG_MATERIAL_DB=a.material_db)
    t0 = time.perf_counter()
    p = subprocess.Popen(["bash", job], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, start_new_session=True)
    try:
        text, _ = p.communicate(timeout=a.max_job_seconds)
        rc = p.returncode
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)                      # also kills the MPI ranks
        p.communicate()
        text, rc = "", -9
    wall = time.perf_counter() - t0
    open(os.path.join(d, "rep%d_stdout.txt" % rep), "w").write(text)
    m = re.search(r"Run time \(s\):\s*([0-9.eE+-]+)", text)
    status = ""
    try:
        status = [l.strip() for l in open(prgf) if l.startswith("status=")][0]
    except (OSError, IndexError):
        pass
    ok = rc == 0 and status == "status=succeeded"
    sums, calls = parse_timing(text)
    if calls != S:                  # e.g. launcher does not forward the output of all ranks, or a different format
        sums = {k: NAN for k in sums}
    return dict(prog=prog, N=N, P=P, T=T, R=R, E=E, S=S, rep=rep, ok=int(ok), time=wall,
                run=float(m.group(1)) if m else NAN, calls=calls, **sums)


def raw_row(r):
    return dict(mpi=int(r["prog"] == "mpi"), N=r["N"], P=r["P"], T=r["T"], R=r["R"], eps=r["E"], S=r["S"], rep=r["rep"], ok=r["ok"],
                wall_s=r["time"], run_s=r["run"], calls=r["calls"], **{k + "_s": r[k] for k, _ in CATS})


def evaluate(runs):
    """simple average / standard deviation over the successful runs of one point."""
    ok = [r for r in runs if r["ok"]]
    n = len(ok)
    mean = lambda key: sum(r[key] for r in ok) / n if n else NAN
    t = [r["time"] for r in ok]
    avg = mean("time")
    std = statistics.stdev(t) if n > 1 else NAN
    res = dict(ok_runs=n, time_s=avg, std_s=std, rel_std=std / avg if n > 1 and avg > 0 else NAN,
               min_s=min(t) if t else NAN, max_s=max(t) if t else NAN, run_s=mean("run"))
    res.update({k + "_s": mean(k) for k, _ in CATS})
    return res


# ---------- main ------------------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here) if os.path.basename(here) in ("run", "tools") else os.getcwd()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--serial", default=os.path.join(root, "pegSerial")); ap.add_argument("--mpi", default=os.path.join(root, "pegMPI"))
    ap.add_argument("--mpi-launcher", default="mpiexec -n {ranks}")
    ap.add_argument("--material-db", default=os.environ.get("PEG_MATERIAL_DB", os.path.join(root, "materialDatabase")))
    ap.add_argument("--out", help="output folder (default: <script folder>/benchmark_results/timing_<host>_<date>)")
    ap.add_argument("--cores", type=int, default=len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count(),
                    help="max T*R of a job (raise it for multi-node runs)")
    ap.add_argument("--rounds", type=int, default=ROUNDS, help="sweep points per rank and job (more = less startup overhead)")
    ap.add_argument("--max-job-seconds", type=float, default=600)
    a = ap.parse_args()
    if not os.path.exists(a.serial):
        print("pegSerial not found: %s  (use --serial)" % a.serial)
        return

    vals, ref, loops, no_mpi = interactive(os.path.exists(a.mpi))
    scans = benchmark(vals, ref, a.cores, no_mpi)
    jobs = jobs_of(scans)
    if not jobs:
        print("nothing to do: no parameter has more than one value")
        return
    n_ser, n_mpi = sum(1 for j in jobs if j[0][0] == "serial"), sum(1 for j in jobs if j[0][0] == "mpi")
    fmt = lambda p: "[" + ", ".join("%g" % v for v in vals[p]) + "]"
    outdir = a.out or os.path.join(here, "benchmark_results", "timing_%s_%s" % (socket.gethostname(), datetime.datetime.now().strftime("%Y%m%d-%H%M")))
    a.meta = ["created   : %s on host %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), socket.gethostname()),
              "programs  : pegSerial = %s ; pegMPI = %s ; launcher '%s'" % (a.serial, "(not used)" if no_mpi else a.mpi, a.mpi_launcher),
              "test case : blazed grating (alpha=2.5 deg, beta=30 deg, rounded corners) as custom profile with P points, period 1, material Au,",
              "            constantIncidence, incidence angle 88, TE, energies 100,101,... eV (S points per job), --measureTiming",
              "ranges    : N in %s ; P in %s ; T in %s ; R in %s ; eps in %s" % tuple(fmt(p) for p in PARAMS),
              "fixed     : %s (values of the other parameters while one is scanned)" % ", ".join("%s=%g" % (p, ref[p]) for p in PARAMS),
              "settings  : %d scans, %d distinct points (%d pegSerial, %d pegMPI), %d loop(s) per point, %d sweep point(s) per rank, max T*R = %d, time limit %g s per job" %
              (len(scans), len(jobs), n_ser, n_mpi, loops, a.rounds, a.cores, a.max_job_seconds)]
    print("\n" + "\n".join(a.meta[4:]))
    print("  => up to %d runs, results in %s" % (len(jobs) * loops, outdir))
    if not yes("start?"):
        return
    os.makedirs(outdir, exist_ok=True)

    depth, fine = rounded_blaze(1.0, 2.5, 30.0)
    cache = {}

    def geo(P):
        if P not in cache:
            cache[P] = geometry(depth, fine, P)
        return cache[P]

    mk = lambda name, title, about, cols: Table(os.path.join(outdir, name), title, about, cols, a.meta)
    raw = mk("raw.xyz", "PEG timing benchmark - raw data: every single run",
             ["One row per run (also failed ones, ok=0), in the order they were run.",
              "Section columns are SUMS over all solver calls (and ranks) of the job, read from the --measureTiming output."], RAW_COLS)
    res = mk("results.xyz", "PEG timing benchmark - results: average over the loops of every point",
             ["One row per point (prog, N, P, T, R, eps). Simple average over the successful loops, nothing else is calculated.",
              "Divide by S yourself if you want times per sweep point. A scan = all rows where all parameters but one equal the 'fixed' values."], RES_COLS)
    fails = mk("fails.xyz", "PEG timing benchmark - unreliable points",
               ["Points whose wall time scatters by more than %d %% (standard deviation / average), points with a failed or timed-out run," % round(100 * MAX_REL_STD),
                "and points that were skipped because a smaller case already failed. Values are over the successful runs. Such points are also in results.xyz if >= 1 run succeeded."],
               FAIL_COLS)

    failed, n_fail, bad_calls = {}, 0, False
    for i, (cfg, N, P, E) in enumerate(jobs, 1):
        prog, T, R = cfg
        point = dict(mpi=int(prog == "mpi"), N=N, P=P, T=T, R=R, eps=E, S=R * a.rounds)
        tag = "[%d/%d] %s N=%d P=%d T=%d R=%d eps=%g" % (i, len(jobs), prog, N, P, T, R, E)
        if any(n <= N and p <= P and e >= E for n, p, e in failed.get(cfg, [])):
            print(tag, "skipped (a smaller case already failed / timed out)")
            fails.row(dict(point, reason=3))
            n_fail += 1
            continue
        runs = []
        for rep in range(1, loops + 1):
            r = run_job(a, geo, outdir, cfg, N, P, E, rep)
            runs.append(r)
            raw.row(raw_row(r))
            bad_calls |= bool(r["ok"] and r["calls"] != r["S"])
            print("%s rep %d: %s %.2f s" % (tag, rep, "ok" if r["ok"] else "FAILED", r["time"]))
            if not r["ok"]:
                failed.setdefault(cfg, []).append((N, P, E))
                break
        stat = evaluate(runs)
        row = dict(point, loops=len(runs), **stat)
        if stat["ok_runs"]:
            res.row(row)
        n_bad = sum(1 for r in runs if not r["ok"])
        reason = 2 if n_bad else (1 if stat["rel_std"] > MAX_REL_STD else 0)
        if reason:
            fails.row(dict(row, reason=reason))
            n_fail += 1
            print("   -> written to fails.xyz (reason %d)" % reason)
    for t in (raw, res, fails):
        t.close()
    if bad_calls:
        print("\nWARNING: the number of 'Timing Profile' blocks != sweep points in some runs -> section times are NaN. "
              "Check peg_results/*/rep*_stdout.txt (output format? launcher forwards all ranks?)")
    print("\n%d point(s) in fails.xyz" % n_fail if n_fail else "\nall points stable")
    print("results in", outdir)


if __name__ == "__main__":
    main()