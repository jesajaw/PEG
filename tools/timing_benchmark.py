"""
AI-Generated

PEG timing benchmark: how do pegSerial / pegMPI scale with N, P, threads and ranks?

    N  truncation index     P  points of the grating profile (blazed, rounded corners, "custom" profile)
    T  threads per process  R  MPI ranks (= "nodes"; for real multi-node use --mpi-launcher and --cores)

HOW THE PARAMETERS ARE PASSED
  1. The ranges of the four parameters are defined once in the block GLOBAL SETTINGS below (python ranges + reference point).
  2. benchmark() turns them into SCANS: in a scan ONE parameter is scaled up through its whole range, the other three
     stay at their reference value (REF).  Run for pegSerial (R=1) and for pegMPI at EVERY rank count R of the range:
         serial:      N-scan, P-scan, T-scan
         mpi, per R:  N-scan, P-scan, T-scan (T>1 = hybrid)      + one R-scan (T=1)
     => pegMPI is always run at least as often as pegSerial (R=1 repeats the serial scans) and more often for R>1.
  3. Every point (prog, N, P, T, R) of a scan is one job. A job calls pegSerial / pegMPI with
         --N <N>   --threads <T> (+ OMP_NUM_THREADS)   --gratingGeometry <P profile points>
         pegMPI:  <launcher with R ranks> pegMPI ...     S = R*rounds sweep points (energies), so every rank has the same work
     and measures wall time + the --measureTiming sections from the terminal output.

OUTPUT  (default  <folder of this script>/benchmark_results/timing_<host>_<date>/ ; no plots, only plottable data)
  All files are plain-text .xyz tables: '#' header (title, context, description of every column), then whitespace
  separated numeric columns (nan = not available): gnuplot, pgfplots, numpy.loadtxt, pandas.read_csv(sep=r"\\s+", comment="#").
    timing_raw.xyz         every single run (in the order they were run)
    summary.xyz            one block per scan (blocks separated by 2 blank lines = gnuplot 'index'), fastest repetition:
                           time, time/point, speed-up, efficiency, time/point of every --measureTiming section
    fits.xyz               per scan: exponent e in  t ~ x^e  (x = scaled parameter), total + every section
    local_exponents.xyz    exponent between neighbouring points of a scan (shows e.g. an N^3 -> N^4 crossover)
    peg_results/<job>/     PEG output, progress file and raw terminal output of every run

--measureTiming is only printed to the terminal, so it is read from the program output: every solver call prints a
"Timing Profile" block; the blocks of all calls (and ranks) are summed per section (alloc, refr, layers, integ, matrix,
rayleigh, effs, solver). time = wall time of the whole process (same definition for pegSerial and pegMPI);
"run" = "Run time" printed by pegMPI rank 0. "per point" = divided by S.
speed-up = throughput / throughput of (pegSerial, T=1) at the same N,P.
Needs only python3.   Run:   python3 tools/timing_benchmark.py [--quick]
"""
import argparse, datetime, math, os, re, signal, socket, statistics, subprocess, time

# ============================== GLOBAL SETTINGS ==============================
N_RANGE = range(8, 33, 4)      # N = 8, 12, 16, ..., 32
P_EXPS  = range(3, 10)         # P = 2**k + 1  ->  9, 17, 33, ..., 513 profile points
T_EXPS  = range(0, 5)          # T = 2**k      ->  1, 2, 4, 8, 16 threads   (only values <= cores are used)
R_EXPS  = range(0, 5)          # R = 2**k      ->  1, 2, 4, 8, 16 MPI ranks (only values <= cores are used)
REF = dict(N=16, P=65, T=1, R=1)   # while one parameter is scaled, the others stay here (the point is added to every scan)
ROUNDS, REPEATS = 2, 3         # sweep points per rank and job / repetitions per point (fastest is used)
QUICK = dict(N=range(8, 25, 8), P=range(3, 9, 2), T=range(0, 3), R=range(0, 3), repeats=1)    # used by --quick
# =============================================================================

PARAMS = ["N", "P", "T", "R"]
CATS = [("alloc", "Allocate Memory"), ("refr", "Look up refractive index"),
        ("layers", "Compute alpha, beta values and layers"), ("integ", "Numerically integrating trial solutions"),
        ("matrix", "Matrix operations"), ("rayleigh", "Computing Rayleigh coeffients B_n"),
        ("effs", "Compute and package efficiencies"), ("solver", "Total (solver) time")]
SECS = [k + "_pp" for k, _ in CATS]
SEC_ID = ["tpp"] + SECS                                  # section index used in fits.xyz / local_exponents.xyz
SEC_TXT = "0 = total time per point (wall)" + "".join(", %d = %s" % (i, lab) for i, (_, lab) in enumerate(CATS, 1))


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
    """'#' header (title, context, one line per column), then whitespace separated numeric columns."""
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
        self.f.write("  " + " ".join(("%*d" % (wd, int(d[n]))) if k == "i" else ("%*.6e" % (wd, d[n])) for n, k, _, wd in self.cols) + "\n")
        self.f.flush()

    def blank(self, n=1):
        self.f.write("\n" * n)

    def close(self):
        self.f.close()


MPI_COL = ("mpi", "i", "program: 0 = pegSerial, 1 = pegMPI")
PAR_COLS = [("N", "i", "truncation index N"), ("P", "i", "number of points of the grating profile P"),
            ("T", "i", "threads per process (--threads / OMP_NUM_THREADS)"), ("R", "i", "number of MPI ranks (1 for pegSerial)")]
RAW_COLS = [MPI_COL] + PAR_COLS + [
    ("S", "i", "sweep points (energies) computed by this job = R*rounds"),
    ("rep", "i", "repetition number of this point"), ("ok", "i", "1 = run succeeded, 0 = failed / timed out"),
    ("wall_s", "f", "wall time of the whole process incl. startup [s]"),
    ("run_s", "f", "'Run time' printed by pegMPI rank 0 [s] (nan for pegSerial)"),
    ("calls", "i", "number of 'Timing Profile' blocks found in the output (should equal S)")] + \
    [(k + "_s", "f", "sum over all solver calls (and ranks) of '%s' [s]" % lab) for k, lab in CATS]
SUM_COLS = [("scan", "i", "scaled parameter of this block: 1 = N, 2 = P, 3 = T, 4 = R (x axis = this column)"), MPI_COL] + PAR_COLS + [
    ("S", "i", "sweep points (energies) computed by this job = R*rounds"),
    ("time_s", "f", "wall time of the job, fastest repetition [s]"), ("time_med_s", "f", "wall time of the job, median over repetitions [s]"),
    ("run_s", "f", "'Run time' printed by pegMPI rank 0 of the fastest repetition [s] (nan for pegSerial)"),
    ("tpp_s", "f", "wall time per sweep point = time_s / S [s]"),
    ("speedup", "f", "throughput (S/time) relative to pegSerial, T=1, same N and P"),
    ("efficiency", "f", "speedup / (T*R)")] + [(k + "_pp", "f", "'%s' per sweep point [s]" % lab) for k, lab in CATS]
SCAN_COLS = [("scan", "i", "scaled parameter: 1 = N, 2 = P, 3 = T, 4 = R"), MPI_COL] + \
            [(n, k, d + " (0 = this is the scaled parameter)") for n, k, d in PAR_COLS] + \
            [("section", "i", "quantity (time per point): " + SEC_TXT)]
FIT_COLS = SCAN_COLS + [("exponent", "f", "exponent e in  t = c * x^e,  x = scaled parameter (for T and R the ideal is e = -1)"),
                        ("r2", "f", "R^2 of the least-squares fit in log-log"), ("npts", "i", "number of points of the scan used")]
LOC_COLS = SCAN_COLS + [("x_lo", "i", "lower value of the scaled parameter of the pair"), ("x_hi", "i", "upper value of the scaled parameter of the pair"),
                        ("exponent", "f", "local exponent = dlog(t) / dlog(x) between x_lo and x_hi")]


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
def parameter_values(cores, quick=False):
    """the values of the four parameters: the global ranges (+ reference value), T and R limited by the cores."""
    q = QUICK if quick else dict(N=N_RANGE, P=P_EXPS, T=T_EXPS, R=R_EXPS)
    vals = {"N": set(q["N"]), "P": {2 ** k + 1 for k in q["P"]}, "T": {2 ** k for k in q["T"]}, "R": {2 ** k for k in q["R"]}}
    return {p: sorted(v | {REF[p]}) if p in ("N", "P") else sorted(x for x in v | {REF[p]} if x <= cores) for p, v in vals.items()}


def benchmark(vals, cores, no_mpi=False):
    """The benchmark = list of scans. In a scan ONE parameter runs through its whole range, the others stay at REF.
       serial: N-, P-, T-scan (R=1).   mpi: for EVERY R: N-, P-, T-scan (T>1 = hybrid), plus one R-scan (T=1)."""
    scans = []

    def add(prog, vary, **fixed):
        base = dict(REF, **fixed)
        pts = [dict(base, **{vary: v}) for v in vals[vary]]
        pts = [p for p in pts if p["T"] * p["R"] <= cores]
        if len(pts) > 1:
            scans.append(dict(prog=prog, vary=vary, base=base, pts=pts))

    for vary in "NPT":
        add("serial", vary, R=1)
    if not no_mpi:
        for R in vals["R"]:
            for vary in "NPT":
                add("mpi", vary, R=R)
        add("mpi", "R", T=1)
    return scans


def jobs_of(scans):
    """all distinct (cfg, N, P) of the scans, cfg = (prog, T, R); ascending along every scan."""
    seen, jobs = set(), []
    for s in scans:
        for p in s["pts"]:
            k = (s["prog"], p["N"], p["P"], p["T"], p["R"])
            if k not in seen:
                seen.add(k)
                jobs.append(((s["prog"], p["T"], p["R"]), p["N"], p["P"]))
    return jobs


# ---------- running ---------------------------------------------------------------------------------
def run_job(a, geo, outdir, cfg, N, P, rep):
    prog, T, R = cfg
    S = R * a.rounds
    d = os.path.join(outdir, "peg_results", "%s_N%d_P%d_T%d_R%d" % (prog, N, P, T, R))
    os.makedirs(d, exist_ok=True)
    out, prgf = os.path.join(d, "rep%d_output.txt" % rep), os.path.join(d, "rep%d_progress.txt" % rep)
    cmd = ["--mode", "constantIncidence", "--incidenceAngle", "88", "--eV", "--min", "100", "--max", str(100 + S - 1),
           "--increment", "1", "--gratingType", "custom", "--gratingMaterial", "Au", "--gratingPeriod", "1",
           "--gratingGeometry", geo(P), "--N", str(N), "--threads", str(T), "--TE", "--measureTiming",
           "--outputFile", out, "--progressFile", prgf]
    cmd = [a.serial] + cmd if prog == "serial" else a.mpi_launcher.format(ranks=R).split() + [a.mpi] + cmd
    env = dict(os.environ, OMP_NUM_THREADS=str(T), PEG_MATERIAL_DB=a.material_db)
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, start_new_session=True)
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
        sums = {k: float("nan") for k in sums}
    return dict(prog=prog, N=N, P=P, T=T, R=R, S=S, rep=rep, ok=int(ok), time=wall,
                run=float(m.group(1)) if m else float("nan"), calls=calls, **sums)


def raw_row(r):
    return dict(mpi=int(r["prog"] == "mpi"), N=r["N"], P=r["P"], T=r["T"], R=r["R"], S=r["S"], rep=r["rep"], ok=r["ok"],
                wall_s=r["time"], run_s=r["run"], calls=r["calls"], **{k + "_s": r[k] for k, _ in CATS})


# ---------- analysis (pure python) ------------------------------------------------------------------
def lstsq(A, z):
    """least squares via normal equations + Gauss-Jordan (few unknowns only)."""
    k = len(A[0])
    M = [[sum(r[i] * r[j] for r in A) for j in range(k)] + [sum(r[i] * zv for r, zv in zip(A, z))] for i in range(k)]
    for i in range(k):
        p = max(range(i, k), key=lambda r: abs(M[r][i]))
        M[i], M[p] = M[p], M[i]
        if abs(M[i][i]) < 1e-12:
            return None
        for r in range(k):
            if r != i:
                f = M[r][i] / M[i][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [M[i][k] / M[i][i] for i in range(k)]


def fit(x, y):
    """log y = c + e log x  ->  (e, R^2, n) or None."""
    ok = [i for i in range(len(y)) if y[i] > 0 and x[i] > 0]
    if len(ok) < 3 or len({x[i] for i in ok}) < 3:
        return None
    A = [[1.0, math.log(x[i])] for i in ok]
    z = [math.log(y[i]) for i in ok]
    beta = lstsq(A, z)
    if beta is None:
        return None
    zm = sum(z) / len(z)
    ss = sum((v - zm) ** 2 for v in z)
    sr = sum((v - beta[0] - beta[1] * row[1]) ** 2 for v, row in zip(z, A))
    return beta[1], (1 - sr / ss if ss > 0 else 1.0), len(ok)


def analyse(outdir, rows, scans, a):
    groups = {}
    for r in rows:
        if r["ok"]:
            groups.setdefault((r["prog"], r["N"], r["P"], r["T"], r["R"]), []).append(r)
    if not groups:
        print("no successful run -> nothing to analyse; see peg_results/*/rep*_stdout.txt")
        return
    best = {k: min(v, key=lambda r: r["time"]) for k, v in groups.items()}
    summ = {}
    for (prog, N, P, T, R), r in best.items():
        t, S = r["time"], r["S"]
        base = best.get(("serial", N, P, 1, 1))
        sp = (S / t) / (base["S"] / base["time"]) if base else float("nan")
        summ[(prog, N, P, T, R)] = dict(
            mpi=int(prog == "mpi"), N=N, P=P, T=T, R=R, S=S, time_s=t, time_med_s=statistics.median(x["time"] for x in groups[(prog, N, P, T, R)]),
            run_s=r["run"], tpp_s=t / S, tpp=t / S, speedup=sp, efficiency=sp / (T * R), **{k + "_pp": r[k] / S for k, _ in CATS})

    t = Table(os.path.join(outdir, "summary.xyz"), "PEG timing benchmark - summary: one block per scan",
              ["In a scan ONE parameter (column 'scan': 1 = N, 2 = P, 3 = T, 4 = R) is scaled up through its range, the other parameters stay at",
               "their reference values. Every block starts with a '# scan:' line; blocks are separated by 2 blank lines (gnuplot: 'index i').",
               "x axis of a block = the column of the scaled parameter. Fastest of the repetitions per point; times per point are per sweep point",
               "(one energy of the spectrum). speedup/efficiency refer to pegSerial with T=1 at the same N and P."], SUM_COLS, a.meta)
    fits, local = [], []
    print("\n--- exponent e in  time per point ~ x^e  (x = the scaled parameter of the scan) ---")
    for si, sc in enumerate(scans):
        prog, vary, base = sc["prog"], sc["vary"], sc["base"]
        g = [summ[k] for k in ((prog, p["N"], p["P"], p["T"], p["R"]) for p in sc["pts"]) if k in summ]
        if not g:
            continue
        fixed = ", ".join("%s=%d" % (q, base[q]) for q in PARAMS if q != vary)
        if si:
            t.blank(2)
        t.f.write("# scan: %s, %s scaled up, fixed: %s\n" % ("pegMPI" if prog == "mpi" else "pegSerial", vary, fixed))
        head = dict(scan=PARAMS.index(vary) + 1, mpi=int(prog == "mpi"), **{q: (0 if q == vary else base[q]) for q in PARAMS})
        x, ex = [s[vary] for s in g], {}
        for s in g:
            t.row(dict(s, scan=head["scan"]))
        for sec_i, sec in enumerate(SEC_ID):
            y = [s[sec] for s in g]
            f = fit(x, y)
            if f:
                ex[sec] = f[0]
                fits.append(dict(head, section=sec_i, exponent=f[0], r2=f[1], npts=f[2]))
            for s0, s1 in zip(g, g[1:]):
                if s0[sec] > 0 and s1[sec] > 0:
                    local.append(dict(head, section=sec_i, x_lo=s0[vary], x_hi=s1[vary], exponent=math.log(s1[sec] / s0[sec]) / math.log(s1[vary] / s0[vary])))
        line = "%-7s %s-scan (%s):" % (prog, vary, fixed)
        if vary in "NP":
            line += "  " + "  ".join("%s %.2f" % (nm, ex[k]) for nm, k in (("total", "tpp"), ("integ", "integ_pp"), ("matrix", "matrix_pp"), ("solver", "solver_pp")) if k in ex)
        else:
            line += "  speed-up %.2f at %s=%d (efficiency %.0f %%)" % (g[-1]["speedup"], vary, g[-1][vary], 100 * g[-1]["efficiency"])
        print(line)
    t.close()
    for name, title, about, cols, data in (
            ("fits.xyz", "PEG timing benchmark - measured scaling exponents (power-law fits in log-log, one per scan)",
             ["Per scan and quantity: least-squares fit  t = c * x^e  of the TIME PER POINT over the scaled parameter x (needs >= 3 points).",
              "section 0 uses wall time incl. process startup, which flattens the exponent at small N; sections 1..8 are the --measureTiming sections.",
              "For T and R scans the exponent of an ideal parallelisation is e = -1."], FIT_COLS, fits),
            ("local_exponents.xyz", "PEG timing benchmark - local exponents between neighbouring scan points",
             ["exponent = dlog(t) / dlog(x) between two neighbouring points of a scan (t = time per point, x = scaled parameter).",
              "Shows whether the scaling changes along the scan (e.g. approaches N^3 or N^4 for large N)."], LOC_COLS, local)):
        tb = Table(os.path.join(outdir, name), title, about, cols, a.meta)
        for x in data:
            tb.row(x)
        tb.close()
    if any(math.isnan(s["solver_pp"]) for s in summ.values()):
        print("\nWARNING: the number of 'Timing Profile' blocks != sweep points in some runs -> section times are NaN. "
              "Check peg_results/*/rep*_stdout.txt (output format? launcher forwards all ranks?)")
    bad = sum(1 for r in rows if not r["ok"])
    if bad:
        print("WARNING: %d run(s) failed or timed out and are excluded" % bad)


# ---------- main ------------------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here) if os.path.basename(here) in ("run", "tools") else os.getcwd()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--serial", default=os.path.join(root, "pegSerial")); ap.add_argument("--mpi", default=os.path.join(root, "pegMPI"))
    ap.add_argument("--no-mpi", action="store_true"); ap.add_argument("--mpi-launcher", default="mpiexec -n {ranks}")
    ap.add_argument("--material-db", default=os.environ.get("PEG_MATERIAL_DB", os.path.join(root, "materialDatabase")))
    ap.add_argument("--out", help="output folder (default: <script folder>/benchmark_results/timing_<host>_<date>)")
    ap.add_argument("--quick", action="store_true", help="small ranges (QUICK in the global settings), 1 repeat")
    ap.add_argument("--cores", type=int, default=len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count(),
                    help="max T*R of a job (raise it for multi-node runs)")
    ap.add_argument("--rounds", type=int, default=ROUNDS, help="sweep points per rank and job (more = less startup overhead)")
    ap.add_argument("--repeats", type=int, default=None); ap.add_argument("--max-job-seconds", type=float, default=600)
    a = ap.parse_args()
    a.repeats = a.repeats or (QUICK["repeats"] if a.quick else REPEATS)
    if not os.path.exists(a.mpi) and not a.no_mpi:
        print("pegMPI not found -> only pegSerial"); a.no_mpi = True
    vals = parameter_values(a.cores, a.quick)
    scans = benchmark(vals, a.cores, a.no_mpi)
    jobs = jobs_of(scans)
    n_ser, n_mpi = sum(1 for j in jobs if j[0][0] == "serial"), sum(1 for j in jobs if j[0][0] == "mpi")
    outdir = a.out or os.path.join(here, "benchmark_results", "timing_%s_%s" % (socket.gethostname(), datetime.datetime.now().strftime("%Y%m%d-%H%M")))
    os.makedirs(outdir, exist_ok=True)
    a.meta = ["created   : %s on host %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), socket.gethostname()),
              "programs  : pegSerial = %s ; pegMPI = %s ; launcher '%s'" % (a.serial, "(not used)" if a.no_mpi else a.mpi, a.mpi_launcher),
              "test case : blazed grating (alpha=2.5 deg, beta=30 deg, rounded corners) as custom profile with P points, period 1, material Au,",
              "            constantIncidence, incidence angle 88, TE, energies 100,101,... eV (S points per job), --measureTiming",
              "ranges    : N in %s ; P in %s" % (vals["N"], vals["P"]),
              "            T in %s ; R in %s" % (vals["T"], vals["R"]),
              "reference : %s (the other parameters stay here while one is scaled)" % ", ".join("%s=%d" % (p, REF[p]) for p in PARAMS),
              "settings  : %d scans, %d distinct points (%d pegSerial, %d pegMPI), %d repetition(s) each, %d sweep point(s) per rank, max T*R = %d, time limit %g s per job" %
              (len(scans), len(jobs), n_ser, n_mpi, a.repeats, a.rounds, a.cores, a.max_job_seconds)]
    print("\n".join(a.meta[4:]))
    depth, fine = rounded_blaze(1.0, 2.5, 30.0)
    cache = {}
    geo = lambda P: cache.setdefault(P, geometry(depth, fine, P)) if P not in cache else cache[P]

    rows, failed = [], {}
    raw = Table(os.path.join(outdir, "timing_raw.xyz"), "PEG timing benchmark - raw data: every single run",
                ["One row per run (also failed ones, ok=0), in the order they were run (scan after scan).",
                 "Section columns are SUMS over all solver calls (and ranks) of the job, read from the --measureTiming output."], RAW_COLS, a.meta)
    for i, (cfg, N, P) in enumerate(jobs, 1):
        if any(n <= N and p <= P for n, p in failed.get(cfg, [])):
            print("[%d/%d] skip %s T=%d R=%d N=%d P=%d (a smaller case already failed / timed out)" % ((i, len(jobs)) + cfg + (N, P)))
            continue
        for rep in range(1, a.repeats + 1):
            r = run_job(a, geo, outdir, cfg, N, P, rep)
            rows.append(r); raw.row(raw_row(r))
            print("[%d/%d] %s N=%d P=%d T=%d R=%d rep %d: %s %.2f s" % (i, len(jobs), cfg[0], N, P, cfg[1], cfg[2], rep, "ok" if r["ok"] else "FAILED", r["time"]))
            if not r["ok"]:
                failed.setdefault(cfg, []).append((N, P)); break
    raw.close()
    analyse(outdir, rows, scans, a)
    print("results in", outdir)


if __name__ == "__main__":
    main()