"""PEG timing benchmark: how do pegSerial / pegMPI scale with N, P, threads and ranks?

    N  truncation index     P  points of the grating profile (blazed, rounded corners, "custom" profile)
    T  threads per process  R  MPI ranks (= "nodes"; for real multi-node use --mpi-launcher + --cores)

Output folder (default  <folder of this script>/benchmark_results/timing_<host>_<date>/ ), NO plots, only plottable data.
All result files are plain-text .xyz tables: '#' header lines (title, context, description of every column),
then whitespace separated numeric columns (nan = not available). Directly readable by gnuplot, pgfplots,
numpy.loadtxt, pandas.read_csv(sep=r"\\s+", comment="#").
    timing_raw.xyz         every single run
    summary.xyz            one row per configuration (fastest repetition): time, time/point, speed-up, efficiency,
                           time/point of every --measureTiming section; blocks separated by blank lines (see header)
    fits.xyz               power-law exponents  t ~ N^a * P^b  (total + every section), joint and per P
    local_exponents.xyz    exponent between neighbouring N values (shows e.g. an N^3 -> N^4 crossover)
    peg_results/<job>/     PEG output, progress file and raw terminal output of every run

--measureTiming is only printed to the terminal, so it is read from the program output: every solver call prints a
"Timing Profile" block; the blocks of all calls (and ranks) are summed per section (alloc, refr, layers, integ, matrix,
rayleigh, effs, solver). time = wall time of the whole process (same definition for pegSerial and pegMPI);
"run" = "Run time" printed by pegMPI rank 0. S = R*rounds sweep points per job, "per point" = divided by S.
speed-up = throughput / throughput of (pegSerial, T=1) at the same N,P.

Designs:  full  every (N,P) for every (T,R) configuration      lean  (N,P) grid for serial T=1 only,
                                                                      all other configurations only over N at P0
Needs only python3 (no numpy / matplotlib).  Run:   python3 tools/timing_benchmark.py [--quick] [--design lean]
"""
import argparse, datetime, math, os, re, signal, socket, statistics, subprocess, time

# ======================= GLOBAL SETTINGS: ranges of the four parameters =======================
# (command line options --N-list --P-list --T-list --R-list override them for a single run)
N_LIST = [6, 8, 12, 16, 24, 32]                        # truncation index N
P_LIST = [9, 17, 33, 65, 129, 257, 513, 1025]          # points of the grating profile P
T_LIST = None                                          # threads per process T   (None = 1,2,4,... up to --cores)
R_LIST = None                                          # MPI ranks R             (None = 1,2,4,... up to --cores)
REF_N, REF_P = 12, 65                                  # reference point (nearest list value is used) for the scans over T, R
QUICK_N, QUICK_P = [6, 10, 16, 24], [9, 65, 513]       # lists used by --quick
# ==============================================================================================

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


SEC_DESC_RAW = [(k + "_s", "f", "sum over all solver calls (and ranks) of '%s' [s]" % lab) for k, lab in CATS]
RAW_COLS = [("mpi", "i", "program: 0 = pegSerial, 1 = pegMPI"), ("N", "i", "truncation index N"),
            ("P", "i", "number of points of the grating profile P"), ("T", "i", "threads per process (--threads / OMP_NUM_THREADS)"),
            ("R", "i", "number of MPI ranks (1 for pegSerial)"), ("S", "i", "sweep points (energies) computed by this job = R*rounds"),
            ("rep", "i", "repetition number of this configuration"), ("ok", "i", "1 = run succeeded, 0 = failed / timed out"),
            ("wall_s", "f", "wall time of the whole process incl. startup [s]"),
            ("run_s", "f", "'Run time' printed by pegMPI rank 0 [s] (nan for pegSerial)"),
            ("calls", "i", "number of 'Timing Profile' blocks found in the output (should equal S)")] + SEC_DESC_RAW
SUM_COLS = [c for c in RAW_COLS[:6]] + [
    ("time_s", "f", "wall time of the job, fastest repetition [s]"), ("time_med_s", "f", "wall time of the job, median over repetitions [s]"),
    ("run_s", "f", "'Run time' printed by pegMPI rank 0 of the fastest repetition [s] (nan for pegSerial)"),
    ("tpp_s", "f", "wall time per sweep point = time_s / S [s]"),
    ("speedup", "f", "throughput (S/time) relative to pegSerial, T=1, same N and P"),
    ("efficiency", "f", "speedup / (T*R)")] + [(k + "_pp", "f", "'%s' per sweep point [s]" % lab) for k, lab in CATS]
FIT_COLS = [("mpi", "i", "program: 0 = pegSerial, 1 = pegMPI"), ("T", "i", "threads per process"), ("R", "i", "MPI ranks"),
            ("section", "i", "fitted quantity (time per point): " + SEC_TXT),
            ("P", "i", "profile points P of the fit; 0 = joint fit over all P (then exp_P is given)"),
            ("exp_N", "f", "exponent a in  t ~ N^a"), ("exp_P", "f", "exponent b in  t ~ N^a * P^b (nan for fits at fixed P)"),
            ("r2", "f", "R^2 of the least-squares fit in log-log"), ("npts", "i", "number of data points used")]
LOC_COLS = [FIT_COLS[0], FIT_COLS[1], FIT_COLS[2], FIT_COLS[3], ("P", "i", "profile points P"),
            ("N_lo", "i", "lower N of the pair"), ("N_hi", "i", "upper N of the pair"),
            ("exp_N", "f", "local exponent = dlog(t) / dlog(N) between N_lo and N_hi")]


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
    """list of (cfg, N, P) with cfg = (prog, T, R); per cfg ascending in N then P (needed for the time-limit skip)."""
    cfgs = [("serial", T, 1) for T in a.T_list]
    if not a.no_mpi:
        cfgs += [("mpi", T, R) for R in a.R_list for T in a.T_list]
    jobs = []
    for cfg in cfgs:
        if cfg[1] * cfg[2] > a.cores:
            continue
        if a.design == "full" or cfg == ("serial", 1, 1):
            pairs = [(N, P) for N in sorted(a.N_list) for P in sorted(a.P_list)]
        else:
            pairs = [(N, a.P0) for N in sorted(a.N_list)]
        jobs += [(cfg, N, P) for N, P in pairs]
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


def fit(cols, y):
    """log y = c + sum_i e_i log x_i  ->  (coefficients [c, e_1, ...], R^2, n) or None."""
    ok = [i for i in range(len(y)) if y[i] > 0 and all(c[i] > 0 for c in cols)]
    if len(ok) < len(cols) + 2 or any(len({c[i] for i in ok}) < 2 for c in cols):
        return None
    A = [[1.0] + [math.log(c[i]) for c in cols] for i in ok]
    z = [math.log(y[i]) for i in ok]
    beta = lstsq(A, z)
    if beta is None:
        return None
    zm = sum(z) / len(z)
    ss = sum((v - zm) ** 2 for v in z)
    sr = sum((v - sum(b * x for b, x in zip(beta, row))) ** 2 for v, row in zip(z, A))
    return beta, (1 - sr / ss if ss > 0 else 1.0), len(ok)


def analyse(outdir, rows, a):
    groups = {}
    for r in rows:
        if r["ok"]:
            groups.setdefault((r["prog"], r["N"], r["P"], r["T"], r["R"]), []).append(r)
    if not groups:
        print("no successful run -> nothing to analyse; see peg_results/*/rep*_stdout.txt")
        return
    best = {k: min(v, key=lambda r: r["time"]) for k, v in groups.items()}
    summ = []
    for (prog, N, P, T, R), r in sorted(best.items()):
        t, S = r["time"], r["S"]
        base = best.get(("serial", N, P, 1, 1))
        sp = (S / t) / (base["S"] / base["time"]) if base else float("nan")
        summ.append(dict(prog=prog, N=N, P=P, T=T, R=R, S=S, time=t, time_med=statistics.median(x["time"] for x in groups[(prog, N, P, T, R)]),
                         run=r["run"], tpp=t / S, speedup=sp, efficiency=sp / (T * R), **{k + "_pp": r[k] / S for k, _ in CATS}))
    cfgs = sorted({(s["prog"], s["T"], s["R"]) for s in summ})
    sel = lambda **kw: sorted([s for s in summ if all(s[k] == v for k, v in kw.items())], key=lambda s: (s["N"], s["P"], s["T"], s["R"]))
    mp = lambda prog: int(prog == "mpi")

    t = Table(os.path.join(outdir, "summary.xyz"), "PEG timing benchmark - summary: one row per configuration (N, P, T, R)",
              ["Fastest of the repetitions per configuration. Times per point are per sweep point (one energy of the spectrum).",
               "Layout for plotting: configurations (mpi,T,R) are separated by 2 blank lines (gnuplot: 'index i');",
               "inside a configuration, rows of equal N form a block separated by 1 blank line (gnuplot: splot / pm3d surface over N,P).",
               "speedup/efficiency refer to pegSerial with T=1 at the same N and P."], SUM_COLS, a.meta)
    for i, (prog, T, R) in enumerate(cfgs):
        if i:
            t.blank(2)
        t.f.write("# configuration: %s, T=%d threads, R=%d MPI ranks\n" % ("pegMPI" if prog == "mpi" else "pegSerial", T, R))
        g, lastN = sel(prog=prog, T=T, R=R), None
        for s in g:
            if lastN is not None and s["N"] != lastN:
                t.blank()
            lastN = s["N"]
            t.row(dict(s, mpi=mp(prog), time_s=s["time"], time_med_s=s["time_med"], run_s=s["run"], tpp_s=s["tpp"]))
    t.close()

    fits, local = [], []
    for prog, T, R in cfgs:
        g = sel(prog=prog, T=T, R=R)
        for si, sec in enumerate(SEC_ID):
            f = fit([[s["N"] for s in g], [s["P"] for s in g]], [s[sec] for s in g])
            if f:
                fits.append(dict(mpi=mp(prog), T=T, R=R, section=si, P=0, exp_N=f[0][1], exp_P=f[0][2], r2=f[1], npts=f[2]))
            for P in sorted({s["P"] for s in g}):
                gp = [s for s in g if s["P"] == P]
                f = fit([[s["N"] for s in gp]], [s[sec] for s in gp])
                if f:
                    fits.append(dict(mpi=mp(prog), T=T, R=R, section=si, P=P, exp_N=f[0][1], exp_P=float("nan"), r2=f[1], npts=f[2]))
                for s0, s1 in zip(gp, gp[1:]):
                    if s0[sec] > 0 and s1[sec] > 0:
                        local.append(dict(mpi=mp(prog), T=T, R=R, section=si, P=P, N_lo=s0["N"], N_hi=s1["N"],
                                          exp_N=math.log(s1[sec] / s0[sec]) / math.log(s1["N"] / s0["N"])))
    for name, title, about, cols, data in (
            ("fits.xyz", "PEG timing benchmark - measured scaling exponents (power-law fits in log-log)",
             ["Least-squares fits of the TIME PER POINT. Joint fit (P=0):  t = c * N^exp_N * P^exp_P  over all (N,P) of one configuration;",
              "P>0: fit  t = c * N^exp_N  at fixed P.  Rows with section 1..8 give the exponents of the individual --measureTiming sections.",
              "Note: section 0 uses wall time incl. process startup, which flattens the exponent at small N; sections 1..8 do not."],
             FIT_COLS, fits),
            ("local_exponents.xyz", "PEG timing benchmark - local exponents between neighbouring N",
             ["exp_N = dlog(t) / dlog(N) between N_lo and N_hi at fixed P (t = time per point of the given section).",
              "Shows whether the scaling changes with N (e.g. approaches N^3 or N^4 for large N); section: " + SEC_TXT],
             LOC_COLS, local)):
        tb = Table(os.path.join(outdir, name), title, about, cols, a.meta)
        for x in data:
            tb.row(x)
        tb.close()

    print("\n--- time per point  t ~ N^a * P^b  (joint fit over the whole (N,P) grid) ---")
    for x in fits:
        if x["P"] == 0 and (x["section"] == 0 or (x["mpi"], x["T"], x["R"]) == (0, 1, 1)):
            print("%s T=%-2d R=%-2d %-11s a=%.2f  b=%.2f  R2=%.3f" % ("mpi   " if x["mpi"] else "serial", x["T"], x["R"], SEC_ID[x["section"]], x["exp_N"], x["exp_P"], x["r2"]))
    Nm = max(s["N"] for s in summ)
    print("\n--- parallel speed-up at N=%d, P=%d ---" % (Nm, a.P0))
    for s in sel(N=Nm, P=a.P0):
        print("%-7s T=%-2d R=%-2d t/pt=%.3g s  speed-up %.2f  efficiency %.0f %%" % (s["prog"], s["T"], s["R"], s["tpp"], s["speedup"], 100 * s["efficiency"]))
    if any(math.isnan(s["solver_pp"]) for s in summ):
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
    ap.add_argument("--quick", action="store_true", help="small lists, 1 repeat")
    ap.add_argument("--design", choices=["full", "lean"], default="full")
    ap.add_argument("--cores", type=int, default=len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count(),
                    help="max T*R of a job (raise it for multi-node runs)")
    ap.add_argument("--N-list"); ap.add_argument("--P-list"); ap.add_argument("--T-list"); ap.add_argument("--R-list")
    ap.add_argument("--rounds", type=int, default=2, help="sweep points per rank and job (more = less startup overhead)")
    ap.add_argument("--repeats", type=int, default=None); ap.add_argument("--max-job-seconds", type=float, default=600)
    a = ap.parse_args()
    ints = lambda s, d: [int(x) for x in s.split(",")] if s else d
    a.N_list = ints(a.N_list, QUICK_N if a.quick else N_LIST)
    a.P_list = ints(a.P_list, QUICK_P if a.quick else P_LIST)
    a.T_list = sorted(set(ints(a.T_list, T_LIST or pow2(a.cores))) | {1})      # T=1, R=1 always included (baseline)
    a.R_list = sorted(set(ints(a.R_list, R_LIST or pow2(a.cores))) | {1})
    a.N0 = min(a.N_list, key=lambda n: abs(n - REF_N)); a.P0 = min(a.P_list, key=lambda p: abs(p - REF_P))
    a.repeats = a.repeats or (1 if a.quick else 3)
    if not os.path.exists(a.mpi) and not a.no_mpi:
        print("pegMPI not found -> only pegSerial"); a.no_mpi = True
    outdir = a.out or os.path.join(here, "benchmark_results", "timing_%s_%s" % (socket.gethostname(), datetime.datetime.now().strftime("%Y%m%d-%H%M")))
    os.makedirs(outdir, exist_ok=True)
    a.meta = ["created   : %s on host %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), socket.gethostname()),
              "programs  : pegSerial = %s ; pegMPI = %s ; launcher '%s'" % (a.serial, "(not used)" if a.no_mpi else a.mpi, a.mpi_launcher),
              "test case : blazed grating (alpha=2.5 deg, beta=30 deg, rounded corners) as custom profile with P points, period 1, material Au,",
              "            constantIncidence, incidence angle 88, TE, energies 100,101,... eV (S points per job), --measureTiming",
              "design    : %s ;  N in %s ;  P in %s" % (a.design, a.N_list, a.P_list),
              "            T in %s ;  R in %s ;  reference point N0=%d, P0=%d" % (a.T_list, a.R_list, a.N0, a.P0),
              "settings  : %d repetition(s) per configuration, %d sweep point(s) per rank, max T*R = %d, time limit %g s per job" %
              (a.repeats, a.rounds, a.cores, a.max_job_seconds)]
    depth, fine = rounded_blaze(1.0, 2.5, 30.0)
    cache = {}
    geo = lambda P: cache.setdefault(P, geometry(depth, fine, P)) if P not in cache else cache[P]

    jobs, rows, failed, last = make_jobs(a), [], {}, None
    raw = Table(os.path.join(outdir, "timing_raw.xyz"), "PEG timing benchmark - raw data: every single run",
                ["One row per run (also failed ones, ok=0), in the order they were run; configurations (mpi,T,R) separated by 2 blank lines.",
                 "Section columns are SUMS over all solver calls (and ranks) of the job, read from the --measureTiming output."], RAW_COLS, a.meta)
    for i, (cfg, N, P) in enumerate(jobs, 1):
        if any(n <= N and p <= P for n, p in failed.get(cfg, [])):
            print("[%d/%d] skip %s T=%d R=%d N=%d P=%d (a smaller case already failed / timed out)" % ((i, len(jobs)) + cfg + (N, P)))
            continue
        if last is not None and cfg != last:
            raw.blank(2)
        last = cfg
        for rep in range(1, a.repeats + 1):
            r = run_job(a, geo, outdir, cfg, N, P, rep)
            rows.append(r); raw.row(raw_row(r))
            print("[%d/%d] %s N=%d P=%d T=%d R=%d rep %d: %s %.2f s" % (i, len(jobs), cfg[0], N, P, cfg[1], cfg[2], rep, "ok" if r["ok"] else "FAILED", r["time"]))
            if not r["ok"]:
                failed.setdefault(cfg, []).append((N, P)); break
    raw.close()
    analyse(outdir, rows, a)
    print("results in", outdir)


if __name__ == "__main__":
    main()