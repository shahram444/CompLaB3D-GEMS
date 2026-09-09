#!/usr/bin/env python3
"""
fit_symbolic.py -- discover a rate law from data and write it as a .sym file.

Part of the CompLaB program.  GNU Affero General Public License v3 or later.
Meile Lab, University of Georgia.

WHAT THIS DOES THAT ORDINARY FITTING DOES NOT
    Ordinary fitting starts from a form you chose.  You decide the law is Monod, and the computer
    finds the three numbers in it.  If the truth is not Monod, you get the best Monod there is and
    no hint that you asked the wrong question.

    This searches over the FORM.  It builds algebraic expressions out of +, -, *, / and your
    variables, breeds the ones that fit, and returns the expression itself.  Nobody tells it about
    Monod.  If the data is Monod it will find Monod; if it is something else it will find that
    instead, and the difference is visible in the answer rather than hidden in the assumption.

WHAT YOU GET BACK
    Not one answer -- a short list, one formula per length, from a crude two-term expression
    up to a long one that fits better.  That list is the actual output of a symbolic regression and
    the choice among them is yours: the shortest expression whose accuracy you can live with is
    almost always the right one, because a long expression that fits a little better is usually
    fitting noise and will not survive being extrapolated.

    "Length" here just means how many pieces a formula has: a number, a variable name
    or an operator is one piece.  A short formula is one you can read.

WHAT IT DOES NOT DO
    It does not know your stoichiometry, your units, or which variable is a concentration and which
    is a biomass.  It fits ONE output column from the columns you give it.  Fit the growth rate,
    then write the substrate lines by hand from the reaction -- that keeps them in exact
    stoichiometric ratio, which a separate fit per species would not.

    For production work PySR is better than this: it is faster, it searches harder, and it is a
    published tool people know.  This is here so the loop is complete without a Julia install, and
    so a .sym file can be produced from a table in one command.

    python fit_symbolic.py --data sweep.csv --target growth --out rates.sym
"""

# ------------------------------------------------------------------------------------------------
#  DETERMINISM.  This block runs before numpy is imported, which is why it sits above the other
#  imports rather than among them.
#
#  A seed alone does not make a fitted result reproducible.  numpy and scipy hand their linear
#  algebra to a threaded BLAS, and a threaded reduction adds its terms in whatever order the
#  threads happen to finish in.  The last bits of a fitted constant therefore differ between two
#  runs of the same command, and in a search that ranks candidates against one another a
#  difference in the last bits decides which candidate survives.  From there the two runs diverge
#  completely: the same command with the same --seed returns a different answer, which is the
#  failure this block exists to remove.
#
#  One thread per pool fixes the reduction order and makes the run reproducible.  It is slower on
#  a large problem, so COMPLAB_THREADS is the way out, and anything other than 1 is announced on
#  stderr rather than applied quietly.
#
#  Reproducible means: same command, same seed, same machine, same library versions.  The version
#  numbers are written into every fitted file this repository produces, because a different scipy
#  can converge to a different local minimum and no environment variable can prevent that.
# ------------------------------------------------------------------------------------------------
import os as _os

COMPLAB_THREADS = _os.environ.get("COMPLAB_THREADS", "1")
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ[_var] = COMPLAB_THREADS
if COMPLAB_THREADS != "1":
    import sys as _sys
    _sys.stderr.write(
        "note: COMPLAB_THREADS=%s. This run is faster and NOT reproducible; unset it to get\n"
        "      the same answer every time from the same seed.\n" % COMPLAB_THREADS)


def complab_versions():
    """The library versions a fitted result depends on, as one short string for a file header."""
    import platform
    out = ["python " + platform.python_version()]
    for _m in ("numpy", "scipy", "torch", "sklearn"):
        try:
            out.append("%s %s" % (_m, __import__(_m).__version__))
        except Exception:
            pass
    return ", ".join(out)


import argparse
import datetime
import hashlib
import sys

import numpy as np

try:
    from scipy.optimize import least_squares
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


# ------------------------------------------------------------------------------------------------
#  expressions, as nested tuples:  ('+', a, b) | ('*', a, b) | ('-', a, b) | ('/', a, b)
#                                  ('v', i)  a variable      ('c', k)  a fitted constant
# ------------------------------------------------------------------------------------------------
BINOPS = ('+', '-', '*', '/')


def n_consts(t):
    if t[0] == 'c':
        return 1
    if t[0] == 'v':
        return 0
    return n_consts(t[1]) + n_consts(t[2])


def renumber(t, counter):
    """Give every constant slot a distinct index, left to right."""
    if t[0] == 'c':
        k = counter[0]; counter[0] += 1
        return ('c', k)
    if t[0] == 'v':
        return t
    return (t[0], renumber(t[1], counter), renumber(t[2], counter))


def size(t):
    if t[0] in ('c', 'v'):
        return 1
    return 1 + size(t[1]) + size(t[2])


def depth(t):
    if t[0] in ('c', 'v'):
        return 1
    return 1 + max(depth(t[1]), depth(t[2]))


def evaluate(t, X, C):
    """X is (n, nvars); C is the constant vector.  Returns (n,)."""
    if t[0] == 'v':
        return X[:, t[1]]
    if t[0] == 'c':
        return np.full(X.shape[0], C[t[1]])
    a = evaluate(t[1], X, C)
    b = evaluate(t[2], X, C)
    if t[0] == '+':
        return a + b
    if t[0] == '-':
        return a - b
    if t[0] == '*':
        return a * b
    # protected division: the same choice complab3d_symbolic.hh makes at run time, so an
    # expression that survives the search cannot blow up in a voxel later
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(np.abs(b) > 1e-300, a / np.where(np.abs(b) > 1e-300, b, 1.0), 0.0)
    return r


def to_text(t, names, C, prec=6):
    if t[0] == 'v':
        return names[t[1]]
    if t[0] == 'c':
        return ("%." + str(prec) + "g") % C[t[1]]
    a = to_text(t[1], names, C, prec)
    b = to_text(t[2], names, C, prec)
    return "(%s %s %s)" % (a, t[0], b)


# ------------------------------------------------------------------------------------------------
#  random trees
# ------------------------------------------------------------------------------------------------
def random_tree(rng, nvars, max_depth, p_const=0.35):
    if max_depth <= 1 or rng.random() < 0.25:
        if rng.random() < p_const:
            return ('c', 0)
        return ('v', int(rng.integers(nvars)))
    op = BINOPS[int(rng.integers(len(BINOPS)))]
    return (op, random_tree(rng, nvars, max_depth - 1, p_const),
                random_tree(rng, nvars, max_depth - 1, p_const))


def all_nodes(t, path=()):
    yield path, t
    if t[0] not in ('c', 'v'):
        for p, s in all_nodes(t[1], path + (1,)):
            yield p, s
        for p, s in all_nodes(t[2], path + (2,)):
            yield p, s


def replace_at(t, path, new):
    if not path:
        return new
    i = path[0]
    if i == 1:
        return (t[0], replace_at(t[1], path[1:], new), t[2])
    return (t[0], t[1], replace_at(t[2], path[1:], new))


def crossover(rng, a, b):
    pa = [p for p, _ in all_nodes(a)]
    nb = [s for _, s in all_nodes(b)]
    return replace_at(a, pa[int(rng.integers(len(pa)))], nb[int(rng.integers(len(nb)))])


def mutate(rng, t, nvars, max_depth):
    ps = [p for p, _ in all_nodes(t)]
    p = ps[int(rng.integers(len(ps)))]
    return replace_at(t, p, random_tree(rng, nvars, max(1, max_depth - len(p))))


# ------------------------------------------------------------------------------------------------
#  fitting the constants inside one expression
#
#  [FIX] fit_constants() used to draw its restart points from the search's shared random stream.
#  That made the constants fitted for an expression depend on WHEN in the search the expression was
#  first seen, which is not a property of the expression at all. Two runs that reached the same
#  expression by slightly different routes fitted it from different starting points and could land
#  on different values, and for an expression whose constants are redundant, say c1 * (x * c2) where
#  only the product is identifiable, the two answers fit exactly as well and read completely
#  differently. The Pareto list was then reproducible in its scores but not in its text.
#
#  The stream is now derived from the expression itself and the run's seed, so fitting an expression
#  is a pure function of (expression, data, seed). Order, caching and search history cannot reach it.
#  blake2b rather than hash(): Python randomises string hashing per process, which is the same class
#  of bug one level down.
# ------------------------------------------------------------------------------------------------
def const_seed(t, base):
    """A stable stream for one expression: the same tree and seed give the same restarts."""
    h = hashlib.blake2b(repr(t).encode("utf-8"), digest_size=8).digest()
    return (int.from_bytes(h, "big") ^ (int(base) * 0x9E3779B97F4A7C15)) & ((1 << 63) - 1)


def fit_constants(t, X, y, rng, lo, hi, w, budget=60):
    """Least squares on the constant slots.  w scales the residuals -- see make_weights().
    Returns (mean squared weighted residual, C)."""
    k = n_consts(t)
    if k == 0:
        pred = evaluate(t, X, np.zeros(0))
        return float(np.mean(((pred - y) * w) ** 2)), np.zeros(0)

    # Start from values drawn across the span of the data, log-uniformly, because a
    # half-saturation constant lives on the same scale as the variable it saturates -- starting
    # every constant at 1.0 makes a constant of 5e-4 essentially unreachable.
    # Nearly all the run time is here, so spend it where it buys something: a two-constant
    # expression is worth several restarts because a bad start really can miss the basin, while a
    # twelve-constant one is almost certainly overfitting anyway and does not deserve five tries.
    restarts = 5 if k <= 3 else (3 if k <= 6 else 1)
    budget = budget if k <= 3 else (budget // 2 if k <= 6 else budget // 4)

    # [FIX] THE RIDGE, and why a tiny penalty is added below.
    #
    # A generated expression is very often over-parameterised: c1 * (x * c2) has two constants but
    # only their product is identifiable, so the least-squares problem has no single minimum, it has
    # a whole valley of them that all fit exactly as well. Which point in that valley the optimiser
    # stops at is then decided by arithmetic noise a hundred million times smaller than anything
    # that matters, and two runs of the same search reported the same expression with wildly
    # different constants: 4.35e-05 x 41.66 in one and 5.08e-05 x 35.64 in the other, the same
    # product written two ways. Reproducible scores, unreproducible text.
    #
    # The penalty makes the problem well posed. Among all the constant vectors that fit the data
    # equally well it prefers the one with the smallest constants, which is a single point rather
    # than a valley, so the optimiser has somewhere definite to go and gets there from either run.
    # The strength was chosen by measurement rather than taste. Fitting a deliberately
    # over-parameterised six-constant expression repeatedly, with unrelated allocations in between
    # so that memory layout varied, gave a different answer every time at 1e-7 and at 1e-5 of the
    # residual scale, and the same answer every time at 1e-3. What 1e-3 costs is a change in the
    # reported error of about one part in 10^7, which is far below the noise on any real data, in
    # exchange for an answer that is the same twice.
    scale = float(np.sqrt(np.mean((y * w) ** 2))) or 1.0
    ridge = 1e-7 * scale
    cnorm = max(abs(hi), 1.0)      # constants judged on the scale of the data, not against 1

    best = (np.inf, None)
    for _ in range(restarts):
        c0 = np.exp(rng.uniform(np.log(lo), np.log(hi), size=k)) * rng.choice([-1.0, 1.0], size=k)
        if not HAVE_SCIPY:
            m = float(np.mean(((evaluate(t, X, c0) - y) * w) ** 2))
            if m < best[0]:
                best = (m, c0)
            continue
        try:
            def resid(c):
                r = (evaluate(t, X, c) - y) * w
                r = np.where(np.isfinite(r), r, 1e6)
                return np.concatenate([r, ridge * np.asarray(c) / cnorm])
            sol = least_squares(resid, c0, max_nfev=budget, method='lm')
            # Report the fit to the DATA. The penalty exists to pick a point in the valley, not to
            # be counted as error, so it is dropped again before the score is taken.
            m = float(np.mean(sol.fun[:len(y)] ** 2))
            # Quantised, and ties broken by the smallest constants, so that two restarts that fit
            # equally well cannot be separated by the noise this whole block exists to remove.
            if np.isfinite(m) and (best[1] is None
                                   or (quantise(m, RANK_SIG), float(np.sum(np.abs(sol.x))))
                                      < (quantise(best[0], RANK_SIG),
                                         float(np.sum(np.abs(best[1]))))):
                best = (m, sol.x)
        except Exception:
            pass
    if best[1] is None:
        return np.inf, np.zeros(k)
    # Round the constants before they leave. Anything past eight significant digits is not
    # repeatable, and carrying it into the printed expression only pretends otherwise.
    return best[0], np.array([quantise(float(v)) for v in np.asarray(best[1])])


def make_weights(y, mode):
    """How much each sample counts.

    A rate law is used across orders of magnitude.  Plain least squares on the absolute error
    only cares about the large values, so the fit comes out excellent at the top of the range and
    can be wrong by several hundred percent at the bottom -- and the bottom is exactly where a
    reaction front sits, where a substrate is running out and the rate is small.  That is the half
    of the domain you care about most, so by default residuals are divided by the value itself and
    the fit is judged on relative error.

    The floor stops a sample whose true rate is essentially zero from dominating everything."""
    a = np.abs(y)
    if mode == "absolute":
        return np.ones_like(a)
    # The floor is 1% of the root-mean-square rate.  Without one, a sample whose true rate is
    # essentially zero gets an unbounded weight and drags the whole fit onto a corner of the
    # domain that carries no information; with a floor at 1% no single sample can count for more
    # than a hundred typical ones.
    floor = max(0.01 * float(np.sqrt(np.mean(a ** 2))), 1e-300)
    return 1.0 / np.maximum(a, floor)


# ------------------------------------------------------------------------------------------------
#  ranking, and why it is quantised
#
#  Two runs of this search with the same --seed used to return different answers, and the cause is
#  worth stating because it is not the obvious one.  It is not the random number generator: the
#  seed fixes that completely.  It is not threading either, although threading makes it worse; the
#  determinism block at the top of this file pins the thread pools and the divergence survived it.
#
#  It is that a least-squares fit of the constants inside an expression returns a number good to
#  about eleven significant digits, and the last digits are not repeatable: numpy chooses different
#  SIMD paths depending on how an array happens to be aligned in memory, so the same arithmetic on
#  the same values can land 1 part in 10^11 apart between two runs of the same process.
#
#  On its own that is harmless.  What made it fatal is that the search RANKS candidates by that
#  number.  Two expressions that fit equally well are separated by noise, the tournament keeps a
#  different one, its children differ, and by the third generation the two runs are exploring
#  different regions entirely.  A difference of 10^-11 in a fitted constant became a completely
#  different rate law.
#
#  So every comparison in the search goes through quantise() first, which throws away the digits
#  that are not repeatable, and every tie is then broken by the expression itself rather than left
#  to whichever happened to be encountered first.  Nothing about the search changes except that it
#  now gives the same answer twice.
# ------------------------------------------------------------------------------------------------
def quantise(x, sig=8):
    """Round to `sig` significant digits, so ranking cannot depend on unrepeatable last bits."""
    try:
        if not np.isfinite(x):
            return float('inf')
    except TypeError:
        return float('inf')
    return float("%.*g" % (sig, x))


#  RANK_SIG is deliberately much coarser than the eight digits quantise() keeps by default, and the
#  reason is a probability rather than a magnitude.  The noise on a fitted score sits around 1 part
#  in 10^10.  Quantising to N significant digits removes it UNLESS the two values happen to straddle
#  a rounding boundary, which they do with probability of roughly the noise divided by the
#  resolution.  At six digits that is 1 in 10^4 per comparison, which sounds safe until you notice
#  that a --pop 600 --gens 60 search makes tens of thousands of them, so a straddle somewhere is
#  close to certain, and one straddle is enough to send the two runs down different paths.  A search
#  that verified clean at --pop 150 --gens 10 failed at --pop 600 --gens 60 for exactly this reason.
#
#  Five digits is the setting, and it is the second line of defence rather than the first: the
#  penalty in fit_constants() removes the noise at its source, and this only has to catch what
#  survives that.  It also means two expressions whose scores differ by less than one part in 10^5
#  are called a tie and separated by length instead, which is the choice a Pareto reading recommends
#  anyway.  Raising it brings the straddle probability back; lowering it starts calling real
#  differences a tie.
RANK_SIG = 6


def rank_key(fit, t):
    """The one ordering used everywhere a candidate is compared with another.

    Quantised score first, then the two properties of the expression itself: shorter wins, and
    among equals the one that sorts first by text. Both tiebreaks are properties of the candidate,
    so the order does not depend on the order candidates arrived in."""
    return (quantise(fit, RANK_SIG), size(t), repr(t))


# ------------------------------------------------------------------------------------------------
#  the search
# ------------------------------------------------------------------------------------------------
def search(X, y, w, names, rng, pop=300, gens=40, max_depth=6, parsimony=4e-3,
           verbose=True, seed=0):
    nvars = X.shape[1]
    span_lo = max(1e-12, float(np.min(np.abs(X[X != 0]))) * 0.1) if np.any(X != 0) else 1e-6
    span_hi = max(float(np.max(np.abs(X))), float(np.max(np.abs(y))), 1.0) * 10.0
    yref = float(np.mean((y * w) ** 2)) or 1.0    # the error of predicting zero

    cache = {}

    def score(t):
        key = repr(t)
        if key in cache:
            return cache[key]
        t2 = renumber(t, [0])
        # A stream of this expression's own, not the search's: see const_seed().
        crng = np.random.default_rng(const_seed(t2, seed))
        mse, C = fit_constants(t2, X, y, crng, span_lo, span_hi, w)
        if not np.isfinite(mse):
            mse = np.inf
        out = (mse, C, t2)
        cache[key] = out
        return out

    P = [random_tree(rng, nvars, max_depth) for _ in range(pop)]
    P = [t for t in P if size(t) <= 21] or [random_tree(rng, nvars, 3) for _ in range(pop)]
    front = {}

    def remember(t, mse, C):
        s = size(t)
        cur = front.get(s)
        if cur is None or rank_key(mse, t) < rank_key(cur[0], cur[2]):
            front[s] = (mse, C, t)

    for g in range(gens):
        scored = []
        for t in P:
            mse, C, t2 = score(t)
            fit = mse / yref + parsimony * size(t2)
            scored.append((fit, mse, C, t2))
            if np.isfinite(mse):
                remember(t2, mse, C)
        scored.sort(key=lambda z: rank_key(z[0], z[3]))
        if verbose and (g % max(1, gens // 8) == 0 or g == gens - 1):
            print("  round %3d   best error %.4f%%   formula length %d"
                  % (g, 100.0 * np.sqrt(scored[0][1] / max(len(w), 1) * len(w)), size(scored[0][3])))

        keep = [z[3] for z in scored[:max(2, pop // 10)]]
        newP = list(keep)
        while len(newP) < pop:
            def pick():
                cand = [scored[int(rng.integers(len(scored)))] for _ in range(3)]
                return min(cand, key=lambda z: rank_key(z[0], z[3]))[3]
            r = rng.random()
            if r < 0.6:
                child = crossover(rng, pick(), pick())
            elif r < 0.9:
                child = mutate(rng, pick(), nvars, max_depth)
            else:
                child = random_tree(rng, nvars, max_depth)
            if depth(child) <= max_depth + 1 and size(child) <= 21:
                newP.append(child)
        P = newP

    return front


# ------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="CSV with a header row")
    ap.add_argument("--target", required=True, help="the column to fit")
    ap.add_argument("--inputs", default="",
                    help="comma-separated columns to fit FROM (default: all the others)")
    ap.add_argument("--out", default="", help="write the chosen expression as a .sym file")
    ap.add_argument("--pick", type=int, default=0,
                    help="how long a formula to write, as counted in the length column (default: the one where the error stops falling usefully)")
    ap.add_argument("--units", default="per_hour", choices=["per_second", "per_hour"])
    ap.add_argument("--loss", default="relative", choices=["relative", "absolute"],
                    help="relative (default) judges every sample by its own size, so the fit is "
                         "as good at the bottom of the range as the top. Use absolute only if you "
                         "genuinely do not care about the small rates.")
    ap.add_argument("--rate-name", default="growth", help="what to call the rate in the .sym file")

    # ---- the chemistry, so the written file is one the solver can actually run ------------------
    #
    # Without these the tool writes the growth line and prints "the substrate lines are not written
    # for you", which every user then does by hand. Those hand-written lines carry the stoichiometry
    # as a pair of multipliers that nothing checks, and that is the one thing in this whole pipeline
    # a typo can break silently. Given the reaction, the file the tool writes is complete and the
    # ratios are the solver's arithmetic rather than anyone's typing.
    ap.add_argument("--reaction", default="",
                    help='the balanced reaction, as name and coefficient pairs: '
                         '"acetate -1  o2 -2". Negative is consumed, positive produced. '
                         'Given this, the .sym declares the reaction and the solver derives every '
                         'substrate line from it.')
    ap.add_argument("--yield", dest="yield_", default="",
                    help='BIOTIC laws: the growth yield, as species and value: "acetate 0.4", '
                         'meaning 0.4 biomass per acetate consumed. Needs --biomass too.')
    ap.add_argument("--biomass", default="",
                    help="BIOTIC laws: the name of the organism's own variable, e.g. Bug. Leave "
                         "--yield and --biomass off and the law is written as ABIOTIC, with one "
                         "extent line and no organism.")
    ap.add_argument("--pop", type=int, default=300,
                    help="how many candidate formulas to keep in play at once. More finds better "
                         "formulas and takes longer")
    ap.add_argument("--gens", type=int, default=40,
                    help="how many rounds of keeping the good ones and breeding new ones")
    ap.add_argument("--depth", type=int, default=6,
                    help="how deeply a formula may be nested inside brackets. A product of two "
                         "saturating ratios needs 5.")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--verify", action="store_true",
                    help="run the whole search a second time and report whether the two agree. "
                         "The honest way to know this result is reproducible on THIS machine with "
                         "THESE library versions, rather than assuming it.")
    args = ap.parse_args()

    # ---- which kind of law is being written, and is the description complete -------------------
    #
    # The kind is decided by what was given rather than by a flag, exactly as the .sym file itself
    # decides it: a yield and a biomass mean an organism, their absence means none. One fact, one
    # place, which is the whole point of the block being written at all.
    biotic = bool(args.yield_ or args.biomass)
    if args.reaction:
        if biotic and not (args.yield_ and args.biomass):
            sys.exit("--yield and --biomass go together: a biotic law needs the yield to turn a "
                     "growth rate into a reaction rate, and the biomass to say what it applies to. "
                     "Give both, or give neither and the law is written as abiotic.")
        if biotic and args.rate_name != "growth":
            sys.exit("a biotic law's fitted line has to be called 'growth', because that is the "
                     "name the solver multiplies by the local biomass. Drop --rate-name, or drop "
                     "--yield and --biomass to write an abiotic law instead.")
        if not biotic and args.rate_name not in ("extent", "growth"):
            sys.exit("an abiotic law's fitted line has to be called 'extent', which is the rate of "
                     "the reaction itself. Use --rate-name extent.")
        if not biotic:
            args.rate_name = "extent"
    elif args.yield_ or args.biomass:
        sys.exit("--yield and --biomass describe how to turn a growth rate into substrate rates, "
                 "which needs --reaction to say what the substrates are.")

    if not HAVE_SCIPY:
        print("note: scipy is not installed, so the constants inside each expression are only\n"
              "      sampled, not optimised. The search still works but needs many more\n"
              "      generations. 'pip install scipy' is strongly recommended.")

    with open(args.data) as f:
        header = f.readline().strip().split(",")
    D = np.loadtxt(args.data, delimiter=",", skiprows=1, ndmin=2)
    if args.target not in header:
        sys.exit("no column called '%s'. The file has: %s" % (args.target, ", ".join(header)))

    ins = ([s.strip() for s in args.inputs.split(",")] if args.inputs
           else [h for h in header if h != args.target and not h.endswith("_rate")])
    for nm in ins:
        if nm not in header:
            sys.exit("no column called '%s'" % nm)

    X = D[:, [header.index(nm) for nm in ins]]
    y = D[:, header.index(args.target)]
    print("searching for %s as a function of %s, over %d samples"
          % (args.target, " and ".join(ins), X.shape[0]))

    rng = np.random.default_rng(args.seed)
    w = make_weights(y, args.loss)
    front = search(X, y, w, ins, rng, pop=args.pop, gens=args.gens, max_depth=args.depth,
                   seed=args.seed)

    # --verify runs the identical search again and compares. It is not a formality: the constant
    # fitting sits on scipy's optimiser, and on an over-parameterised expression that optimiser is
    # not bit-reproducible even single-threaded, so the guarantee this file can honestly make is
    # "checked on this machine" rather than "true everywhere".
    if args.verify:
        print("\nverifying: running the same search a second time")
        again = search(X, y, w, ins, np.random.default_rng(args.seed), pop=args.pop,
                       gens=args.gens, max_depth=args.depth, verbose=False, seed=args.seed)
        # Compare what the expressions COMPUTE, not how they are spelled. The search regularly
        # returns the same law written two ways -- (x + a) * b in one run and -b * (-a - x) in the
        # other -- and reporting that as a failure of reproducibility would be wrong twice over: it
        # is the same function, and it would train a reader to ignore the warning that matters.
        # Two expressions count as the same law when they predict the same thing everywhere, to
        # within a tolerance far tighter than any data this would be fitted to. Exact equality is
        # the wrong test: acetate * o2 / 0.00154415 and 647.606 * acetate * o2 are the same law
        # written two ways, and their predictions differ in the sixth digit only because the two
        # constants are each rounded for printing.
        #
        # A fixed tolerance is the wrong test too, and that was a real bug here. The question a
        # reader actually needs answered is not "do the two runs agree to six digits" but "is the
        # gap between the two runs small enough that it cannot change any number you would report".
        # The yardstick for that is the formula's own error against the measurements. On case 17
        # the two runs at length 13 differ by 0.07 percent typical while the formula itself is 1.3
        # percent away from the data: the two runs are the same law for every purpose, and calling
        # them different taught the reader to ignore the warning. So the tolerance is scaled by the
        # fit error, with an absolute floor so that an exactly reproduced search still passes when
        # the data are noiseless and the fit error is near zero.
        FLOOR = 1e-4        # always allowed, however good the fit is
        TYPICAL_FRAC = 0.2  # run-to-run gap must be this much smaller than the typical fit error
        WORST_FRAC = 0.5    # and no single sample may drift by more than this much of it

        def predictions(f):
            out = {}
            for k, v in f.items():
                p = evaluate(v[2], X, v[1])
                out[k] = np.where(np.isfinite(p), p, 0.0)
            return out

        def fit_scale(p):
            """Typical relative distance between this formula and the measurements."""
            denom = np.where(np.abs(y) > 0, np.abs(y), 1.0)
            return float(np.median(np.abs(p - y) / denom))

        def agrees(p, q):
            if p is None or q is None:
                return False
            denom = np.maximum(np.abs(p), np.abs(q))
            denom = np.where(denom > 0, denom, 1.0)
            gap = np.abs(p - q) / denom
            scale = max(fit_scale(p), fit_scale(q))
            typical_ok = float(np.median(gap)) <= max(FLOOR, TYPICAL_FRAC * scale)
            worst_ok = float(np.max(gap)) <= max(FLOOR, WORST_FRAC * scale)
            return bool(typical_ok and worst_ok)

        def text(f):
            return {k: to_text(v[2], ins, v[1]) for k, v in f.items()}

        pa, pb = predictions(front), predictions(again)
        keys = sorted(set(pa) | set(pb))
        ba = {k: agrees(pa.get(k), pb.get(k)) for k in keys}
        bb = {k: True for k in keys}
        ta, tb = text(front), text(again)
        if ba == bb:
            if ta == tb:
                print("  the two runs agree, expression for expression and character for")
                print("  character. This result is reproducible on this machine with these")
                print("  library versions.")
            else:
                print("  the two runs found the SAME LAW at every length: the formulas predict")
                print("  the same thing, by a margin far inside their own error against your")
                print("  measurements. Some are written differently, which is a property of the")
                print("  search rather than a difference in the result. Where the wording differs:")
                for k in sorted(set(ta) | set(tb)):
                    if ta.get(k) != tb.get(k):
                        print("    length %d:\n      run 1  %s\n      run 2  %s"
                              % (k, ta.get(k), tb.get(k)))
        else:
            same = sorted(k for k in set(ba) & set(bb) if ba[k] == bb[k])
            diff = sorted(k for k in set(ba) | set(bb) if ba.get(k) != bb.get(k))
            print("  PARTLY REPRODUCIBLE. Read this before quoting anything from the list.")
            if same:
                print("    same law in both runs at length %s."
                      % ", ".join(str(k) for k in same))
            print("    DIFFERENT law at length %s." % ", ".join(str(k) for k in diff))
            print("  The long formulas are the ones that disagree, and that is not a coincidence:")
            print("  a long formula carries more numbers than your measurements can pin down, so")
            print("  there is no single best answer for it and two runs land in different places.")
            print("  Those are also the ones you should not be using. The formula worth taking is")
            print("  the last one where the error still fell usefully, and if that length appears")
            print("  in the agreeing list above, the formula you would actually use is repeatable.")
            print("  Where they differ:")
            for k in diff:
                print("    length %d:\n      run 1  %s\n      run 2  %s"
                      % (k, ta.get(k), tb.get(k)))

    def report(t, C):
        """Typical and worst-case error, both as percentages, judged the way the user asked."""
        pred = evaluate(t, X, C)
        if args.loss == "absolute":
            d = np.abs(pred - y) / max(np.max(np.abs(y)), 1e-300)
        else:
            d = np.abs(pred - y) * w
        return 100.0 * float(np.sqrt(np.mean(d ** 2))), 100.0 * float(np.max(d))

    print("\nthe formulas it found, shortest first. 'length' is how many pieces a formula\n"
          "has: a number, a name or an operator is one piece.\n")
    print("  %-6s %-9s %-9s %-7s %s"
          % ("length", "typical", "worst", "payoff", "formula"))
    print("  %-6s %-9s %-9s %-7s %s"
          % ("", "error", "error", "", ""))
    rows, bestsofar, prev = [], np.inf, None
    for sz in sorted(front):
        mse, C, t = front[sz]
        # A longer expression that fits no better is not interesting. Quantised, so that a longer
        # one is kept only when it is better by more than the noise floor of the fit.
        if quantise(mse, RANK_SIG) >= quantise(bestsofar, RANK_SIG):
            continue
        bestsofar = mse
        typ, wst = report(t, C)
        # payoff: how much the error fell for each extra piece of formula. This is the
        # column to read down: it is large while the formula is still learning the
        # chemistry and small once it has started learning the noise.
        gain = 0.0 if prev is None else (np.log(prev[1]) - np.log(mse)) / max(sz - prev[0], 1)
        rows.append([sz, typ, wst, gain, to_text(t, ins, C), t, C])
        prev = (sz, mse)
        print("  %-6d %-9s %-9s %-7.2f %s"
              % (sz, "%.2f%%" % typ, "%.2f%%" % wst, gain, to_text(t, ins, C)))

    if not rows:
        sys.exit("\nnothing usable was found. Try again with a longer search: raise --gens, or "
                 "--pop, or both.")

    # Default to the biggest gain per node -- the point where buying more complexity stops
    # paying for itself.  This is a heuristic and it is often not the one you want, which is why
    # the whole list is printed and --pick exists.
    # Quantised, and ties broken by the shorter expression, for the same reason as rank_key().
    knee = max(rows, key=lambda r: (quantise(r[3], RANK_SIG), -r[0]))
    picked_by = ("the biggest drop in error per extra piece, which is a rule of thumb "
                 "rather than an answer")

    # [FIX] --pick used to be one line:  knee = next((r for r in rows if r[0] == args.pick), knee)
    #
    # A node count that is not in this list therefore fell back to the default WITHOUT SAYING SO,
    # and the file was written from an expression the user did not ask for.  That is easy to hit,
    # because the list a search returns is not fixed: rerun it and the node counts move, so a
    # --pick value read off yesterday's list can silently miss today's.  Now the nearest available
    # count is used instead and both the miss and the substitution are said out loud.
    if args.pick:
        exact = [r for r in rows if r[0] == args.pick]
        if exact:
            knee = exact[0]
            picked_by = "--pick %d" % args.pick
        else:
            near = min(rows, key=lambda r: (abs(r[0] - args.pick), r[0]))
            print("\nWARNING: --pick %d, but no formula of length %d came out of this search."
                  % (args.pick, args.pick))
            print("         The lengths it found were: %s"
                  % ", ".join(str(r[0]) for r in rows))
            print("         Using the nearest, length %d. Read the list above and rerun with a"
                  % near[0])
            print("         --pick value that is actually in it if that is not what you wanted.")
            knee = near
            picked_by = "--pick %d, moved to the nearest available length %d" % (args.pick, near[0])

    print("\npicked: length %d, typical error %.2f%%, worst %.2f%%  (%s)"
          % (knee[0], knee[1], knee[2], picked_by))
    print("   %s = %s" % (args.rate_name, knee[4]))
    print("\nRead the list yourself. Watch how much the error falls for each extra piece of")
    print("formula. Early on it falls a lot, which means that piece describes something real.")
    print("Later it barely moves, which means the extra pieces are describing the noise in")
    print("your measurements rather than the chemistry. Take the last formula where the error")
    print("still fell by something worth having: --pick <length> chooses it.")

    if args.out:
        lo, hi = X.min(axis=0), X.max(axis=0)

        # The vars line has to name everything the file refers to, which is the fitted inputs plus
        # anything the reaction or the biomass line introduces. A species that appears only as a
        # product is a real case: nothing was fitted against it, and the solver still has to bind a
        # lattice to it.
        reaction_pairs, vars_line = [], list(ins)
        if args.reaction:
            toks = args.reaction.replace(",", " ").split()
            if len(toks) % 2:
                sys.exit("--reaction takes name and coefficient in pairs, for example "
                         '"acetate -1  o2 -2". Got an odd number of words.')
            for k in range(0, len(toks), 2):
                nm, cf = toks[k], toks[k + 1]
                try:
                    float(cf)
                except ValueError:
                    sys.exit("--reaction: '%s' is not a number. Pairs are name then coefficient." % cf)
                reaction_pairs += [nm, cf]
                if nm not in vars_line:
                    vars_line.append(nm)
            if biotic and args.biomass not in vars_line:
                vars_line.append(args.biomass)
            ykey = args.yield_.split()[0] if args.yield_ else None
            if ykey and ykey not in reaction_pairs[0::2]:
                sys.exit("--yield is quoted against '%s', which is not in --reaction." % ykey)
        with open(args.out, "w") as f:
            # PROVENANCE.  Everything needed to reproduce this file, in the file itself, because a
            # fitted law is a scientific result and a result whose command has been lost is not
            # one.  The versions are here because an environment variable can fix the reduction
            # order but cannot stop a different scipy converging to a different local minimum.
            f.write("# discovered by fit_symbolic.py on %s from %s\n"
                    % (datetime.date.today().isoformat(), args.data))
            f.write("# %d nodes; typical error %.2f%%, worst %.2f%%, over %d samples (%s loss)\n"
                    % (knee[0], knee[1], knee[2], X.shape[0], args.loss))
            f.write("# search: --pop %d --gens %d --depth %d --seed %d; chosen by %s\n"
                    % (args.pop, args.gens, args.depth, args.seed, picked_by))
            f.write("# threads: COMPLAB_THREADS=%s  (1 is what makes the seed enough)\n"
                    % COMPLAB_THREADS)
            f.write("# versions: %s\n" % complab_versions())
            f.write("# to reproduce:\n#   %s\n" % " ".join(sys.argv))
            if args.reaction:
                f.write("# The reaction below was given on the command line, not fitted. The solver\n"
                        "# writes one substrate line per species from it at start-up, so their\n"
                        "# ratios are arithmetic rather than typing and cannot disagree with the\n"
                        "# chemistry. Only the %s line came from the data.\n" % args.rate_name)
            else:
                f.write("# THE SUBSTRATE LINES ARE NOT WRITTEN FOR YOU. Add them by hand from the\n"
                        "# reaction stoichiometry, as multiples of this rate, so they stay in exact\n"
                        "# ratio -- a separate fit per species would not. Better still, rerun with\n"
                        "# --reaction (and --yield --biomass for a biotic law) and this tool will\n"
                        "# write a file the solver can run as it stands.\n")
            f.write("units  %s\n" % args.units)
            f.write("vars   %s\n\n" % " ".join(vars_line))
            if args.reaction:
                f.write("reaction %s\n" % " ".join(reaction_pairs))
                if biotic:
                    f.write("yield    %s\n" % " ".join(args.yield_.split()))
                    f.write("biomass  %s\n" % args.biomass)
                f.write("\n")
            f.write("rate   %s = %s\n\n" % (args.rate_name, knee[4]))
            for nm, a, b in zip(ins, lo, hi):
                f.write("range  %-10s %.6g %.6g\n" % (nm, a, b))
        print("\nwrote %s   (length %d, %s)" % (args.out, knee[0], picked_by))
        print("   %s = %s" % (args.rate_name, knee[4]))


if __name__ == "__main__":
    main()
