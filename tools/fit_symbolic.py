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
    Not one answer -- a short list, one per level of complexity, from a crude two-term expression
    up to a long one that fits better.  That list is the actual output of a symbolic regression and
    the choice among them is yours: the shortest expression whose accuracy you can live with is
    almost always the right one, because a long expression that fits a little better is usually
    fitting noise and will not survive being extrapolated.

    Complexity here is just the node count of the expression tree.

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

import argparse
import datetime
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
# ------------------------------------------------------------------------------------------------
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
                return np.where(np.isfinite(r), r, 1e6)
            sol = least_squares(resid, c0, max_nfev=budget, method='lm')
            m = float(np.mean(sol.fun ** 2))
            if np.isfinite(m) and m < best[0]:
                best = (m, sol.x)
        except Exception:
            pass
    if best[1] is None:
        return np.inf, np.zeros(k)
    return best


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
#  the search
# ------------------------------------------------------------------------------------------------
def search(X, y, w, names, rng, pop=300, gens=40, max_depth=6, parsimony=4e-3, verbose=True):
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
        mse, C = fit_constants(t2, X, y, rng, span_lo, span_hi, w)
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
        if s not in front or mse < front[s][0]:
            front[s] = (mse, C, t)

    for g in range(gens):
        scored = []
        for t in P:
            mse, C, t2 = score(t)
            fit = mse / yref + parsimony * size(t2)
            scored.append((fit, mse, C, t2))
            if np.isfinite(mse):
                remember(t2, mse, C)
        scored.sort(key=lambda z: z[0])
        if verbose and (g % max(1, gens // 8) == 0 or g == gens - 1):
            print("  generation %3d   best error %.4f%%   size %d"
                  % (g, 100.0 * np.sqrt(scored[0][1] / max(len(w), 1) * len(w)), size(scored[0][3])))

        keep = [z[3] for z in scored[:max(2, pop // 10)]]
        newP = list(keep)
        while len(newP) < pop:
            def pick():
                cand = [scored[int(rng.integers(len(scored)))] for _ in range(3)]
                return min(cand, key=lambda z: z[0])[3]
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
                    help="node count of the expression to write (default: the knee of the list)")
    ap.add_argument("--units", default="per_hour", choices=["per_second", "per_hour"])
    ap.add_argument("--loss", default="relative", choices=["relative", "absolute"],
                    help="relative (default) judges every sample by its own size, so the fit is "
                         "as good at the bottom of the range as the top. Use absolute only if you "
                         "genuinely do not care about the small rates.")
    ap.add_argument("--rate-name", default="growth", help="what to call the rate in the .sym file")
    ap.add_argument("--pop", type=int, default=300)
    ap.add_argument("--gens", type=int, default=40)
    ap.add_argument("--depth", type=int, default=6,
                    help="deepest expression tree. A product of two saturating ratios needs 5.")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

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
    front = search(X, y, w, ins, rng, pop=args.pop, gens=args.gens, max_depth=args.depth)

    def report(t, C):
        """Typical and worst-case error, both as percentages, judged the way the user asked."""
        pred = evaluate(t, X, C)
        if args.loss == "absolute":
            d = np.abs(pred - y) / max(np.max(np.abs(y)), 1e-300)
        else:
            d = np.abs(pred - y) * w
        return 100.0 * float(np.sqrt(np.mean(d ** 2))), 100.0 * float(np.max(d))

    print("\nthe expressions it found, shortest first:\n")
    print("  %-6s %-9s %-9s %-7s %s" % ("nodes", "typical", "worst", "gain", "expression"))
    rows, bestsofar, prev = [], np.inf, None
    for sz in sorted(front):
        mse, C, t = front[sz]
        if mse >= bestsofar:      # a longer expression that fits no better is not interesting
            continue
        bestsofar = mse
        typ, wst = report(t, C)
        # gain: how much the error fell per extra node, the standard way of reading this list
        gain = 0.0 if prev is None else (np.log(prev[1]) - np.log(mse)) / max(sz - prev[0], 1)
        rows.append([sz, typ, wst, gain, to_text(t, ins, C), t, C])
        prev = (sz, mse)
        print("  %-6d %-9s %-9s %-7.2f %s"
              % (sz, "%.2f%%" % typ, "%.2f%%" % wst, gain, to_text(t, ins, C)))

    if not rows:
        sys.exit("\nnothing usable was found; try more generations or a larger population")

    # Default to the biggest gain per node -- the point where buying more complexity stops
    # paying for itself.  This is a heuristic and it is often not the one you want, which is why
    # the whole list is printed and --pick exists.
    knee = max(rows, key=lambda r: r[3])
    if args.pick:
        knee = next((r for r in rows if r[0] == args.pick), knee)

    print("\npicked: %d nodes, typical error %.2f%%, worst %.2f%%" % (knee[0], knee[1], knee[2]))
    print("   %s = %s" % (args.rate_name, knee[4]))
    print("\nThis pick is the steepest gain per node, which is a rule of thumb, not an answer.")
    print("Read the list yourself: the shortest expression whose accuracy you can live with is")
    print("usually the right one, and --pick <nodes> takes a different one.")

    if args.out:
        lo, hi = X.min(axis=0), X.max(axis=0)
        with open(args.out, "w") as f:
            # 'provenance' lines are kept by the loader and printed in the run log, so a
            # result carries the fit that produced it.  '#' lines are plain comments and are
            # dropped -- which is why the explanation below is a comment and these two are not.
            f.write("provenance discovered by fit_symbolic.py on %s from %s\n"
                    % (datetime.date.today().isoformat(), args.data))
            f.write("provenance %d nodes; typical error %.2f%%, worst %.2f%%, over %d samples (%s loss)\n"
                    % (knee[0], knee[1], knee[2], X.shape[0], args.loss))
            f.write("# THE SUBSTRATE LINES ARE NOT WRITTEN FOR YOU. Add them by hand from the\n"
                    "# reaction stoichiometry, as multiples of this rate, so they stay in exact\n"
                    "# ratio -- a separate fit per species would not.\n")
            f.write("units  %s\n" % args.units)
            f.write("vars   %s\n\n" % " ".join(ins))
            f.write("rate   %s = %s\n\n" % (args.rate_name, knee[4]))
            for nm, a, b in zip(ins, lo, hi):
                f.write("range  %-10s %.6g %.6g\n" % (nm, a, b))
        print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
