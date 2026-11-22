"""Scalarization methods for continuous multi-criteria problems [F: U -> R^N] -> min.

S1 – linear scalarization:            min  sum_i lambda_i F_i(x)
S2 – epsilon-constraint method:       min  F_j(x)  s.t.  F_i(x) <= a_i, i != j
S3 – distance from a dominating point: min  || F(x) - r ||_{lambda, p}

Every scalarized problem is solved with scipy.optimize.minimize (SLSQP – the counterpart of Matlab's fmincon,
constraints as scipy Bounds / NonlinearConstraint) from Latin hypercube starting points (scipy.stats.qmc);
the solution for the previous parameter value is added as a warm start.
The criteria are entered as formulas in x1..xn and parsed with sympy, which also gives the exact Jacobian.
"""
import itertools
import math
import time
from dataclasses import dataclass, field

import numpy as np
import sympy as sp
from scipy.interpolate import PchipInterpolator, make_interp_spline
from scipy.optimize import Bounds, NonlinearConstraint, minimize
from scipy.spatial import Delaunay
from scipy.stats import qmc, uniform_direction

from owd.algorithms.with_filtration import owd_z_filtracja

DOMAINS = {
    "square": "Kwadrat / sześcian [-1, 1]ⁿ",
    "ball": "Koło / kula ‖x‖ ≤ 1",
}
METHODS = {
    "S1": "S1 – funkcja liniowa",
    "S2": "S2 – ε-ograniczenia",
    "S3": "S3 – odległość od punktu",
}
PRESETS = {
    2: {
        "Wypukły (domyślny)": ["x1**2 + x2**2", "(x1 - 1)**2 + (x2 - 1)**2"],
        "Niewypukły front (Fonseca–Fleming)": ["1 - exp(-((x1 - sqrt(2)/2)**2 + (x2 - sqrt(2)/2)**2))",
                                               "1 - exp(-((x1 + sqrt(2)/2)**2 + (x2 + sqrt(2)/2)**2))"],
        "Wielomiany 3. stopnia": ["x1**2 + x2**2 + x1*x2", "(x1 - 1)**2 + x2**3 + x2"],
    },
    3: {
        "Mieszany – jedno liniowe (domyślny)": ["x1**2 + x2**2 + x3**2", "(x1 - 1)**2 + x2**2 + (x3 + 1)**2",
                                                "x1 + x2 + x3"],
        "Trzy paraboloidy": ["(x1 - 1)**2 + x2**2 + x3**2", "x1**2 + (x2 - 1)**2 + x3**2",
                             "x1**2 + x2**2 + (x3 - 1)**2"],
    },
}
DEFAULT_PRESET = {2: "Wypukły (domyślny)", 3: "Mieszany – jedno liniowe (domyślny)"}

GAP_FACTOR = 4.0  # a spline / surface is not drawn across gaps longer than GAP_FACTOR × the median spacing


# ------------------------------------------------------------- formulas
def variables(n):
    return sp.symbols(f"x1:{n + 1}", real=True)


def parse_function(text, n):
    """Parse a formula in x1..xn into a sympy expression; raises ValueError with a message for the user."""
    symbols = variables(n)
    if not text or not text.strip():
        raise ValueError("pusta formuła")
    try:
        expr = sp.sympify(text, locals={s.name: s for s in symbols})
    except Exception as error:  # sympify raises many kinds of errors for invalid input
        raise ValueError(f"niepoprawna formuła ({type(error).__name__})") from None
    if not isinstance(expr, sp.Expr):
        raise ValueError("formuła nie jest wyrażeniem liczbowym")
    unknown = expr.free_symbols - set(symbols)
    if unknown:
        raise ValueError(f"nieznane symbole: {', '.join(sorted(map(str, unknown)))} "
                         f"(dozwolone: {', '.join(s.name for s in symbols)})")
    undefined = expr.atoms(sp.core.function.AppliedUndef)
    if undefined:
        raise ValueError(f"nieznane funkcje: {', '.join(sorted(map(str, undefined)))}")
    if expr.has(sp.I):
        raise ValueError("funkcja musi mieć wartości rzeczywiste")
    return expr


def is_nonlinear(expr, n):
    """A function is affine if and only if its Hessian is the zero matrix."""
    return sp.hessian(expr, variables(n)).applyfunc(sp.simplify).is_zero_matrix is not True


def to_latex(expr):
    return sp.latex(expr)


# -------------------------------------------------------------- problem
class Problem:
    """Vector function F = (F_1, ..., F_N) of x in U ⊂ R^n with its exact Jacobian."""

    def __init__(self, expressions, n_vars, domain):
        symbols = variables(n_vars)
        vector = sp.Matrix(expressions)
        self.expressions = list(expressions)
        self.n_vars = n_vars
        self.n_criteria = len(expressions)
        self.domain = domain
        self._values = sp.lambdify([symbols], vector, "numpy")
        self._jacobian = sp.lambdify([symbols], vector.jacobian(symbols), "numpy")

    def values(self, x):
        return np.asarray(self._values(x), dtype=float).reshape(self.n_criteria)

    def jacobian(self, x):
        return np.asarray(self._jacobian(x), dtype=float).reshape(self.n_criteria, self.n_vars)

    def random_points(self, count, rng):
        """Starting points of the multistart spread over U: Latin hypercube in the square, uniform in the ball."""
        if self.domain == "square":
            return qmc.scale(qmc.LatinHypercube(d=self.n_vars, rng=rng).random(count), -1.0, 1.0)
        directions = uniform_direction(self.n_vars).rvs(count, random_state=rng).reshape(count, self.n_vars)
        return directions * rng.uniform(size=(count, 1)) ** (1.0 / self.n_vars)


@dataclass
class Run:
    """One scalarized problem: its parameters and the solution (None if no feasible solution was found)."""
    params: dict
    grid: tuple                       # position on the parameter grid – used to triangulate 3D surfaces
    x: np.ndarray | None = None
    F: np.ndarray | None = None
    objective: float | None = None    # value of the scalarizing function S(x) at the solution


@dataclass
class MethodResult:
    method: str
    runs: list
    time_s: float
    nondominated: np.ndarray = field(default=None)  # mask over the successful runs

    @property
    def solved(self):
        return [run for run in self.runs if run.F is not None]


def _minimize(problem, fun, jac, constraints, starts, extra_start=None):
    """Minimize fun(z) over z = (x, [t]), x ∈ U, from every start; return the best solution found or None.

    `constraints` are scipy NonlinearConstraint objects; `extra_start(x)` gives the initial value of the auxiliary
    variable t (Chebyshev metric). Runs that SLSQP does not report as successful (e.g. infeasible) are skipped.
    """
    n = problem.n_vars
    n_extra = 0 if extra_start is None else 1
    pad = np.zeros(n_extra)
    bounds = Bounds(np.r_[np.full(n, -1.0), np.full(n_extra, -np.inf)], np.r_[np.ones(n), np.full(n_extra, np.inf)])
    constraints = list(constraints)
    if problem.domain == "ball":
        constraints.append(NonlinearConstraint(lambda z: z[:n] @ z[:n], -np.inf, 1.0,
                                               jac=lambda z: np.r_[2.0 * z[:n], pad]))
    best = None
    for x0 in starts:
        z0 = x0 if extra_start is None else np.append(x0, extra_start(x0))
        with np.errstate(all="ignore"):
            try:
                result = minimize(fun, z0, jac=jac, method="SLSQP", bounds=bounds, constraints=constraints,
                                  options={"maxiter": 300, "ftol": 1e-10})
            except (ValueError, ArithmeticError):
                continue
        if result.success and np.all(np.isfinite(result.x)) and (best is None or result.fun < best.fun):
            best = result
    return None if best is None else best.x


def _with_previous(starts, previous):
    return starts if previous is None else [previous, *starts]


def _store(problem, run, z):
    if z is not None:
        run.x = z[:problem.n_vars]
        run.F = problem.values(run.x)
        if not np.all(np.isfinite(run.F)):
            run.x = run.F = None


# ---------------------------------------------------------------- grids
def grid_size(step):
    """Number of intervals on [0, 1] – the step is rounded so that it divides 1."""
    return max(1, round(1.0 / step))


def simplex_grid(m, step):
    """Weights lambda with lambda_i > 0 and sum lambda_i = 1, with the given step (interior of the simplex)."""
    k = grid_size(step)
    for counts in itertools.product(range(1, k), repeat=m - 1):
        if sum(counts) < k:
            yield counts, np.array([c / k for c in counts] + [(k - sum(counts)) / k])


def box_grid(m, step):
    """Normalized levels z ∈ [0, 1]^m with the given step."""
    k = grid_size(step)
    for counts in itertools.product(range(k + 1), repeat=m):
        yield counts, np.array(counts, dtype=float) / k


def n_runs(method, m, step):
    k = grid_size(step)
    if method == "S2":
        return (k + 1) ** (m - 1)
    return math.comb(k - 1, m - 1)


# ------------------------------------------------------------- payoff table
def payoff_table(problem, starts):
    """Row i = F at the minimizer of F_i (ties broken by the sum of the other criteria).

    Returns (table, ideal, nadir); the nadir is estimated from the table, as usual.
    """
    m = problem.n_criteria
    table = np.empty((m, m))
    for i in range(m):
        z = _minimize(problem, lambda x, i=i: problem.values(x)[i], lambda x, i=i: problem.jacobian(x)[i], [],
                      starts)
        if z is None:
            raise RuntimeError(f"Nie udało się zminimalizować kryterium {i + 1}.")
        best_i = problem.values(z)[i]
        others = [j for j in range(m) if j != i]
        # Second stage: among the minimizers of F_i take one that is nondominated
        level = best_i + 1e-7 * (1.0 + abs(best_i))
        constraint = NonlinearConstraint(lambda x, i=i: problem.values(x)[i], -np.inf, level,
                                         jac=lambda x, i=i: problem.jacobian(x)[i][None, :])
        refined = _minimize(problem, lambda x: problem.values(x)[others].sum(),
                            lambda x: problem.jacobian(x)[others].sum(axis=0), [constraint], [z, *starts])
        table[i] = problem.values(z if refined is None else refined)
    return table, np.diag(table).copy(), table.max(axis=0)


# --------------------------------------------------------------- methods
def linear_scalarization(problem, step, starts, progress=None):
    """S1: min sum lambda_i F_i(x) for lambda on the interior of the simplex."""
    runs, previous = [], None
    grid = list(simplex_grid(problem.n_criteria, step))
    for index, (counts, weights) in enumerate(grid):
        run = Run({f"λ{i + 1}": w for i, w in enumerate(weights)}, counts)
        z = _minimize(problem, lambda x: weights @ problem.values(x), lambda x: weights @ problem.jacobian(x), [],
                      _with_previous(starts, previous))
        _store(problem, run, z)
        previous = run.x if run.x is not None else previous
        runs.append(run)
        if progress:
            progress((index + 1) / len(grid))
    return runs


def epsilon_constraint(problem, step, j, ideal, nadir, starts, progress=None):
    """S2: min F_j(x) s.t. F_i(x) <= a_i, a_i = ideal_i + z_i (nadir_i - ideal_i), z_i ∈ [0, 1], i != j."""
    others = [i for i in range(problem.n_criteria) if i != j]
    runs, previous = [], None
    grid = list(box_grid(len(others), step))
    for index, (counts, levels) in enumerate(grid):
        bounds = ideal[others] + levels * (nadir[others] - ideal[others])
        params = {f"z{i + 1}": z for i, z in zip(others, levels)}
        params.update({f"a{i + 1}": a for i, a in zip(others, bounds)})
        run = Run(params, counts)
        constraint = NonlinearConstraint(lambda x: problem.values(x)[others], -np.inf, bounds,
                                         jac=lambda x: problem.jacobian(x)[others])
        z = _minimize(problem, lambda x: problem.values(x)[j], lambda x: problem.jacobian(x)[j], [constraint],
                      _with_previous(starts, previous))
        _store(problem, run, z)
        previous = run.x if run.x is not None else previous
        runs.append(run)
        if progress:
            progress((index + 1) / len(grid))
    return runs


def distance_scalarization(problem, step, reference, p, starts, progress=None):
    """S3: min || lambda ∘ (F(x) - r) ||_p (numpy.linalg.norm; the gradient is approximated by scipy).

    For p = inf the max is not differentiable, which SLSQP handles poorly, so it is replaced by an auxiliary
    variable t: min t s.t. lambda_i (F_i(x) - r_i) <= t.
    """
    n, m = problem.n_vars, problem.n_criteria
    runs, previous = [], None
    grid = list(simplex_grid(m, step))
    for index, (counts, weights) in enumerate(grid):
        run = Run({f"λ{i + 1}": w for i, w in enumerate(weights)}, counts)
        starts_now = _with_previous(starts, previous)
        if math.isinf(p):
            constraint = NonlinearConstraint(
                lambda z: weights * (problem.values(z[:n]) - reference) - z[-1], -np.inf, 0.0,
                jac=lambda z: np.column_stack([weights[:, None] * problem.jacobian(z[:n]), -np.ones(m)]))
            z = _minimize(problem, lambda z: z[-1], lambda z: np.r_[np.zeros(n), 1.0], [constraint], starts_now,
                          extra_start=lambda x: np.max(weights * (problem.values(x) - reference)))
        else:
            z = _minimize(problem, lambda x: np.linalg.norm(weights * (problem.values(x) - reference), ord=p),
                          "2-point", [], starts_now)
        _store(problem, run, z)
        previous = run.x if run.x is not None else previous
        runs.append(run)
        if progress:
            progress((index + 1) / len(grid))
    return runs


def _rounded(F, scale):
    """F rounded to 1e-6 × scale, so solver noise does not decide between (weakly) equal points."""
    tolerance = 1e-6 * scale if scale > 0 else 1e-9
    return np.round(np.asarray(F, dtype=float) / tolerance) * tolerance


def nondominated_mask(F, scale):
    """Mask of nondominated rows of F (our Pareto filtering algorithm)."""
    if len(F) == 0:
        return np.zeros(0, dtype=bool)
    rounded = [tuple(row) for row in _rounded(F, scale).tolist()]
    nondominated, _ = owd_z_filtracja(rounded)
    nondominated = set(nondominated)
    return np.array([row in nondominated for row in rounded])


def first_occurrence(F, scale):
    """Index of the first row equal to each row of F – many parameter values may give the same point."""
    if len(F) == 0:
        return np.zeros(0, dtype=int)
    _, first, inverse = np.unique(_rounded(F, scale), axis=0, return_index=True, return_inverse=True)
    return first[inverse.ravel()]


def scalarized_value(method, run, settings):
    """Value of the scalarizing function that was minimized in this run."""
    if method == "S2":
        return float(run.F[settings["eps_criterion"]])
    weights = np.array([value for name, value in run.params.items() if name.startswith("λ")])
    if method == "S1":
        return float(weights @ run.F)
    return float(np.linalg.norm(weights * (run.F - settings["reference"]), ord=settings["p"]))


def distance_to_ideal(F, ideal, nadir, p):
    """Distance of each row of F from the ideal point after normalizing every criterion by (nadir - ideal)."""
    ranges = np.where(nadir - ideal > 0, nadir - ideal, 1.0)
    normalized = np.abs((np.atleast_2d(F) - ideal) / ranges)
    return np.linalg.norm(normalized, ord=p, axis=1)


def run_method(method, problem, settings, ideal, nadir, starts, progress=None):
    started = time.perf_counter()
    if method == "S1":
        runs = linear_scalarization(problem, settings["step"], starts, progress)
    elif method == "S2":
        runs = epsilon_constraint(problem, settings["step"], settings["eps_criterion"], ideal, nadir, starts,
                                  progress)
    else:
        runs = distance_scalarization(problem, settings["step"], settings["reference"], settings["p"], starts,
                                      progress)
    result = MethodResult(method, runs, time.perf_counter() - started)
    solved = result.solved
    for run in solved:
        run.objective = scalarized_value(method, run, settings)
    result.nondominated = nondominated_mask(np.array([run.F for run in solved]),
                                            float(np.max(nadir - ideal)))
    return result


# --------------------------------------------------------------- splines
def _split_at_gaps(points):
    """Split an ordered polyline where a segment is much longer than the median one."""
    if len(points) < 4:
        return [points]
    gaps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    limit = GAP_FACTOR * max(np.median(gaps), 1e-12)
    cuts = np.flatnonzero(gaps > limit) + 1
    return np.split(points, cuts)


def _dedupe(points, tolerance):
    """Mask without consecutive repeated points (a spline needs distinct nodes)."""
    return np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1) > tolerance]


def front_order(F):
    """Order of the points of a 2D front – by the first criterion."""
    return np.lexsort((F[:, 1], F[:, 0]))


def spline_segments(points, samples=200):
    """Interpolating cubic spline through ordered points (chord-length parametrization), split at gaps.

    Returns a list of arrays (samples × dim) – one per continuous piece.
    """
    points = np.asarray(points, dtype=float)
    scale = float(np.ptp(points, axis=0).max()) if len(points) else 0.0
    points = points[_dedupe(points, 1e-9 * max(scale, 1.0))]
    pieces = []
    for piece in _split_at_gaps(points):
        if len(piece) < 2:
            pieces.append(piece)
            continue
        t = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(piece, axis=0), axis=1))])
        spline = make_interp_spline(t, piece, k=min(3, len(piece) - 1))
        pieces.append(spline(np.linspace(0.0, t[-1], samples)))
    return pieces


def front_spline_segments(F, samples=200):
    """Monotone (PCHIP) spline F2(F1) through a 2D front – no overshoot, the front stays monotone."""
    F = np.asarray(F, dtype=float)[front_order(np.asarray(F, dtype=float))]
    scale = float(np.ptp(F, axis=0).max()) if len(F) else 0.0
    F = F[_dedupe(F, 1e-9 * max(scale, 1.0))]
    pieces = []
    for piece in _split_at_gaps(F):
        _, unique = np.unique(piece[:, 0], return_index=True)
        piece = piece[np.sort(unique)]
        if len(piece) < 2:
            pieces.append(piece)
            continue
        f1 = np.linspace(piece[0, 0], piece[-1, 0], samples)
        pieces.append(np.column_stack([f1, PchipInterpolator(piece[:, 0], piece[:, 1])(f1)]))
    return pieces


def surface_triangles(grid, values):
    """Triangles of a 3D surface: Delaunay triangulation of the parameter grid, without triangles over gaps."""
    grid = np.asarray(grid, dtype=float)
    values = np.asarray(values, dtype=float)
    if len(grid) < 3:
        return np.zeros((0, 3), dtype=int)
    try:
        triangles = Delaunay(grid).simplices
    except Exception:  # collinear parameter points (e.g. one row of the grid)
        return np.zeros((0, 3), dtype=int)
    edges = np.stack([np.linalg.norm(values[triangles[:, a]] - values[triangles[:, b]], axis=1)
                      for a, b in ((0, 1), (1, 2), (0, 2))], axis=1)
    longest = edges.max(axis=1)
    limit = GAP_FACTOR * max(np.median(longest), 1e-12)
    return triangles[longest <= limit]
