"""Page: scalarization methods (S1, S2, S3) for continuous problems with 2 or 3 criteria."""
import math

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from owd.scalarization import (DEFAULT_PRESET, DOMAINS, METHODS, PRESETS, Problem, distance_to_ideal, front_order,
                           front_spline_segments, is_nonlinear, n_runs, parse_function, payoff_table, run_method,
                           spline_segments, surface_triangles, to_latex, first_occurrence)

VARIANTS = {2: "(a) dwa kryteria F, U ⊂ R²", 3: "(b) trzy kryteria G, U ⊂ R³"}
LETTERS = {2: "F", 3: "G"}
P_OPTIONS = {"1": 1.0, "2": 2.0, "∞": math.inf, "inna": None}
METHOD_COLORS = {"S1": "#1f77b4", "S2": "#ff7f0e", "S3": "#2ca02c"}
COMPROMISE_COLOR = "#e6a800"
TABLE_FILTERS = {
    "all": "Wszystkie",
    "nondominated": "Niezdominowane",
    "unique": "Unikalne niezdominowane",
    "dominated": "Zdominowane",
    "failed": "Brak rozwiązania",
}
DOMAIN_LATEX = {
    "square": {2: r"U = [-1, 1]^2", 3: r"U = [-1, 1]^3"},
    "ball": {2: r"U = \{x \in \mathbb{R}^2 : x_1^2 + x_2^2 \le 1\}",
             3: r"U = \{x \in \mathbb{R}^3 : x_1^2 + x_2^2 + x_3^2 \le 1\}"},
}


# ------------------------------------------------------------------ state
def init_state():
    defaults = {
        "sc_variant": 2,
        "sc_domain": "square",
        "sc_methods": list(METHODS),
        "sc_step_2": 0.05,
        "sc_step_3": 0.1,
        "sc_p": "2",
        "sc_p_custom": 3.0,
        "sc_delta": 0.1,
        "sc_starts": 3,
        "sc_seed": 0,
        "sc_results": None,
        "sc_ideal_p": "2",
        "sc_show_dominated": True,
        "sc_show_curve": True,
    }
    for n in PRESETS:
        defaults[f"sc_preset_{n}"] = DEFAULT_PRESET[n]
        defaults[f"sc_eps_{n}"] = 0
        for i, text in enumerate(PRESETS[n][DEFAULT_PRESET[n]]):
            defaults[f"sc_expr_{n}_{i}"] = text
    for method in METHODS:
        defaults[f"sc_filter_{method}"] = "all"
    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


def load_preset(n):
    for i, text in enumerate(PRESETS[n][st.session_state[f"sc_preset_{n}"]]):
        st.session_state[f"sc_expr_{n}_{i}"] = text


def restore_defaults(n):
    st.session_state[f"sc_preset_{n}"] = DEFAULT_PRESET[n]
    load_preset(n)


# ---------------------------------------------------------------- formulas
def method_latex(method, letter, settings):
    if method == "S1":
        return (rf"S_1:\quad \min_{{x \in U}}\ \sum_{{i=1}}^{{N}} \lambda_i\, {letter}_i(x), \qquad "
                r"\lambda_i > 0,\ \ \sum_{i=1}^{N} \lambda_i = 1")
    if method == "S2":
        j = settings["eps_criterion"] + 1
        return (rf"S_2:\quad \min_{{x \in U}}\ {letter}_{{{j}}}(x) \quad \text{{p.o.}} \quad "
                rf"{letter}_i(x) \le a_i = {letter}_i^{{*}} + z_i\,({letter}_i^{{\,\mathrm{{nad}}}} - {letter}_i^{{*}}),"
                rf"\ \ i \ne {j},\ \ z_i \in [0, 1]")
    p = settings["p"]
    reference = (rf"r = {letter}^{{*}} - \delta\,({letter}^{{\,\mathrm{{nad}}}} - {letter}^{{*}}),"
                 rf"\ \ \delta = {settings['delta']:g}")
    if math.isinf(p):
        return (rf"S_3:\quad \min_{{x \in U}}\ \max_{{i}}\ \lambda_i\, |{letter}_i(x) - r_i|"
                rf"\qquad p = \infty \text{{ (Czebyszew)}},\ \ {reference}")
    return (rf"S_3:\quad \min_{{x \in U}}\ \Big(\sum_{{i=1}}^{{N}} \big(\lambda_i\, |{letter}_i(x) - r_i|\big)^{{p}}"
            rf"\Big)^{{1/p}} \qquad p = {p:g},\ \ {reference}")


def distance_latex(letter, p_label):
    return (rf"d(x) = \left\| \frac{{{letter}(x) - {letter}^{{*}}}}{{{letter}^{{\,\mathrm{{nad}}}} - {letter}^{{*}}}} "
            rf"\right\|_{{{p_label}}}")


def read_functions(n):
    """Parsed criteria and a list of errors that block the computation."""
    letter = LETTERS[n]
    expressions, errors = [], []
    for i in range(n):
        try:
            expressions.append(parse_function(st.session_state[f"sc_expr_{n}_{i}"], n))
        except ValueError as error:
            expressions.append(None)
            errors.append(f"{letter}{i + 1}: {error}")
    nonlinear = [None if e is None else is_nonlinear(e, n) for e in expressions]
    if not errors:
        if n == 2 and not all(nonlinear):
            linear = ", ".join(f"{letter}{i + 1}" for i, flag in enumerate(nonlinear) if not flag)
            errors.append(f"W wariancie (a) obie funkcje celu muszą być nieliniowe – liniowe: {linear}.")
        if n == 3 and not any(nonlinear):
            errors.append("W wariancie (b) co najmniej jedno kryterium musi być nieliniowe – wszystkie są liniowe.")
    return expressions, nonlinear, errors


# ------------------------------------------------------------------ plots
def domain_traces(domain, n):
    style = {"mode": "lines", "line": {"color": "#999999", "width": 1}, "name": "Brzeg U", "hoverinfo": "skip"}
    if n == 2:
        if domain == "square":
            return [go.Scatter(x=[-1, 1, 1, -1, -1], y=[-1, -1, 1, 1, -1], **style)]
        theta = np.linspace(0, 2 * np.pi, 200)
        return [go.Scatter(x=np.cos(theta), y=np.sin(theta), **style)]
    xs, ys, zs = [], [], []
    if domain == "square":
        corners = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)])
        for a in range(8):
            for b in range(a + 1, 8):
                if np.sum(corners[a] != corners[b]) == 1:
                    for coordinate, values in zip((xs, ys, zs), zip(corners[a], corners[b])):
                        coordinate.extend([*values, None])
    else:
        theta = np.linspace(0, 2 * np.pi, 100)
        for height in (-0.5, 0.0, 0.5):  # latitude circles
            radius = math.sqrt(1 - height ** 2)
            xs.extend([*radius * np.cos(theta), None])
            ys.extend([*radius * np.sin(theta), None])
            zs.extend([*np.full_like(theta, height), None])
        for angle in (0.0, np.pi / 2):  # meridians
            xs.extend([*np.cos(angle) * np.cos(theta), None])
            ys.extend([*np.sin(angle) * np.cos(theta), None])
            zs.extend([*np.sin(theta), None])
    return [go.Scatter3d(x=xs, y=ys, z=zs, **style)]


def run_hover(run, letter):
    params = ", ".join(f"{name}={value:.3g}" for name, value in run.params.items())
    x = ", ".join(f"{v:.4g}" for v in run.x)
    F = ", ".join(f"{v:.4g}" for v in run.F)
    return f"{params}<br>x = ({x})<br>{letter} = ({F})"


def add_method_traces(fig, result, space, n, letter, label, color, show_dominated, show_curve):
    """Points of one method in the objective (space="F") or decision (space="x") space with a spline / surface."""
    solved = result.solved
    if not solved:
        return
    values = np.array([run.F if space == "F" else run.x for run in solved])
    hover = np.array([run_hover(run, letter) for run in solved])
    mask = result.nondominated
    use_3d = n == 3
    scatter = go.Scatter3d if use_3d else go.Scatter
    size = 4 if use_3d else 7

    def points(selection, name, **marker):
        if not selection.any():
            return
        coordinates = dict(zip("xyz", values[selection].T))
        fig.add_trace(scatter(mode="markers", name=f"{name} ({int(selection.sum())})", text=hover[selection],
                              hovertemplate="%{text}<extra>%{fullData.name}</extra>", marker=marker,
                              legendgroup=label, **coordinates))

    if show_curve and mask.sum() >= 2:
        nondominated = values[mask]
        if not use_3d:
            if space == "F":
                pieces = front_spline_segments(nondominated)
            else:  # the decision-space curve follows the order of the front
                order = front_order(np.array([run.F for run in solved])[mask])
                pieces = spline_segments(nondominated[order])
            for index, piece in enumerate(pieces):
                if len(piece) >= 2:
                    fig.add_trace(go.Scatter(x=piece[:, 0], y=piece[:, 1], mode="lines", name=f"{label} – spline",
                                             line={"color": color, "width": 2}, opacity=0.7, hoverinfo="skip",
                                             legendgroup=label, showlegend=index == 0))
        else:
            grid = np.array([run.grid for run in solved])[mask]
            triangles = surface_triangles(grid, nondominated)
            if len(triangles):
                fig.add_trace(go.Mesh3d(x=nondominated[:, 0], y=nondominated[:, 1], z=nondominated[:, 2],
                                        i=triangles[:, 0], j=triangles[:, 1], k=triangles[:, 2], color=color,
                                        opacity=0.35, name=f"{label} – powierzchnia", hoverinfo="skip",
                                        legendgroup=label, showlegend=True))
    points(mask, f"{label} – niezdominowane", size=size, color=color, line={"color": "black", "width": 0.5})
    if show_dominated:
        points(~mask, f"{label} – zdominowane", size=size, color="#d62728", symbol="x")


def add_marker(fig, n, name, points, marker, hover=None):
    """Extra points (ideal point, closest solution, selected rows) on top of a figure."""
    points = np.atleast_2d(points)
    if hover is None:
        hover = [f"{name}: (" + ", ".join(f"{v:.4g}" for v in point) + ")" for point in points]
    fig.add_trace((go.Scatter3d if n == 3 else go.Scatter)(
        mode="markers", name=name, marker=marker, text=hover, hovertemplate="%{text}<extra></extra>",
        **dict(zip("xyz", points.T))))


def space_figure(results, space, n, domain, letter, show_dominated, show_curve, markers=(), height=560):
    fig = go.Figure()
    if space == "x":
        for trace in domain_traces(domain, n):
            fig.add_trace(trace)
    for method, result in results.items():
        label = method if len(results) > 1 else METHODS[method].split(" – ")[0]
        add_method_traces(fig, result, space, n, letter, label, METHOD_COLORS[method], show_dominated, show_curve)
    for name, point, symbol, color in markers:
        add_marker(fig, n, name, point, {"symbol": symbol, "size": 9 if n == 3 else 14, "color": color})

    names = [f"{letter}{i + 1}" for i in range(n)] if space == "F" else [f"x{i + 1}" for i in range(n)]
    title = "Oceny niezdominowane F_P(U)" if space == "F" else "Punkty niezdominowane P(U, F)"
    layout = {"title": {"text": title}, "height": height, "margin": {"l": 10, "r": 10, "t": 50, "b": 10},
              "legend": {"orientation": "h", "y": -0.12, "yanchor": "top", "x": 0}}
    if n == 3:
        layout["scene"] = {f"{axis}axis_title": name for axis, name in zip("xyz", names)}
        if space == "x":
            layout["scene"]["aspectmode"] = "cube"
    else:
        layout["xaxis_title"], layout["yaxis_title"] = names
        if space == "x":
            layout["yaxis"] = {"scaleanchor": "x", "scaleratio": 1}
    fig.update_layout(**layout)
    return fig


def analyze(result, stored, p):
    """Per-run data of one method: distance from the ideal point, duplicates and the solution closest to it."""
    solved_ids = [i for i, run in enumerate(result.runs) if run.F is not None]
    F = np.array([result.runs[i].F for i in solved_ids]).reshape(len(solved_ids), stored["n"])
    distances = distance_to_ideal(F, stored["ideal"], stored["nadir"], p) if len(F) else np.zeros(0)
    first = first_occurrence(F, float(np.max(stored["nadir"] - stored["ideal"])))
    unique = first == np.arange(len(F))
    nondominated = result.nondominated
    best = None
    if nondominated.any():
        best = solved_ids[int(np.argmin(np.where(nondominated, distances, np.inf)))]
    # Duplicates point to the first run that gave the same point
    duplicate_of = {solved_ids[i]: solved_ids[j] for i, j in enumerate(first) if i != j}
    return {
        "solved_ids": solved_ids,
        "distance": dict(zip(solved_ids, distances)),
        "nondominated": dict(zip(solved_ids, nondominated)),
        "unique": dict(zip(solved_ids, unique)),
        "duplicate_of": duplicate_of,
        "best": best,
    }


def result_table(result, info, n, letter, method, p_label):
    rows = []
    for run_id, run in enumerate(result.runs):
        row = {"nr": run_id, **{name: float(value) for name, value in run.params.items()}}
        if run.F is None:
            row.update({f"x{i + 1}": np.nan for i in range(n)})
            row.update({f"{letter}{i + 1}": np.nan for i in range(n)})
            row.update({f"{method}(x)": np.nan, f"d (p={p_label})": np.nan, "jak nr": ""})
            row["status"] = "brak rozwiązania"
        else:
            row.update({f"x{i + 1}": float(v) for i, v in enumerate(run.x)})
            row.update({f"{letter}{i + 1}": float(v) for i, v in enumerate(run.F)})
            row[f"{method}(x)"] = run.objective
            row[f"d (p={p_label})"] = float(info["distance"][run_id])
            row["jak nr"] = str(info["duplicate_of"].get(run_id, ""))
            if run_id == info["best"]:
                row["status"] = "★ najbliżej punktu idealnego"
            else:
                row["status"] = "niezdominowany" if info["nondominated"][run_id] else "zdominowany"
        rows.append(row)
    return pd.DataFrame(rows)


def filter_table(table, info, table_filter):
    ids = table["nr"]
    if table_filter == "nondominated":
        return table[ids.map(lambda i: info["nondominated"].get(i, False))]
    if table_filter == "unique":
        return table[ids.map(lambda i: info["nondominated"].get(i, False) and info["unique"][i])]
    if table_filter == "dominated":
        return table[ids.map(lambda i: i in info["nondominated"] and not info["nondominated"][i])]
    if table_filter == "failed":
        return table[~ids.isin(info["solved_ids"])]
    return table


def table_column_config(table, letter, method, p_label):
    config = {"nr": st.column_config.NumberColumn("nr", help="Numer zadania (wartości parametru)", format="%d"),
              "jak nr": st.column_config.TextColumn(
                  "jak nr", help="Ten sam punkt co w zadaniu o podanym numerze (puste – pierwszy raz)"),
              f"{method}(x)": st.column_config.NumberColumn(
                  format="%.5f", help="Wartość funkcji skalaryzującej w rozwiązaniu – minimalizowana wielkość"),
              f"d (p={p_label})": st.column_config.NumberColumn(
                  format="%.4f", help="Znormalizowana odległość ocen od punktu idealnego")}
    for column in table.columns:
        if column not in config and column not in ("status", "jak nr"):
            parameter = column[0] in "λz"
            config[column] = st.column_config.NumberColumn(format="%.3f" if parameter else "%.5f")
    return config


def run_label(run):
    return ", ".join(f"{name} = {value:.3g}" for name, value in run.params.items())


def vector(values):
    return "(" + ", ".join(f"{v:.4g}" for v in values) + ")"


def show_plots(results, stored, show_dominated, show_curve, key, highlights=()):
    """Objective and decision space side by side.

    `highlights`: (name, runs, color, size) – runs marked with open circles on both plots; for the solution closest
    to the ideal point (name starting with "★") a dotted line from the ideal point is drawn too.
    """
    n, letter = stored["n"], stored["letter"]
    markers = [("Punkt idealny", stored["ideal"], "diamond" if n == 3 else "star", "#9467bd")]
    if "S3" in results:
        markers.append(("Punkt odniesienia r", stored["reference"], "cross", "#8c564b"))
    figures = {"F": space_figure(results, "F", n, stored["domain"], letter, show_dominated, show_curve, markers),
               "x": space_figure(results, "x", n, stored["domain"], letter, show_dominated, show_curve)}
    for name, runs, color, size in highlights:
        if not runs:
            continue
        if name.startswith("★"):
            link = np.array([stored["ideal"], runs[0].F])
            figures["F"].add_trace((go.Scatter3d if n == 3 else go.Scatter)(
                mode="lines", line={"color": color, "width": 2, "dash": "dot"}, hoverinfo="skip", showlegend=False,
                **dict(zip("xyz", link.T))))
        hover = [f"{name}<br>{run_hover(run, letter)}" for run in runs]
        marker = {"symbol": "circle-open", "size": size * (0.6 if n == 3 else 1), "color": color,
                  "line": {"width": 3, "color": color}}
        for space in figures:
            add_marker(figures[space], n, name, [run.F if space == "F" else run.x for run in runs], marker, hover)
    objective_column, decision_column = st.columns(2, gap="medium")
    objective_column.plotly_chart(figures["F"], width="stretch", key=f"{key}_F", config={"scrollZoom": False})
    decision_column.plotly_chart(figures["x"], width="stretch", key=f"{key}_x", config={"scrollZoom": False})


# ------------------------------------------------------------------- page
init_state()

with st.sidebar:
    st.header("Problem")
    n = st.radio("Wariant", list(VARIANTS), format_func=VARIANTS.get, key="sc_variant")
    letter = LETTERS[n]
    domain = st.radio("Zbiór U", list(DOMAINS), format_func=DOMAINS.get, key="sc_domain")

    st.header("Funkcje celu")
    preset_column, reset_column = st.columns([4, 1], vertical_alignment="bottom")
    preset_column.selectbox("Zestaw funkcji", list(PRESETS[n]), key=f"sc_preset_{n}", on_change=load_preset,
                            args=(n,))
    reset_column.button("↺", on_click=restore_defaults, args=(n,), help="Przywróć domyślne funkcje")
    for i in range(n):
        st.text_input(f"{letter}{i + 1}(x)", key=f"sc_expr_{n}_{i}",
                      help="Zmienne: " + ", ".join(f"x{j + 1}" for j in range(n)) +
                           ". Potęga: ** albo ^, funkcje: exp, log, sin, cos, sqrt, Abs, …")

    st.header("Metody")
    methods = st.pills("Metody", list(METHODS), selection_mode="multi", format_func=lambda m: METHODS[m],
                       key="sc_methods", label_visibility="collapsed") or []
    step = st.number_input("Krok ρ parametrów λ / z", min_value=0.01, max_value=0.5, step=0.01, format="%.2f",
                           key=f"sc_step_{n}", help="Zaokrąglany tak, żeby 1/ρ było liczbą całkowitą.")
    if "S2" in methods:
        eps_criterion = st.selectbox("S2: minimalizowane kryterium j", range(n),
                                     format_func=lambda i: f"{letter}{i + 1}", key=f"sc_eps_{n}")
    else:
        eps_criterion = st.session_state[f"sc_eps_{n}"]
    if "S3" in methods:
        p_column, custom_column = st.columns(2, vertical_alignment="bottom")
        p_choice = p_column.selectbox("S3: metryka p", list(P_OPTIONS), key="sc_p")
        p_custom = custom_column.number_input("p", min_value=1.0, step=0.5, key="sc_p_custom",
                                              disabled=p_choice != "inna")
        st.number_input("S3: przesunięcie δ punktu r", min_value=0.0, max_value=10.0, step=0.05, format="%.2f",
                        key="sc_delta", help="r = punkt idealny − δ·(nadir − punkt idealny); r dominuje cały FP(U).")
    p = P_OPTIONS[st.session_state.sc_p] or float(st.session_state.sc_p_custom)

    st.header("Solver (SLSQP)")
    starts_column, seed_column = st.columns(2)
    starts_column.number_input("Punkty startowe", min_value=1, max_value=30, key="sc_starts",
                               help="Multistart: liczba losowych punktów startowych dla każdego zadania "
                                    "(plus rozwiązanie dla poprzedniej wartości parametru).")
    seed_column.number_input("Seed", min_value=0, max_value=2**31 - 1, key="sc_seed")

expressions, nonlinear, errors = read_functions(n)
settings = {
    "n": n, "domain": domain, "texts": tuple(st.session_state[f"sc_expr_{n}_{i}"] for i in range(n)),
    "methods": tuple(m for m in METHODS if m in methods), "step": 1.0 / max(1, round(1.0 / step)),
    "eps_criterion": int(eps_criterion), "p": p, "delta": float(st.session_state.sc_delta),
    "starts": int(st.session_state.sc_starts), "seed": int(st.session_state.sc_seed),
}
if not methods:
    errors.append("Wybierz co najmniej jedną metodę.")
for method in settings["methods"]:
    if n_runs(method, n, settings["step"]) == 0:
        errors.append(f"{METHODS[method]}: krok ρ = {settings['step']:g} jest za duży – brak wartości λ > 0.")

with st.sidebar:
    st.divider()
    total = sum(n_runs(m, n, settings["step"]) for m in settings["methods"])
    st.caption(f"Liczba zadań optymalizacji: {total} × {settings['starts'] + 1} startów.")
    compute_clicked = st.button("Oblicz", type="primary", width="stretch", disabled=bool(errors))

st.title("Skalaryzacja")
st.caption("Wyznaczanie zbiorów P(U, F) i F_P(U) metodami skalaryzacji: każda wartość parametru λ (lub z) "
           "daje jedno zadanie optymalizacji jednokryterialnej rozwiązywane metodą SLSQP (odpowiednik fmincon).")

# ---------------------------------------------------------------- problem
st.subheader("Problem")
st.latex(rf"\left[{letter}: U \to \mathbb{{R}}^{n}\right] \to \min, \qquad " + DOMAIN_LATEX[domain][n])
for i, (expr, flag) in enumerate(zip(expressions, nonlinear)):
    formula_column, badge_column = st.columns([6, 1], vertical_alignment="center")
    if expr is None:
        formula_column.markdown(f"**{letter}{i + 1}** – niepoprawna formuła: "
                                f"`{st.session_state[f'sc_expr_{n}_{i}']}`")
        badge_column.badge("błąd", color="red")
    else:
        formula_column.latex(rf"{letter}_{{{i + 1}}}(x) = {to_latex(expr)}")
        badge_column.badge("nieliniowa" if flag else "liniowa", color="green" if flag else "orange")
for error in errors:
    st.error(error)

# ---------------------------------------------------------------- compute
if compute_clicked and not errors:
    problem = Problem(expressions, n, domain)
    rng = np.random.default_rng(settings["seed"])
    starts = list(problem.random_points(settings["starts"], rng))
    progress = st.progress(0.0, text="Tablica wypłat…")
    try:
        table, ideal, nadir = payoff_table(problem, starts)
    except RuntimeError as error:
        progress.empty()
        st.error(str(error))
        st.stop()
    reference = ideal - settings["delta"] * (nadir - ideal)
    solver_settings = {**settings, "reference": reference}
    results = {}
    for index, method in enumerate(settings["methods"]):
        def report(fraction, index=index, method=method):
            progress.progress((index + fraction) / len(settings["methods"]), text=f"{METHODS[method]}…")
        results[method] = run_method(method, problem, solver_settings, ideal, nadir, starts, report)
    progress.empty()
    previous = st.session_state.sc_results
    st.session_state.sc_results = {
        "version": 0 if previous is None else previous["version"] + 1,  # resets the table selections
        "settings": settings, "n": n, "letter": letter, "domain": domain, "payoff": table, "ideal": ideal,
        "nadir": nadir, "reference": reference, "results": results,
        "latex": [to_latex(e) for e in expressions],
    }

stored = st.session_state.sc_results
if stored is None:
    st.info("Ustaw problem i metody w panelu bocznym, a potem kliknij „Oblicz”.")
    st.stop()
if stored["settings"] != settings:
    st.warning("Ustawienia zmieniły się od ostatniego obliczenia – poniżej wyniki dla poprzednich ustawień. "
               "Kliknij „Oblicz”, aby je przeliczyć.")

# ---------------------------------------------------------------- results
n_stored, letter_stored, results = stored["n"], stored["letter"], stored["results"]
st.subheader("Wyniki")
with st.expander("Tablica wypłat, punkt idealny i nadir"):
    criteria = [f"{letter_stored}{i + 1}" for i in range(n_stored)]
    payoff = pd.DataFrame(stored["payoff"], columns=criteria,
                          index=[f"min {name}" for name in criteria])
    payoff.loc["punkt idealny"] = stored["ideal"]
    payoff.loc["nadir"] = stored["nadir"]
    if "S3" in results:
        payoff.loc["punkt odniesienia r"] = stored["reference"]
    st.dataframe(payoff.style.format("{:.5g}"), width="stretch")
    st.caption("Wiersz „min Fᵢ” to oceny w punkcie minimalizującym i-te kryterium (remisy rozstrzygane "
               "sumą pozostałych). Nadir jest szacowany z tablicy wypłat – używa go S2 (zakres aᵢ) i S3 (punkt r).")

options_column, curve_column, distance_column = st.columns([2, 3, 2], vertical_alignment="center")
show_dominated = options_column.toggle("Pokaż punkty zdominowane", key="sc_show_dominated")
show_curve = curve_column.toggle("Spline / powierzchnia przez punkty niezdominowane", key="sc_show_curve")
ideal_p_label = distance_column.selectbox("Odległość od punktu idealnego: p", list(P_OPTIONS)[:3], key="sc_ideal_p",
                                          help="Metryka, w której szukane jest rozwiązanie najbliższe punktowi "
                                               "idealnemu (każde kryterium znormalizowane zakresem nadir − idealny).")
ideal_p = P_OPTIONS[ideal_p_label]
analyses = {method: analyze(result, stored, ideal_p) for method, result in results.items()}

tab_names = [METHODS[m] for m in results] + (["Porównanie metod"] if len(results) > 1 else [])
tabs = st.tabs(tab_names)
for tab, (method, result) in zip(tabs, results.items()):
    info = analyses[method]
    with tab:
        st.latex(method_latex(method, letter_stored, stored["settings"]))
        nondominated_ids = [i for i in info["solved_ids"] if info["nondominated"][i]]
        metrics = st.columns(5)
        metrics[0].metric("Zadania", len(result.runs))
        metrics[1].metric("Rozwiązane", len(result.solved))
        metrics[2].metric("Niezdominowane", len(nondominated_ids))
        metrics[3].metric("Unikalne niezdominowane", sum(info["unique"][i] for i in nondominated_ids),
                          help="Różne punkty – wiele wartości parametru może dać ten sam punkt "
                               "(np. nieaktywne ograniczenie w S2, niewypukły front w S1).")
        metrics[4].metric("Czas [s]", f"{result.time_s:.2f}")

        best_runs = []
        if info["best"] is not None:
            best = result.runs[info["best"]]
            best_runs = [best]
            st.success(f"**★ Najbliżej punktu idealnego** – zadanie nr {info['best']} ({run_label(best)}):  \n"
                       f"x = {vector(best.x)}, {letter_stored}(x) = {vector(best.F)}, "
                       f"d = {info['distance'][info['best']]:.4f}")
        st.latex(distance_latex(letter_stored, ideal_p_label))

        # The table is prepared first – its selection is highlighted on the plots
        table = result_table(result, info, n_stored, letter_stored, method, ideal_p_label)
        table_filter = st.session_state[f"sc_filter_{method}"]
        shown = filter_table(table, info, table_filter)
        table_key = f"sc_table_{method}_{table_filter}_{stored['version']}"
        selection = st.session_state.get(table_key)
        selected_ids = shown["nr"].iloc[selection["selection"]["rows"]].tolist() if selection else []
        selected_runs = [result.runs[i] for i in selected_ids if result.runs[i].F is not None]

        show_plots({method: result}, stored, show_dominated, show_curve, key=f"plot_{method}",
                   highlights=[("★ najbliżej punktu idealnego", best_runs, COMPROMISE_COLOR, 18),
                               ("Zaznaczone w tabeli", selected_runs, "magenta", 22)])

        st.radio("Pokaż", list(TABLE_FILTERS), format_func=TABLE_FILTERS.get, horizontal=True,
                 key=f"sc_filter_{method}")
        st.dataframe(shown, key=table_key, on_select="rerun", selection_mode="multi-row", hide_index=True,
                     width="stretch", height=320,
                     column_config=table_column_config(table, letter_stored, method, ideal_p_label))
        st.caption(f"Wierszy: {len(shown)} / {len(table)}. Zaznacz wiersze (pole po lewej), aby wyróżnić punkty "
                   "na wykresach; kliknięcie nagłówka sortuje.")
        st.download_button("Pobierz wyniki (CSV)", table.to_csv(index=False).encode("utf-8"),
                           file_name=f"skalaryzacja_{method}_{letter_stored}.csv", mime="text/csv",
                           key=f"download_{method}")

if len(results) > 1:
    with tabs[-1]:
        rows, best_overall = [], None
        for method, result in results.items():
            info = analyses[method]
            nondominated_ids = [i for i in info["solved_ids"] if info["nondominated"][i]]
            row = {"metoda": METHODS[method], "zadania": len(result.runs), "rozwiązane": len(result.solved),
                   "niezdominowane": len(nondominated_ids),
                   "unikalne": sum(info["unique"][i] for i in nondominated_ids)}
            if info["best"] is not None:
                best = result.runs[info["best"]]
                distance = float(info["distance"][info["best"]])
                row.update({"najbliższy – parametry": run_label(best), "x": vector(best.x),
                            f"{letter_stored}(x)": vector(best.F), f"d (p={ideal_p_label})": distance})
                if best_overall is None or distance < best_overall[1]:
                    best_overall = (method, distance)
            row["czas [s]"] = round(result.time_s, 3)
            rows.append(row)
        summary = pd.DataFrame(rows)
        if best_overall:
            summary.insert(0, "★", summary["metoda"] == METHODS[best_overall[0]])
        st.dataframe(summary, hide_index=True, width="stretch",
                     column_config={f"d (p={ideal_p_label})": st.column_config.NumberColumn(format="%.4f"),
                                    "★": st.column_config.CheckboxColumn(
                                        help="Metoda, która znalazła punkt najbliższy punktowi idealnemu")})
        st.latex(distance_latex(letter_stored, ideal_p_label))
        highlights = [(f"{'★' if best_overall and m == best_overall[0] else '○'} {m} – najbliżej punktu idealnego",
                       [] if analyses[m]["best"] is None else [results[m].runs[analyses[m]["best"]]],
                       COMPROMISE_COLOR if best_overall and m == best_overall[0] else METHOD_COLORS[m], 18)
                      for m in results]
        highlights.sort(key=lambda highlight: highlight[0].startswith("★"))  # the overall best is drawn on top
        show_plots(results, stored, False, show_curve, key="plot_all", highlights=highlights)

st.caption("Spline: w 2D krzywa PCHIP F₂(F₁) (monotoniczna, bez oscylacji) w przestrzeni ocen i krzywa sklejana "
           "3. stopnia w przestrzeni decyzji; w 3D powierzchnia z triangulacji siatki parametrów. Krzywa jest "
           "przerywana w miejscach dużych luk – np. tam, gdzie S1 nie znajduje punktów niewypukłej części frontu.")
