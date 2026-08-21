"""Interactive Streamlit dashboard for the Teiko cell-count analysis."""

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from scipy.stats import mannwhitneyu


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "teiko.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

st.set_page_config(page_title="Teiko immune cell analysis", layout="wide")
st.title("Immune Cell Population Analysis")
st.caption("Clinical trial cell-count explorer")


@st.cache_data
def load_summary() -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as connection:
        return pd.read_sql_query("SELECT * FROM population_summary", connection)


@st.cache_data
def load_baseline_subset() -> pd.DataFrame:
    query = """
        SELECT sample, subject_id, condition, sex, treatment, response,
               time_from_treatment_start, sample_type
        FROM population_summary
        WHERE condition = 'melanoma'
          AND sample_type = 'PBMC'
          AND treatment = 'miraclib'
          AND time_from_treatment_start = 0
        GROUP BY sample, subject_id, condition, sex, treatment, response,
                 time_from_treatment_start, sample_type
        ORDER BY sample
    """
    with sqlite3.connect(DB_PATH) as connection:
        return pd.read_sql_query(query, connection)


@st.cache_data
def baseline_male_responder_average() -> float | None:
    query = """
        SELECT AVG(c.count) AS average_b_cells
        FROM cell_counts AS c
        JOIN samples AS s ON s.sample_id = c.sample_id
        JOIN subjects AS subject ON subject.subject_id = s.subject_id
        WHERE c.population = 'b_cell'
          AND subject.condition = 'melanoma'
          AND subject.sex = 'M'
          AND s.response = 'yes'
          AND s.time_from_treatment_start = 0
    """
    with sqlite3.connect(DB_PATH) as connection:
        value = connection.execute(query).fetchone()[0]
    return value


def comparison_data(summary: pd.DataFrame) -> pd.DataFrame:
    filtered = summary[
        (summary["condition"] == "melanoma")
        & (summary["sample_type"] == "PBMC")
        & (summary["treatment"] == "miraclib")
        & (summary["response"].isin(["yes", "no"]))
    ].copy()
    filtered["group"] = filtered["response"].map(
        {"yes": "Responder", "no": "Non-responder"}
    )
    return filtered


def statistical_results(comparison: pd.DataFrame) -> pd.DataFrame:
    results = []
    for population in POPULATIONS:
        responder = comparison.loc[
            comparison["population"] == population, "percentage"
        ]
        responder = responder[comparison.loc[responder.index, "response"] == "yes"]
        non_responder = comparison.loc[
            comparison["population"] == population, "percentage"
        ]
        non_responder = non_responder[
            comparison.loc[non_responder.index, "response"] == "no"
        ]
        if responder.empty or non_responder.empty:
            statistic, p_value = float("nan"), float("nan")
        else:
            statistic, p_value = mannwhitneyu(
                responder, non_responder, alternative="two-sided"
            )
        results.append(
            {
                "population": population,
                "responder_n": len(responder),
                "non_responder_n": len(non_responder),
                "responder_median_pct": responder.median(),
                "non_responder_median_pct": non_responder.median(),
                "U_statistic": statistic,
                "p_value": p_value,
                "significant_at_0.05": p_value < 0.05,
            }
        )
    return pd.DataFrame(results)


if not DB_PATH.exists():
    st.error("teiko.db is missing. Run `python load_data.py` from the repository root.")
    st.stop()

summary = load_summary()

st.subheader("Part 2: Relative frequency by sample")
st.dataframe(
    summary[["sample", "total_count", "population", "count", "percentage"]],
    use_container_width=True,
    hide_index=True,
)

st.subheader("Part 3: Miraclib response comparison")
comparison = comparison_data(summary)
if comparison.empty:
    st.warning("No matching melanoma PBMC miraclib samples were found.")
else:
    figure = px.box(
        comparison,
        x="population",
        y="percentage",
        color="group",
        points="all",
        category_orders={"population": POPULATIONS},
        labels={"percentage": "Relative frequency (%)", "population": "Population"},
        color_discrete_map={"Responder": "#0f766e", "Non-responder": "#d97706"},
    )
    figure.update_layout(legend_title_text="Response")
    st.plotly_chart(figure, use_container_width=True)

    st.markdown("**Mann-Whitney U tests**")
    st.caption("Two-sided tests compare sample-level relative frequencies; alpha = 0.05.")
    results = statistical_results(comparison)
    st.dataframe(
        results.style.format(
            {
                "responder_median_pct": "{:.2f}",
                "non_responder_median_pct": "{:.2f}",
                "U_statistic": "{:.2f}",
                "p_value": "{:.4g}",
            }
        ),
        use_container_width=True,
        hide_index=True,
    )

st.subheader("Part 4: Baseline melanoma PBMC miraclib samples")
baseline = load_baseline_subset()
st.write(f"Matching samples: {len(baseline):,}")
st.dataframe(baseline, use_container_width=True, hide_index=True)

average = baseline_male_responder_average()
st.metric(
    "Average B cells: melanoma males, responders, time = 0",
    "No matching samples" if average is None else f"{average:.2f}",
)
st.caption("The average uses all sample and treatment types, as requested, and only filters the listed patient attributes.")
