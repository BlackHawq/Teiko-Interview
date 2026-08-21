"""Interactive Streamlit dashboard for the Teiko cell-count analysis."""

import hashlib
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
st.title("Immune Cell Population Analysis by Hammud Haq")
st.caption("Clinical trial cell-count explorer")


@st.cache_data
def load_summary() -> pd.DataFrame:
    with sqlite3.connect(DB_PATH) as connection:
        return pd.read_sql_query("SELECT * FROM population_summary", connection)


def summary_from_csv(uploaded_file) -> pd.DataFrame:
    raw = pd.read_csv(uploaded_file)
    required = {
        "project", "subject", "condition", "age", "sex", "treatment",
        "response", "sample", "sample_type", "time_from_treatment_start",
        *POPULATIONS,
    }
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    metadata = [
        "project", "subject", "condition", "age", "sex", "treatment",
        "response", "sample", "sample_type", "time_from_treatment_start",
    ]
    long = raw.melt(
        id_vars=metadata,
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long["count"] = pd.to_numeric(long["count"], errors="raise")
    long["time_from_treatment_start"] = pd.to_numeric(
        long["time_from_treatment_start"], errors="raise"
    )
    long["total_count"] = long.groupby("sample")["count"].transform("sum")
    long["percentage"] = 100 * long["count"] / long["total_count"]
    return long.rename(columns={"subject": "subject_id"})


def baseline_subset(summary: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "sample", "project", "subject_id", "condition", "sex", "treatment", "response",
        "time_from_treatment_start", "sample_type",
    ]
    return (
        summary.loc[
            (summary["condition"] == "melanoma")
            & (summary["sample_type"] == "PBMC")
            & (summary["treatment"] == "miraclib")
            & (summary["time_from_treatment_start"] == 0),
            columns,
        ]
        .drop_duplicates()
        .sort_values("sample")
    )


def baseline_male_responder_average(summary: pd.DataFrame) -> float | None:
    matching = male_responder_baseline_counts(summary)
    return None if matching.empty else float(matching["count"].mean())


def male_responder_baseline_counts(summary: pd.DataFrame) -> pd.DataFrame:
    return summary.loc[
        (summary["population"] == "b_cell")
        & (summary["condition"] == "melanoma")
        & (summary["sex"] == "M")
        & (summary["response"] == "yes")
        & (summary["time_from_treatment_start"] == 0),
        ["sample", "sample_type", "treatment", "count"],
    ].sort_values("sample")


def sample_color(sample: str) -> str:
    palette = px.colors.qualitative.Set2
    digest = hashlib.md5(sample.encode("utf-8")).hexdigest()
    return palette[int(digest[:8], 16) % len(palette)]


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


with st.sidebar:
    st.header("Dataset")
    uploaded_file = st.file_uploader(
        "Upload a cell-count CSV",
        type="csv",
        help="The CSV must contain the same metadata and five population columns as the supplied file.",
    )

if uploaded_file is not None:
    try:
        summary = summary_from_csv(uploaded_file)
        st.sidebar.success(f"Using uploaded data: {summary['sample'].nunique():,} samples")
    except (ValueError, pd.errors.ParserError) as error:
        st.error(f"Could not load the uploaded CSV: {error}")
        st.stop()
elif DB_PATH.exists():
    summary = load_summary()
else:
    st.error("teiko.db is missing. Run `python load_data.py` from the repository root.")
    st.stop()

st.sidebar.divider()
st.sidebar.header("Sample explorer")
sample_options = ["All samples"] + sorted(summary["sample"].unique().tolist())
selected_sample = st.sidebar.selectbox("Show composition for", sample_options)

st.subheader("Sample composition")
if selected_sample == "All samples":
    pie_data = (
        summary.groupby("population", as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )
    pie_title = "Aggregate cell composition across all samples"
else:
    pie_data = summary.loc[
        summary["sample"] == selected_sample, ["population", "count"]
    ]
    pie_title = f"Cell composition: {selected_sample}"
pie = px.pie(
    pie_data,
    names="population",
    values="count",
    hole=0.35,
    title=pie_title,
    color="population",
    color_discrete_sequence=px.colors.qualitative.Bold,
)
pie.update_traces(textposition="inside", textinfo="percent+label")
st.plotly_chart(pie, use_container_width=True)

st.subheader("Part 2: Relative frequency by sample")
displayed_summary = summary[["sample", "total_count", "population", "count", "percentage"]].copy()
displayed_summary["sample_color"] = displayed_summary["sample"].map(sample_color)
if selected_sample == "All samples":
    rows_per_page = 500
    page_count = (len(displayed_summary) + rows_per_page - 1) // rows_per_page
    page = min(st.session_state.get("summary_page", 1), page_count)
    start = (page - 1) * rows_per_page
    page_rows = displayed_summary.iloc[start : start + rows_per_page]
    st.dataframe(
        page_rows.style.apply(
            lambda row: [
                f"background-color: {row.sample_color}; color: #111827" if column == "sample" else ""
                for column in page_rows.columns
            ],
            axis=1,
        ),
        use_container_width=True,
        hide_index=True,
    )
    page = st.number_input(
        "Table page",
        min_value=1,
        max_value=page_count,
        value=page,
        step=1,
        key="summary_page",
        help=f"Showing {rows_per_page:,} color-coded rows per page.",
    )
    st.caption(
        f"Page {page:,} of {page_count:,}. Every sample entry is color coded; "
        "use the sample selector above to focus on one sample."
    )
else:
    selected_rows = displayed_summary[displayed_summary["sample"] == selected_sample]
    st.dataframe(
        selected_rows.style.apply(
            lambda row: [
                f"background-color: {row.sample_color}; color: #111827" if column == "sample" else ""
                for column in selected_rows.columns
            ],
            axis=1,
        ),
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
baseline = baseline_subset(summary)
st.write(f"Matching samples: {len(baseline):,}")
st.dataframe(baseline, use_container_width=True, hide_index=True)

st.markdown("**Baseline cohort counts**")
st.caption(
    "These counts describe the samples shown above: melanoma, PBMC, miraclib, "
    "and time_from_treatment_start = 0. Subjects are counted once per category."
)
project_counts = (
    baseline.groupby("project")["sample"]
    .nunique()
    .rename("sample_count")
    .reset_index()
    .sort_values("project")
)
response_counts = (
    baseline.loc[baseline["response"].isin(["yes", "no"])]
    .drop_duplicates(["subject_id", "response"])
    .groupby("response")["subject_id"]
    .nunique()
    .rename("subject_count")
    .reset_index()
    .replace({"response": {"yes": "Responder", "no": "Non-responder"}})
)
sex_counts = (
    baseline.drop_duplicates(["subject_id", "sex"])
    .groupby("sex")["subject_id"]
    .nunique()
    .rename("subject_count")
    .reset_index()
    .replace({"sex": {"M": "Male", "F": "Female"}})
)
count_columns = st.columns(3)
with count_columns[0]:
    st.markdown("*Samples by project*")
    st.dataframe(project_counts, use_container_width=True, hide_index=True)
with count_columns[1]:
    st.markdown("*Subjects by response*")
    st.dataframe(response_counts, use_container_width=True, hide_index=True)
with count_columns[2]:
    st.markdown("*Subjects by sex*")
    st.dataframe(sex_counts, use_container_width=True, hide_index=True)

male_responder_counts = male_responder_baseline_counts(summary)
average = None if male_responder_counts.empty else float(male_responder_counts["count"].mean())
st.metric(
    "Average B cells: melanoma males, responders, time = 0",
    "No matching samples" if average is None else f"{average:.2f}",
)
if male_responder_counts.empty:
    st.error("No melanoma male responder samples were found at time = 0.")
else:
    matching_samples = len(male_responder_counts)
    total_b_cells = int(male_responder_counts["count"].sum())
    st.success("Calculation complete: matching melanoma male responder samples found.")
    st.markdown(
        f"**Calculation:** {total_b_cells:,} total B cells / "
        f"{matching_samples:,} matching samples = "
        f"{total_b_cells:,} / {matching_samples:,} = **{average:.2f}**"
    )
    calculation_table = pd.DataFrame(
        [
            {
                "condition": "melanoma",
                "sex": "M",
                "response": "yes",
                "time_from_treatment_start": 0,
                "sample_types": "all",
                "treatments": "all",
                "matching_samples": matching_samples,
                "total_b_cells": total_b_cells,
                "average_b_cells": average,
            }
        ]
    )
    st.dataframe(
        calculation_table.style.format({"average_b_cells": "{:.2f}"}),
        use_container_width=True,
        hide_index=True,
    )
    with st.expander("Show individual B-cell counts used in the calculation"):
        st.dataframe(male_responder_counts, use_container_width=True, hide_index=True)
st.caption(
    "The average uses all sample and treatment types, as requested. "
    "Only melanoma, male, responder = yes, and time = 0 are filters."
)
