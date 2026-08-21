"""Run the complete non-interactive analysis pipeline."""

from pathlib import Path
import sqlite3

import pandas as pd
import plotly.express as px
from scipy.stats import mannwhitneyu

from load_data import DB_PATH, POPULATIONS, create_database


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "outputs"


def run_pipeline() -> None:
    create_database()
    with sqlite3.connect(DB_PATH) as connection:
        summary = pd.read_sql_query(
            "SELECT * FROM population_summary",
            connection,
        )
    OUTPUT_DIR.mkdir(exist_ok=True)

    summary_columns = ["sample", "total_count", "population", "count", "percentage"]
    summary[summary_columns].to_csv(OUTPUT_DIR / "population_summary.csv", index=False)

    comparison = summary[
        (summary["condition"] == "melanoma")
        & (summary["sample_type"] == "PBMC")
        & (summary["treatment"] == "miraclib")
        & (summary["response"].isin(["yes", "no"]))
    ].copy()
    comparison["group"] = comparison["response"].map(
        {"yes": "Responder", "no": "Non-responder"}
    )
    comparison[summary_columns + ["group"]].to_csv(
        OUTPUT_DIR / "miraclib_response_comparison.csv", index=False
    )

    statistics = []
    for population in POPULATIONS:
        responder = comparison.loc[
            (comparison["population"] == population)
            & (comparison["response"] == "yes"),
            "percentage",
        ]
        non_responder = comparison.loc[
            (comparison["population"] == population)
            & (comparison["response"] == "no"),
            "percentage",
        ]
        if responder.empty or non_responder.empty:
            u_statistic, p_value = float("nan"), float("nan")
        else:
            u_statistic, p_value = mannwhitneyu(
                responder, non_responder, alternative="two-sided"
            )
        statistics.append(
            {
                "population": population,
                "responder_n": len(responder),
                "non_responder_n": len(non_responder),
                "responder_median_pct": responder.median(),
                "non_responder_median_pct": non_responder.median(),
                "U_statistic": u_statistic,
                "p_value": p_value,
                "significant_at_0.05": p_value < 0.05,
            }
        )
    pd.DataFrame(statistics).to_csv(
        OUTPUT_DIR / "miraclib_response_statistics.csv", index=False
    )

    boxplot = px.box(
        comparison,
        x="population",
        y="percentage",
        color="group",
        points="all",
        category_orders={"population": POPULATIONS},
        labels={"percentage": "Relative frequency (%)", "population": "Population"},
        color_discrete_map={"Responder": "#0f766e", "Non-responder": "#d97706"},
        title="Melanoma PBMC miraclib responders versus non-responders",
    )
    boxplot.write_html(OUTPUT_DIR / "miraclib_response_boxplot.html", include_plotlyjs="cdn")

    composition = (
        summary.groupby("population", as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )
    pie = px.pie(
        composition,
        names="population",
        values="count",
        hole=0.35,
        title="Aggregate cell composition across all samples",
    )
    pie.write_html(OUTPUT_DIR / "all_samples_composition_pie.html", include_plotlyjs="cdn")

    baseline = summary[
        (summary["condition"] == "melanoma")
        & (summary["sample_type"] == "PBMC")
        & (summary["treatment"] == "miraclib")
        & (summary["time_from_treatment_start"] == 0)
    ].drop_duplicates(
        [
            "sample", "project", "subject_id", "condition", "sex", "treatment",
            "response", "time_from_treatment_start", "sample_type",
        ]
    )
    baseline_columns = [
        "sample", "project", "subject_id", "condition", "sex", "treatment",
        "response", "time_from_treatment_start", "sample_type",
    ]
    baseline[baseline_columns].to_csv(
        OUTPUT_DIR / "baseline_melanoma_pbmc_miraclib.csv", index=False
    )

    project_counts = (
        baseline.groupby("project")["sample"].nunique()
        .rename("sample_count").reset_index().sort_values("project")
    )
    response_counts = (
        baseline.loc[baseline["response"].isin(["yes", "no"])]
        .drop_duplicates(["subject_id", "response"])
        .groupby("response")["subject_id"].nunique()
        .rename("subject_count").reset_index()
    )
    sex_counts = (
        baseline.drop_duplicates(["subject_id", "sex"])
        .groupby("sex")["subject_id"].nunique()
        .rename("subject_count").reset_index()
    )
    pd.concat(
        [
            project_counts.rename(columns={"project": "category", "sample_count": "count"})
            .assign(category_type="project"),
            response_counts.rename(columns={"response": "category", "subject_count": "count"})
            .assign(category_type="response"),
            sex_counts.rename(columns={"sex": "category", "subject_count": "count"})
            .assign(category_type="sex"),
        ],
        ignore_index=True,
    )[["category_type", "category", "count"]].to_csv(
        OUTPUT_DIR / "baseline_cohort_counts.csv", index=False
    )

    male_responder_b_cells = summary[
        (summary["population"] == "b_cell")
        & (summary["condition"] == "melanoma")
        & (summary["sex"] == "M")
        & (summary["response"] == "yes")
        & (summary["time_from_treatment_start"] == 0)
    ][["sample", "sample_type", "treatment", "count"]].sort_values("sample")
    total_b_cells = int(male_responder_b_cells["count"].sum())
    matching_samples = len(male_responder_b_cells)
    average_b_cells = (
        None if matching_samples == 0 else round(total_b_cells / matching_samples, 2)
    )
    pd.DataFrame(
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
                "average_b_cells": average_b_cells,
            }
        ]
    ).to_csv(OUTPUT_DIR / "baseline_male_responder_b_cell_average.csv", index=False)
    male_responder_b_cells.to_csv(
        OUTPUT_DIR / "baseline_male_responder_b_cell_counts.csv", index=False
    )

    print(f"Generated analysis outputs in {OUTPUT_DIR}")


if __name__ == "__main__":
    run_pipeline()
