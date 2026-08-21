# Teiko-Interview

## Run the analysis

From the repository root:

```text
python load_data.py
```

This creates `teiko.db` with normalized `subjects`, `samples`, and `cell_counts`
tables and a `population_summary` view containing the Part 2 frequency table.

Install dashboard dependencies with `pip install -r requirements.txt`, then run:

```text
streamlit run dashboard.py
```

The dashboard provides the full sample summary, melanoma PBMC miraclib responder
versus non-responder boxplots, Mann-Whitney U test results, the requested baseline
subset, and the average baseline B-cell count for responding melanoma males across
all treatments.

## Teiko Interview Analysis

This project loads the supplied immune-cell CSV into SQLite, calculates
sample-level population frequencies, compares miraclib responders with
non-responders, and provides baseline subset analysis through a Streamlit
dashboard.

## Automated grading commands

Run these commands from the repository root. Python 3.10 or newer is required.

```text
make setup
make pipeline
make dashboard
```

`make setup` installs the dependencies in `requirements.txt`.

`make pipeline` runs the complete non-interactive workflow. It rebuilds
`teiko.db` from `inputs/cell-count.csv` and writes these artifacts to `outputs/`:

- `population_summary.csv`: Part 2 frequency table.
- `miraclib_response_comparison.csv`: Part 3 filtered sample-level data.
- `miraclib_response_statistics.csv`: Mann-Whitney U statistics and p-values.
- `miraclib_response_boxplot.html`: Part 3 responder/non-responder boxplot.
- `all_samples_composition_pie.html`: Aggregate cell-composition pie chart.
- `baseline_melanoma_pbmc_miraclib.csv`: Part 4 baseline subset.
- `baseline_cohort_counts.csv`: Part 4 project, response, and sex counts.
- `baseline_male_responder_b_cell_average.csv`: Part 4 average and calculation totals.
- `baseline_male_responder_b_cell_counts.csv`: Individual B-cell values used in that average.

`make dashboard` starts the interactive dashboard at the local Streamlit URL.
The dashboard supports CSV upload, sample selection, color-coded frequency
tables, composition pie charts, response boxplots, statistical results, and
the Part 4 cohort counts and calculation.

## Direct commands

The equivalent commands without Make are:

```text
python3 load_data.py
python3 pipeline.py
streamlit run dashboard.py
```

On Windows, use `py -3.12` in place of `python3` if needed.
