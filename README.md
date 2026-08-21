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
