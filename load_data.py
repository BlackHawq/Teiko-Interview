"""Load the supplied cell-count CSV into a normalized SQLite database."""

import csv
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "inputs" / "cell-count.csv"
DB_PATH = ROOT / "teiko.db"
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")


def create_database() -> None:
    """Create a fresh database and load all CSV records."""
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {CSV_PATH}")

    with sqlite3.connect(DB_PATH) as connection:
        connection.executescript(
            """
            PRAGMA foreign_keys = ON;
            DROP VIEW IF EXISTS population_summary;
            DROP TABLE IF EXISTS cell_counts;
            DROP TABLE IF EXISTS samples;
            DROP TABLE IF EXISTS subjects;

            CREATE TABLE subjects (
                subject_id TEXT PRIMARY KEY,
                project TEXT NOT NULL,
                condition TEXT NOT NULL,
                age INTEGER,
                sex TEXT NOT NULL
            );

            CREATE TABLE samples (
                sample_id TEXT PRIMARY KEY,
                subject_id TEXT NOT NULL REFERENCES subjects(subject_id),
                sample_type TEXT NOT NULL,
                treatment TEXT NOT NULL,
                response TEXT,
                time_from_treatment_start INTEGER NOT NULL
            );

            CREATE TABLE cell_counts (
                sample_id TEXT NOT NULL REFERENCES samples(sample_id),
                population TEXT NOT NULL,
                count INTEGER NOT NULL CHECK (count >= 0),
                PRIMARY KEY (sample_id, population)
            );

            CREATE VIEW population_summary AS
            WITH totals AS (
                SELECT sample_id, SUM(count) AS total_count
                FROM cell_counts
                GROUP BY sample_id
            )
            SELECT
                s.sample_id AS sample,
                totals.total_count,
                c.population,
                c.count,
                (100.0 * c.count / totals.total_count) AS percentage,
                sb.project,
                sb.subject_id,
                sb.condition,
                sb.age,
                sb.sex,
                s.sample_type,
                s.treatment,
                s.response,
                s.time_from_treatment_start
            FROM cell_counts AS c
            JOIN samples AS s ON s.sample_id = c.sample_id
            JOIN subjects AS sb ON sb.subject_id = s.subject_id
            JOIN totals ON totals.sample_id = c.sample_id;

            CREATE INDEX idx_samples_filters
                ON samples(sample_type, treatment, response, time_from_treatment_start);
            CREATE INDEX idx_subjects_condition_sex
                ON subjects(condition, sex);
            """
        )

        with CSV_PATH.open(newline="", encoding="utf-8") as csv_file:
            rows = csv.DictReader(csv_file)
            for row in rows:
                connection.execute(
                    """
                    INSERT INTO subjects(subject_id, project, condition, age, sex)
                    VALUES (?, ?, ?, ?, ?)
                    ON CONFLICT(subject_id) DO NOTHING
                    """,
                    (row["subject"], row["project"], row["condition"], row["age"], row["sex"]),
                )
                connection.execute(
                    """
                    INSERT INTO samples(
                        sample_id, subject_id, sample_type, treatment, response,
                        time_from_treatment_start
                    ) VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        row["sample"],
                        row["subject"],
                        row["sample_type"],
                        row["treatment"],
                        row["response"] or None,
                        row["time_from_treatment_start"],
                    ),
                )
                connection.executemany(
                    """
                    INSERT INTO cell_counts(sample_id, population, count)
                    VALUES (?, ?, ?)
                    """,
                    [(row["sample"], population, row[population]) for population in POPULATIONS],
                )

    print(f"Loaded {DB_PATH}")


if __name__ == "__main__":
    create_database()