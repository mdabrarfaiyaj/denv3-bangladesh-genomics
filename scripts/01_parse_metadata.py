"""
01_parse_metadata.py
--------------------
Turn NCBI Datasets' data_report.jsonl into one clean metadata table.

Input : data/raw/denv3_all/ncbi_dataset/data/data_report.jsonl
Output: results/metadata/metadata.tsv
        one row per sequence: accession, country, division, date, date_precision,
        year, length, completeness, host, release_date, submitter_affiliation

Run from the project folder:
  python scripts/01_parse_metadata.py
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

IN = Path("data/raw/denv3_all/ncbi_dataset/data/data_report.jsonl")
OUT = Path("results/metadata/metadata.tsv")


def get(d, *keys):
    """Safely walk nested dicts: get(rec, 'isolate', 'collectionDate')."""
    for k in keys:
        if not isinstance(d, dict):
            return None
        d = d.get(k)
    return d


def clean_date(raw):
    """Return (YYYY-MM-DD or YYYY-MM or YYYY, precision) from NCBI date strings."""
    if not raw:
        return None, "missing"
    raw = raw.strip()
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})", raw)
    if m:
        return m.group(0), "day"
    m = re.match(r"^(\d{4})-(\d{2})$", raw)
    if m:
        return m.group(0), "month"
    m = re.match(r"^(\d{4})$", raw)
    if m:
        return m.group(0), "year"
    m = re.search(r"(19|20)\d{2}", raw)          # e.g. "Jul-2019" or "2019/2020"
    if m:
        return m.group(0), "year"
    return None, "missing"


def main():
    if not IN.exists():
        sys.exit(f"Input not found: {IN} (download the data first, see README)")
    rows = []
    with IN.open() as fh:
        for line in fh:
            rec = json.loads(line)
            loc = get(rec, "location", "geographicLocation") or ""
            country, _, division = loc.partition(":")
            date, precision = clean_date(get(rec, "isolate", "collectionDate"))
            rows.append({
                "accession": rec.get("accession"),
                "country": country.strip() or None,
                "division": division.strip() or None,
                "date": date,
                "date_precision": precision,
                "year": int(date[:4]) if date else None,
                "length": rec.get("length"),
                "completeness": rec.get("completeness"),
                "host": get(rec, "host", "organismName"),
                "release_date": (rec.get("releaseDate") or "")[:10] or None,
                "submitter_affiliation": get(rec, "submitter", "affiliation"),
            })

    df = pd.DataFrame(rows)
    df["year"] = df["year"].astype("Int64")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)

    # ---- short report: check these numbers against the earlier awk counts ----
    print(f"Sequences: {len(df)}  ->  {OUT}")
    print(f"Missing country: {df['country'].isna().sum()} | missing date: {(df['date_precision'] == 'missing').sum()}")
    print("\nDate precision:\n" + df["date_precision"].value_counts().to_string())
    bd = df[df["country"] == "Bangladesh"]
    bd_kind = pd.cut(bd["length"], [0, 1399, 9999, 20000], labels=["short", "E-gene-size", "complete"])
    print(f"\nBangladesh: {len(bd)} sequences")
    print(pd.crosstab(bd["year"].fillna(0).astype(int), bd_kind).to_string())
    print("\nHosts (top 5):\n" + df["host"].value_counts().head(5).to_string())


if __name__ == "__main__":
    main()
