#!/usr/bin/env python3
"""
Heavily downsample SARS-CoV-2 metadata by lineage x country x month.

SARS-CoV-2 open data is ~9M genomes, which is far beyond what a
maximum-likelihood tree and a mugration model can handle, and the sampling is
wildly uneven between countries. This takes at most `--per-cell` genomes from
each (lineage, country, year-month) cell, which flattens both the temporal and
the per-country sequencing effort at once, then trims to a total budget.

Because the analysis focuses on Singapore and Southeast Asia, a fixed share of
the budget (`--focus-fraction`) is reserved for the focus countries, and the
remaining global backbone is trimmed round-robin across countries so that no
single heavily-sequenced country (USA, UK) crowds out the rest.

Reads a TSV on stdin (`--metadata -`) so it can be fed straight from
`zstd -dc`, avoiding ever writing the decompressed file to disk.
"""
import argparse
import sys
from collections import Counter, defaultdict

# Values meaning "unknown"
MISSING = {"", "?", "na", "n/a", "unknown", "nan", "none", "-", "missing",
           "not collected", "not applicable", "not provided"}


def is_missing(value):
    return value.strip().lower() in MISSING


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--metadata", required=True,
                        help="metadata TSV, or '-' for stdin")
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-strains", required=True)
    parser.add_argument("--lineage-column", default="Nextclade_pango")
    parser.add_argument("--fallback-lineage-column", default="pango_lineage")
    parser.add_argument("--per-cell", type=int, default=1,
                        help="max genomes per (lineage, country, month)")
    parser.add_argument("--max-sequences", type=int, default=4000)
    parser.add_argument("--focus-countries",
                        help="file of countries to reserve budget for")
    parser.add_argument("--focus-fraction", type=float, default=0.5)
    parser.add_argument(
        "--force-include-countries",
        help="file of countries to take EVERY usable genome from, bypassing "
             "the per-cell cap and the QC/coverage filters. Intended for a "
             "focus country that is barely represented in the data.")
    parser.add_argument(
        "--max-per-country", type=int, default=0,
        help="cap any single country's contribution (0 = no cap). Stops one "
             "heavily-sequenced country dominating the map.")
    parser.add_argument("--min-date", default="2019-12",
                        help="earliest YYYY-MM to keep")
    parser.add_argument("--max-date", default="9999-12")
    parser.add_argument("--min-coverage", type=float, default=0.93)
    parser.add_argument("--qc-status", default="good",
                        help="comma-separated QC_overall_status values to keep")
    parser.add_argument("--host", default="Homo sapiens")
    args = parser.parse_args()

    focus = set()
    if args.focus_countries:
        with open(args.focus_countries, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    focus.add(line)

    forced = set()
    if args.force_include_countries:
        with open(args.force_include_countries, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    forced.add(line)

    keep_qc = {s.strip().lower() for s in args.qc_status.split(",") if s.strip()}

    stream = sys.stdin if args.metadata == "-" else open(args.metadata, encoding="utf-8")
    header = stream.readline().rstrip("\n").rstrip("\r").split("\t")
    idx = {name: i for i, name in enumerate(header)}

    def get(fields, name):
        i = idx.get(name)
        return fields[i].strip() if i is not None and i < len(fields) else ""

    required = ["strain", "date", "region", "country"]
    for col in required:
        if col not in idx:
            sys.exit(f"metadata has no {col!r} column")

    cell_counts = Counter()
    selected = []          # (strain, country, lineage, month, raw_line)
    stats = Counter()
    total = 0

    for line in stream:
        line = line.rstrip("\n").rstrip("\r")
        if not line:
            continue
        total += 1
        fields = line.split("\t")

        date = get(fields, "date")
        # need an unambiguous year-month
        if len(date) < 7 or "X" in date[:7] or "x" in date[:7]:
            stats["bad date"] += 1
            continue
        month = date[:7]
        if month < args.min_date or month > args.max_date:
            stats["out of date range"] += 1
            continue

        region, country = get(fields, "region"), get(fields, "country")
        if is_missing(region) or is_missing(country):
            stats["missing geography"] += 1
            continue

        is_forced = country in forced

        host = get(fields, "host")
        if not is_forced and args.host and host and host != args.host:
            stats["non-human host"] += 1
            continue

        qc = get(fields, "QC_overall_status").lower()
        if not is_forced and keep_qc and qc and qc not in keep_qc:
            stats["failed QC"] += 1
            continue

        coverage = get(fields, "coverage")
        if coverage and not is_forced:
            try:
                if float(coverage) < args.min_coverage:
                    stats["low coverage"] += 1
                    continue
            except ValueError:
                pass

        lineage = get(fields, args.lineage_column)
        if is_missing(lineage):
            lineage = get(fields, args.fallback_lineage_column)
        if is_missing(lineage):
            stats["missing lineage"] += 1
            continue

        cell = (lineage, country, month)
        if not is_forced and cell_counts[cell] >= args.per_cell:
            stats["cell already full"] += 1
            continue
        cell_counts[cell] += 1
        selected.append((get(fields, "strain"), country, lineage, month, line))

    if stream is not sys.stdin:
        stream.close()

    print(f"subsample: read {total} rows", file=sys.stderr)
    for reason, n in stats.most_common():
        print(f"  dropped {n:>9}  {reason}", file=sys.stderr)
    print(f"subsample: {len(selected)} after 1-per-"
          f"(lineage,country,month) [{len(cell_counts)} cells]", file=sys.stderr)

    # --- trim to budget, protecting the focus countries -------------------
    def round_robin(rows, budget):
        """Take `budget` rows, cycling across countries for even coverage."""
        if len(rows) <= budget:
            return rows
        by_country = defaultdict(list)
        for row in rows:
            by_country[row[1]].append(row)
        # within a country, spread across months
        for country in by_country:
            by_country[country].sort(key=lambda r: (r[3], r[2], r[0]))
        order = sorted(by_country, key=lambda c: (-len(by_country[c]), c))
        out, exhausted, i = [], 0, 0
        while len(out) < budget and exhausted < len(order):
            exhausted = 0
            for country in order:
                bucket = by_country[country]
                if i < len(bucket):
                    out.append(bucket[i])
                    if len(out) >= budget:
                        break
                else:
                    exhausted += 1
            i += 1
        return out

    if args.max_per_country:
        capped, per_country = [], Counter()
        # spread each country's quota across months rather than taking the
        # earliest ones
        for row in sorted(selected, key=lambda r: (r[1], r[3], r[2], r[0])):
            if row[1] in forced or per_country[row[1]] < args.max_per_country:
                capped.append(row)
                per_country[row[1]] += 1
        dropped = len(selected) - len(capped)
        selected = capped
        print(f"subsample: per-country cap {args.max_per_country} dropped "
              f"{dropped}, leaving {len(selected)}", file=sys.stderr)

    focus_rows = [r for r in selected if r[1] in focus]
    other_rows = [r for r in selected if r[1] not in focus]

    focus_budget = min(len(focus_rows), int(args.max_sequences * args.focus_fraction))
    other_budget = max(0, args.max_sequences - focus_budget)
    # if the backbone cannot fill its budget, give the slack back to the focus
    if len(other_rows) < other_budget:
        focus_budget = min(len(focus_rows),
                           args.max_sequences - len(other_rows))
        other_budget = len(other_rows)

    kept = round_robin(focus_rows, focus_budget) + round_robin(other_rows, other_budget)
    kept.sort(key=lambda r: (r[3], r[1], r[0]))

    with open(args.output_metadata, "w", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in kept:
            fh.write(row[4] + "\n")
    with open(args.output_strains, "w", encoding="utf-8") as fh:
        for row in kept:
            fh.write(row[0] + "\n")

    by_country = Counter(r[1] for r in kept)
    print(f"\nsubsample: kept {len(kept)} genomes "
          f"({len([r for r in kept if r[1] in focus])} from focus countries)",
          file=sys.stderr)
    print(f"  distinct lineages {len({r[2] for r in kept})} | "
          f"countries {len(by_country)} | months {len({r[3] for r in kept})}",
          file=sys.stderr)
    print("  top countries: " + ", ".join(
        f"{c}={n}" for c, n in by_country.most_common(12)), file=sys.stderr)


if __name__ == "__main__":
    main()
