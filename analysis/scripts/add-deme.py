#!/usr/bin/env python3
"""
Add a `deme` column to a Nextstrain metadata TSV for phylogeographic analysis.

The deme is the geographic unit whose ancestral states `augur traits`
reconstructs, and hence the unit between which Auspice draws transmission
lines on the map. It is resolved in this order:

  1. country listed in --country-level-countries  -> deme = country
  2. region == --country-level-region             -> deme = country
  3. country listed in --overrides                -> deme = the override value
  4. otherwise                                    -> deme = region

so you can have country-level resolution for a named set of countries (e.g.
Southeast Asia) or for a whole region (e.g. all of Asia), with region-level
resolution everywhere else. The overrides file additionally lets a coarse
region be split, e.g. Nextstrain's "North America" into North America /
Central America / Caribbean.

Missing region, missing country where a country is required, and
non-geographic regions all yield deme="?" so `augur filter` drops them.

Input is read and written as plain tab-separated text with no quote handling,
which matches how Augur writes metadata and guarantees that every other column
is passed through byte-for-byte.
"""
import argparse
import sys
from collections import Counter

MISSING = "?"
MISSING_VALUES = {
    "", "?", "na", "n/a", "unknown", "nan", "none", "-",
    "missing", "not collected", "not applicable", "not provided",
}

# Nextstrain regions that are not a place anyone can be infected in.
NON_GEOGRAPHIC_REGIONS = {"oceans and seas"}


def is_missing(value):
    return value.strip().lower() in MISSING_VALUES


def is_non_geographic(region):
    return region.strip().lower() in NON_GEOGRAPHIC_REGIONS


def read_list(path):
    """Read a one-value-per-line file, ignoring blanks and # comments."""
    values = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r").strip()
            if line and not line.startswith("#"):
                values.append(line)
    return values


def read_overrides(path):
    """Read a `country<TAB>deme` TSV, ignoring blank lines and # comments."""
    overrides = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}: expected 2 tab-separated fields, got {len(fields)}: {line!r}"
                )
            country, deme = fields[0].strip(), fields[1].strip()
            if not country or not deme:
                raise ValueError(f"{path}: empty country or deme in line: {line!r}")
            overrides[country] = deme
    return overrides


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--metadata", required=True, help="input metadata TSV")
    parser.add_argument("--overrides", help="TSV of country<TAB>deme overrides")
    parser.add_argument(
        "--country-level-countries",
        help="file of countries that each become their own deme",
    )
    parser.add_argument(
        "--country-level-region", default=None,
        help="region whose countries each become their own deme (e.g. Asia)",
    )
    parser.add_argument(
        "--all-countries", action="store_true",
        help="every country is its own deme (the country-everywhere scheme)",
    )
    parser.add_argument("--region-column", default="region")
    parser.add_argument("--country-column", default="country")
    parser.add_argument("--deme-column", default="deme")
    parser.add_argument("--output", required=True, help="output metadata TSV")
    args = parser.parse_args()

    overrides = read_overrides(args.overrides) if args.overrides else {}
    country_level = set(
        read_list(args.country_level_countries)
        if args.country_level_countries else []
    )

    counts = Counter()
    n_rows = 0

    with open(args.metadata, encoding="utf-8") as fin, \
         open(args.output, "w", encoding="utf-8") as fout:

        header_line = fin.readline().rstrip("\n").rstrip("\r")
        if not header_line:
            raise ValueError(f"{args.metadata} is empty")
        header = header_line.split("\t")

        for required in (args.region_column, args.country_column):
            if required not in header:
                raise ValueError(
                    f"{args.metadata} has no {required!r} column; found: {header}"
                )
        if args.deme_column in header:
            raise ValueError(
                f"{args.metadata} already has a {args.deme_column!r} column"
            )

        region_idx = header.index(args.region_column)
        country_idx = header.index(args.country_column)
        fout.write("\t".join(header + [args.deme_column]) + "\n")

        for line in fin:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))

            region = fields[region_idx].strip()
            country = fields[country_idx].strip()

            if is_missing(region) or is_non_geographic(region):
                deme = MISSING
            elif args.all_countries or country in country_level or (
                    args.country_level_region
                    and region == args.country_level_region):
                deme = MISSING if is_missing(country) else country
            else:
                deme = overrides.get(country, region)

            counts[deme] += 1
            n_rows += 1
            fout.write("\t".join(fields + [deme]) + "\n")

    print(f"add-deme: wrote {n_rows} rows to {args.output}", file=sys.stderr)
    print(f"add-deme: {len(counts)} distinct deme values", file=sys.stderr)
    for deme, count in counts.most_common():
        label = "(missing - will be filtered out)" if deme == MISSING else ""
        print(f"  {count:7d}  {deme} {label}".rstrip(), file=sys.stderr)


if __name__ == "__main__":
    main()
