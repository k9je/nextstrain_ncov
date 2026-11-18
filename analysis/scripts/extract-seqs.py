#!/usr/bin/env python3
"""
Pull a named subset of sequences out of a FASTA stream.

The full aligned SARS-CoV-2 open dataset is ~1.2 GB compressed and well over
100 GB decompressed, so it is never written to disk: this reads the stream on
stdin (fed by `zstd -dc`) and emits only the wanted records.

Strain names in the Nextstrain FASTA carry the same name as the metadata
`strain` column, but may have a `hCoV-19/` style prefix, so those are stripped
before matching.
"""
import argparse
import sys

PREFIXES = ("hCoV-19/", "SARS-CoV-2/", "BetaCoV/")


def normalise(name):
    name = name.strip().split()[0] if name.strip() else ""
    for prefix in PREFIXES:
        if name.startswith(prefix):
            name = name[len(prefix):]
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strains", required=True, help="one strain per line")
    parser.add_argument("--output", required=True)
    parser.add_argument("--input", default="-", help="FASTA, or '-' for stdin")
    args = parser.parse_args()

    with open(args.strains, encoding="utf-8") as fh:
        wanted = {normalise(line) for line in fh if line.strip()}
    print(f"extract-seqs: looking for {len(wanted)} strains", file=sys.stderr)

    stream = sys.stdin if args.input == "-" else open(args.input, encoding="utf-8")
    found = set()
    keeping = False
    n_records = 0

    with open(args.output, "w", encoding="utf-8") as out:
        for line in stream:
            if line.startswith(">"):
                n_records += 1
                name = normalise(line[1:])
                keeping = name in wanted and name not in found
                if keeping:
                    found.add(name)
                    # write the normalised name so it matches the metadata
                    out.write(f">{name}\n")
                if n_records % 1000000 == 0:
                    print(f"  scanned {n_records:,} records, "
                          f"found {len(found)}/{len(wanted)}", file=sys.stderr)
            elif keeping:
                out.write(line)

    if stream is not sys.stdin:
        stream.close()

    missing = wanted - found
    print(f"extract-seqs: scanned {n_records:,} records; "
          f"wrote {len(found)}/{len(wanted)}", file=sys.stderr)
    if missing:
        print(f"extract-seqs: {len(missing)} strains not found in the FASTA, "
              f"e.g. {sorted(missing)[:5]}", file=sys.stderr)


if __name__ == "__main__":
    main()
