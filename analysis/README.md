# SARS-CoV-2 phylogeography: Southeast Asia in a global context

A lean, self-contained pipeline that downsamples the global SARS-CoV-2 open
dataset hard, builds one time-resolved tree, reconstructs ancestral geography,
and renders phylogeographic maps with **transmission lines**.

It deliberately does **not** use the full `ncov` workflow, which is built for
routine global builds with GISAID credentials and a large compute budget. The
question here is narrower, so the pipeline is correspondingly smaller and
easier to audit.

## Demes

| Where | Deme |
|---|---|
| Southeast Asia | the **country** (see `defaults/sea_countries.txt`) |
| everywhere else | the **region** (Asia, Europe, Africa, Oceania, …) |
| Nextstrain's "North America" | split into North America / Central America / Caribbean |

That yields **17 demes**: 9 SE Asian countries with data, plus 8 world
regions. Note that "Asia" as a deme therefore means *non-Southeast* Asia.

The companion country-everywhere view needs no extra work: `country` is
reconstructed by the same `augur traits` run, so one tree serves both maps.

## Downsampling

SARS-CoV-2 open data is ~9.4M genomes, which no ML tree or mugration model can
take, and the sampling is wildly uneven between countries. `scripts/subsample.py`:

1. drops non-human hosts, failed Nextclade QC, low coverage, and ambiguous dates;
2. keeps at most **one genome per (Pango lineage, country, month)** — this
   flattens temporal *and* per-country sequencing effort at once;
3. caps any single country (`--max-per-country`), so Thailand and the USA do
   not crowd the map;
4. reserves a share of the budget for the focus countries, and fills the rest
   round-robin across countries.

```
9,446,229 rows  ->  70,974 after 1-per-cell  ->  4,000 genomes
                    773 lineages | 130 countries | 74 months
```

### A hard limit on the Singapore focus

**Singapore has only 32 SARS-CoV-2 genomes in the entire Nextstrain open
dataset** (against Thailand's 6,478 and Vietnam's 3,304), because Singapore
submitted to GISAID rather than to the INSDC/GenBank-derived open set. All 32
are force-included (`--force-include-countries`, which also bypasses the
per-cell cap and the QC filters for that country), but 32 tips is a thin basis
for estimating import/export, and the Singapore numbers on these figures
should be read as indicative only. Running this against GISAID data would be
the fix, and needs credentials.

## Figures

`scripts/plot-deme-map.py` is shared with the dengue analysis. Pie sizes and
line widths are **binned**, not continuous, because Auspice's continuous
scaling is unreadable when one deme has 900 tips and another has 2.

Three options exist specifically to keep a dense SARS-CoV-2 map legible:

| Option | Effect |
|---|---|
| `--max-lines N` | draw only the N strongest transitions (focus lines are never dropped) |
| `--focus-only` | draw *only* transitions touching the focus deme |
| `--min-slice-fraction F` | merge pie slices below F into one grey "other" wedge |

```bash
snakemake --configfile config.yaml -c 32   # build
bash scripts/make-figures.sh Singapore     # render
```

Outputs, each as a 330 dpi PNG and a vector PDF:

| Figure | What |
|---|---|
| `ncov_deme_global` | global, SE Asian countries + world regions |
| `ncov_country_global` | global, country everywhere |
| `ncov_deme_sea` | cropped to South & Southeast Asia |
| `ncov_singapore_global` | only the transitions touching Singapore |
| `ncov_singapore_sea` | the same, cropped |

## Caveats

- `augur traits` is a discrete mugration model. The lines are inferred
  ancestral-state changes along branches, not migration-rate estimates with
  support values. For publication-grade claims use BEAST/DTA.
- Subsampling reduces but cannot remove sampling bias; a deme's apparent
  importance still partly reflects how much it sequenced and deposited openly.
- Pie positions on crowded maps may be displaced to stop them overlapping; a
  grey leader line joins a displaced pie to its true centroid.
