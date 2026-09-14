## CheckGV
A completeness and redundancy estimator for giant viruses of the phylum Nucleocytoviricota. 
## Install
```bash
git clone https://github.com/BenMinch/CheckGV
```
```bash
pip install -r requirements.txt --break-system-packages
```

No external HMMER or Prodigal binary needed — `pyhmmer` and `pyrodigal` are
pure pip installs.

## Run

```bash
# One protein FASTA
python CheckGV.py --faa genome.faa -o out.tsv

# One nucleotide FASTA (genes called internally)
python CheckGV.py --fna genome.fna -o out.tsv

# A whole folder of bins (.faa/.fna/.fa/.fasta, mixed OK)
python CheckGV.py --input-dir bins/ -o summary.tsv

# Attach TIGTOG lineage predictions as context columns
python CheckGV.py --input-dir bins/ \
    --tigtog TIGTOG_prediction_result.tsv -o summary.tsv
```

## What it computes

For each genome, proteins are searched (via `pyhmmer`) against a public
35-marker panel of conserved Nucleocytoviricota Orthologous Groups (GVOGs),
sourced from Frank Aylward's `ncldv_markersearch` repository
(https://github.com/faylward/ncldv_markersearch).

For each marker panel (GVOG7 / GVOG8 / GVOG9 "core" — the same
supermatrix-marker sets named in the GVClass paper — and the full
35-marker "extended" panel):

```
completeness (%) = distinct markers detected / panel size × 100
duplication factor = total marker hits / distinct markers detected
```

Completeness and duplication are then binned into low/medium/high using thresholds:

| Metric | Low | Medium | High |
|---|---|---|---|
| Completeness | <30% | 30–70% | >70% |
| Duplication factor | <1.5 | 1.5–3 | >3 |

If a `--tigtog` prediction table is supplied, the matching row's predicted
order/family and confidence are reported alongside the quality metrics as
**context only** — see the caveat below.

## Important caveat — read before trusting the numbers


This tool instead uses a single **pan-order** marker panel as a
lineage-agnostic proxy. That means:

- It will *systematically underestimate* completeness for genomes in
  lineages that legitimately lack one or more of the 35 markers (some
  giant virus lineages have lost individual core genes — Medusavirus, for
  example, correctly scores low here because it genuinely lacks several of
  the panel's RNA polymerase/transcription markers).

## Output columns

- `genome`, `input_file`, `ttable_selected`, `n_contigs`, `total_length_bp`,
  `gc_percent`, `n_proteins`
- Per panel (`gvog7`, `gvog8`, `gvog9`, `extended`):
  `{panel}_completeness_pct`, `{panel}_present` (x/N), `{panel}_duplication_factor`
- `gvog9_markers_present` / `gvog9_markers_missing` — named list of the 9
  core markers (MCP, PolB, RNAPL, RNAPS, SFII, TFIIB, TopoII, A32, VLTF3)
- `estimated_completeness`, `completeness_category`
- `estimated_duplication_factor`, `contamination_category`
- (if `--tigtog` given) `tigtog_predicted_order`, `tigtog_order_confidence`,
  `tigtog_predicted_family`, `tigtog_family_confidence`
