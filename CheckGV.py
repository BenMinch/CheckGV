#!/usr/bin/env python3
"""
gv_completeness.py
===================

A standalone giant-virus (Nucleocytoviricota / NCLDV) genome completeness and
contamination estimator, built to replicate the *core logic* of the
completeness/contamination module inside GVClass (Pitot, Brůna & Schulz 2024,
npj Viruses 2:60, https://doi.org/10.1038/s44298-024-00069-7), using only
public marker HMMs so it can run without the internal GVClass database.

WHAT IT DOES
------------
For each input genome (a protein FASTA, or a nucleotide FASTA that will be
gene-called with pyrodigal first):

  1. Searches the proteins against a panel of conserved, mostly single-copy
     NCLDV marker HMMs (GVOGs) with pyhmmer.
  2. Counts, for each panel (GVOG7 / GVOG8 / GVOG9 / the 35-marker "extended"
     panel), how many distinct markers are detected and how many total hits
     there are.
  3. From that:
        completeness (%) = distinct_markers_detected / panel_size * 100
        duplication factor = total_hits / distinct_markers_detected
     -- exactly the formula GVClass's paper describes for its order-level
     panel ("the total number of hits is divided by the unique hits"), just
     applied to a public, pan-order marker set instead of GVClass's internal
     per-order panel (see marker_panels.py for why, and for the caveat this
     implies).
  4. Bins completeness/duplication into low/medium/high using the paper's
     own thresholds (completeness: <30/30-70/>70; duplication: <1.5/1.5-2/>3).
  5. Optionally merges in a TIGTOG (Ha & Aylward 2024, npj Viruses 2:9)
     prediction table to report the predicted order/family alongside the
     quality metrics as context (TIGTOG is a fast random-forest classifier
     that assigns lineage from protein-family content; it is NOT used here
     to reweight the marker panel -- see the module docstring in
     marker_panels.py).

WHAT IT DELIBERATELY DOES NOT DO
---------------------------------
It does not reproduce GVClass's order-specific 576-marker panel or its
trained ExtraTrees contamination model (`extra_trees_v1`), because both ship
only inside GVClass's internal resource bundle. It also does not flag
eukaryote-like HGT genes as contamination (matching GVClass's own design
choice) -- it only ever measures duplication/absence of the core viral
single-copy markers themselves.

USAGE
-----
    # Protein input, single genome
    python CheckGV.py --faa genome.faa -o out.tsv

    # Nucleotide input (gene-called internally with pyrodigal)
    python CheckGV.py --fna genome.fna -o out.tsv

    # A whole directory of bins (protein or nucleotide, mixed OK)
    python CheckGV.py --input-dir bins/ -o summary.tsv

    # Attach TIGTOG lineage predictions as context columns
    python CheckGV.py --input-dir bins/ \\
        --tigtog Good_Bins_TIGTOG_prediction_result.tsv -o summary.tsv

Requires: pyhmmer, pyrodigal (both pure-Python-installable via pip; no
external HMMER or Prodigal binary needed).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import pyhmmer

from marker_panels import PANELS, CORE_GVOG_NAMES, COMPLETENESS_BINS, DUPLICATION_BINS, bin_value

def _as_str(x):
    """pyhmmer returns str or bytes for names depending on version; normalize."""
    return x.decode() if isinstance(x, bytes) else x


HMM_DIR = Path(__file__).parent / "hmm"
COMMON_HMM = HMM_DIR / "gvogs.common.hmm"

# pyhmmer / HMMER search thresholds. GVClass's own "sensitive mode" (its
# default) accepts hits at E <= 1e-5 for both the full-sequence and the
# best-domain E-value; we mirror that here rather than using the curated
# per-model gathering (GA) cutoffs, since GA cutoffs are only meaningful for
# the exact HMM build GVClass ships.
E_VALUE_CUTOFF = 1e-5

# On top of the E-value cutoff, GVClass's gene-calling step separately
# requires ">60% of model coverage" for a hit to count as a "complete
# profile hit". We apply the same idea here at the marker-counting stage:
# without it, a permissive E-value cutoff lets short, weak partial
# domain-only matches (common between distantly related helicases,
# ATPases, etc.) count as spurious extra "copies" of a marker and inflate
# the duplication factor. A hit only counts if its best domain covers at
# least this fraction of the marker's HMM length.
MIN_MODEL_COVERAGE = 0.4

FNA_EXTS = {".fna", ".fa", ".fasta"}
FAA_EXTS = {".faa", ".fasta.faa", ".pep"}


# --------------------------------------------------------------------------
# Gene calling (nucleotide input only)
# --------------------------------------------------------------------------
def call_genes(fna_path: Path, translation_tables=(11, 4, 1, 6, 15, 29)):
    """
    Gene-call a nucleotide FASTA with pyrodigal, trying meta mode plus a
    handful of single-genome-trained genetic codes, and keeping whichever
    call set yields the most conserved-marker hits downstream.

    This mirrors (in a simplified form) GVClass's own "opgecall" step, which
    tests multiple genetic codes and picks the winner by marker-hit count,
    average bit score, and coding density. We only use marker-hit count here
    to keep the tool self-contained; see the how-it-works notes in the
    module docstring.

    Returns a dict: {translation_table_label: [(seq_id, protein_seq), ...]}
    so the caller can hmmsearch each candidate and pick the best.
    """
    import pyrodigal

    seqs = list(_read_fasta(fna_path))
    if not seqs:
        raise ValueError(f"No sequences found in {fna_path}")

    candidates = {}

    # Meta mode (pretrained models, closest to GVClass's default "code 0").
    try:
        gf = pyrodigal.GeneFinder(meta=True)
        proteins = []
        for seq_id, seq in seqs:
            genes = gf.find_genes(seq.encode())
            for i, gene in enumerate(genes, 1):
                proteins.append((f"{seq_id}_{i}", gene.translate()))
        candidates["codemeta"] = proteins
    except Exception as exc:  # pragma: no cover
        print(f"  [warn] pyrodigal meta mode failed: {exc}", file=sys.stderr)

    # Single-genome training under specific genetic codes. Requires enough
    # sequence to train on; skip gracefully if too short.
    total_len = sum(len(s) for _, s in seqs)
    if total_len >= 20_000:
        for table in translation_tables:
            if table == 11:
                continue  # meta mode already covers the standard code well
            try:
                gf = pyrodigal.GeneFinder(meta=False)
                gf.train(*[s.encode() for _, s in seqs], translation_table=table)
                proteins = []
                for seq_id, seq in seqs:
                    genes = gf.find_genes(seq.encode())
                    for i, gene in enumerate(genes, 1):
                        proteins.append((f"{seq_id}_{i}", gene.translate()))
                candidates[f"code{table}"] = proteins
            except Exception as exc:
                print(f"  [warn] pyrodigal table {table} failed: {exc}", file=sys.stderr)

    return candidates, seqs


def _read_fasta(path: Path):
    seq_id, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if seq_id is not None:
            yield seq_id, "".join(chunks)


def gc_percent(seqs):
    seq = "".join(s for _, s in seqs).upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return round(100.0 * gc / len(seq), 2)


# --------------------------------------------------------------------------
# Marker search
# --------------------------------------------------------------------------
def load_hmms():
    with pyhmmer.plan7.HMMFile(str(COMMON_HMM)) as hf:
        hmms = list(hf)
    names = [_as_str(h.name).replace(".trim", "") for h in hmms]
    return hmms, names


def load_extended_panel(names):
    """The full common-GVOG panel, minus the 9-marker core (avoid double
    counting) -- used as PANELS['extended'] covering all 35 markers total
    (core + extra)."""
    PANELS["extended"] = list(names)


def search_proteins(proteins, hmms, alphabet):
    """
    proteins: list of (id, sequence) tuples (amino acids)
    hmms: list of pyhmmer.plan7.HMM
    Returns: dict marker_name -> list of (protein_id, bitscore, evalue) hits
    passing E_VALUE_CUTOFF, keeping only the best hit per protein per marker.
    """
    digital_seqs = []
    for i, (seq_id, seq) in enumerate(proteins):
        clean = re.sub(r"[^A-Za-z]", "", seq).upper().replace("*", "")
        if not clean:
            continue
        bio_seq = pyhmmer.easel.TextSequence(name=str(i).encode(), sequence=clean)
        digital_seqs.append(bio_seq.digitize(alphabet))

    id_lookup = [seq_id for seq_id, seq in proteins if re.sub(r"[^A-Za-z]", "", seq)]

    hits_by_marker = defaultdict(list)
    if not digital_seqs:
        return hits_by_marker

    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, digital_seqs)
    for hits in pyhmmer.hmmsearch(hmms, block, E=E_VALUE_CUTOFF, domE=E_VALUE_CUTOFF, cpus=0):
        marker_name = _as_str(hits.query.name).replace(".trim", "")
        model_length = hits.query.M
        for hit in hits:
            if hit.evalue > E_VALUE_CUTOFF:
                continue
            dom = hit.best_domain
            aln = dom.alignment
            coverage = (aln.hmm_to - aln.hmm_from + 1) / model_length if model_length else 0
            if coverage < MIN_MODEL_COVERAGE:
                continue
            protein_id = id_lookup[int(_as_str(hit.name))]
            hits_by_marker[marker_name].append((protein_id, hit.score, hit.evalue))
    return hits_by_marker


# --------------------------------------------------------------------------
# Completeness / duplication computation
# --------------------------------------------------------------------------
def compute_panel_stats(hits_by_marker, panel_markers):
    """
    hits_by_marker: marker_name -> list of (protein_id, score, evalue)
    panel_markers: list of marker names making up this panel

    Returns dict with distinct_present, panel_size, total_hits,
    completeness_pct, duplication_factor, markers_present (sorted list),
    markers_missing (sorted list).
    """
    present = [m for m in panel_markers if hits_by_marker.get(m)]
    total_hits = sum(len(hits_by_marker.get(m, [])) for m in panel_markers)
    n_present = len(present)
    panel_size = len(panel_markers)
    completeness_pct = round(100.0 * n_present / panel_size, 2) if panel_size else 0.0
    dup = round(total_hits / n_present, 2) if n_present else 0.0
    missing = [m for m in panel_markers if m not in present]
    return {
        "panel_size": panel_size,
        "distinct_present": n_present,
        "total_hits": total_hits,
        "completeness_pct": completeness_pct,
        "duplication_factor": dup,
        "markers_present": present,
        "markers_missing": missing,
    }


# --------------------------------------------------------------------------
# TIGTOG merge
# --------------------------------------------------------------------------
def load_tigtog(path: Path):
    table = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = row.get("Sequence") or row.get("sequence") or row.get("genome")
            if key:
                table[key] = row
    return table


def match_tigtog(genome_name, tigtog_table):
    if genome_name in tigtog_table:
        return tigtog_table[genome_name]
    # Fall back to a loose substring match, since TIGTOG's "Sequence" column
    # often carries a "-contigs" suffix that may not match the input filename
    # exactly.
    for key, row in tigtog_table.items():
        if key.startswith(genome_name) or genome_name.startswith(key):
            return row
    return None


# --------------------------------------------------------------------------
# Per-genome pipeline
# --------------------------------------------------------------------------
def process_genome(path: Path, hmms, hmm_names, alphabet, tigtog_table=None):
    is_nucleotide = path.suffix.lower() in FNA_EXTS
    genome_name = path.stem

    result = {"genome": genome_name, "input_file": path.name}

    if is_nucleotide:
        candidates, nt_seqs = call_genes(path)
        result["n_contigs"] = len(nt_seqs)
        result["total_length_bp"] = sum(len(s) for _, s in nt_seqs)
        result["gc_percent"] = gc_percent(nt_seqs)

        best_label, best_hits, best_proteins, best_score = None, None, None, -1
        for label, proteins in candidates.items():
            if not proteins:
                continue
            hits = search_proteins(proteins, hmms, alphabet)
            score = sum(len(v) for v in hits.values())
            if score > best_score:
                best_label, best_hits, best_proteins, best_score = label, hits, proteins, score
        if best_hits is None:
            best_label, best_hits, best_proteins = "codemeta", defaultdict(list), []

        result["ttable_selected"] = best_label
        result["n_proteins"] = len(best_proteins)
        hits_by_marker = best_hits
    else:
        proteins = list(_read_fasta(path))
        result["ttable_selected"] = "no_fna"
        result["n_contigs"] = "NA"
        result["total_length_bp"] = "NA"
        result["gc_percent"] = "NA"
        result["n_proteins"] = len(proteins)
        hits_by_marker = search_proteins(proteins, hmms, alphabet)

    # Per-panel stats
    for panel_name, panel_markers in PANELS.items():
        stats = compute_panel_stats(hits_by_marker, panel_markers)
        prefix = panel_name
        result[f"{prefix}_completeness_pct"] = stats["completeness_pct"]
        result[f"{prefix}_present"] = f"{stats['distinct_present']}/{stats['panel_size']}"
        result[f"{prefix}_duplication_factor"] = stats["duplication_factor"]
        if panel_name == "gvog9":
            result["gvog9_markers_present"] = ";".join(
                f"{m}({CORE_GVOG_NAMES.get(m, m)})" for m in stats["markers_present"]
            )
            result["gvog9_markers_missing"] = ";".join(
                f"{m}({CORE_GVOG_NAMES.get(m, m)})" for m in stats["markers_missing"]
            )

    # Headline numbers: the 35-marker "extended" panel is our best available
    # proxy for a genome-wide completeness/duplication estimate (broadest
    # coverage of conserved single-copy NCLDV genes short of GVClass's
    # internal order-specific panel -- see marker_panels.py docstring).
    result["estimated_completeness"] = result["extended_completeness_pct"]
    result["completeness_category"] = bin_value(result["estimated_completeness"], COMPLETENESS_BINS)
    result["estimated_duplication_factor"] = result["extended_duplication_factor"]
    result["contamination_category"] = bin_value(result["estimated_duplication_factor"], DUPLICATION_BINS)

    # TIGTOG context
    if tigtog_table:
        row = match_tigtog(genome_name, tigtog_table)
        if row:
            result["tigtog_predicted_order"] = row.get("Predicted_Order", "")
            result["tigtog_order_confidence"] = row.get("Confidence_Order_Pred", "")
            result["tigtog_predicted_family"] = row.get("Predicted_Family", "")
            result["tigtog_family_confidence"] = row.get("Confidence_Family_Pred", "")
        else:
            result["tigtog_predicted_order"] = ""
            result["tigtog_order_confidence"] = ""
            result["tigtog_predicted_family"] = ""
            result["tigtog_family_confidence"] = ""

    return result


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def gather_inputs(args):
    paths = []
    if args.faa:
        paths.append(Path(args.faa))
    if args.fna:
        paths.append(Path(args.fna))
    if args.input_dir:
        d = Path(args.input_dir)
        for p in sorted(d.iterdir()):
            if p.suffix.lower() in FNA_EXTS | {".faa"}:
                paths.append(p)
    return paths


def main():
    ap = argparse.ArgumentParser(
        description="Standalone NCLDV/giant-virus completeness & contamination estimator "
                    "(GVClass-style marker-panel approach, public HMMs only).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--faa", help="Single protein FASTA input")
    ap.add_argument("--fna", help="Single nucleotide FASTA input (gene-called internally)")
    ap.add_argument("--input-dir", help="Directory of .faa/.fna/.fa/.fasta genome files (one per genome/bin)")
    ap.add_argument("--tigtog", help="TIGTOG prediction TSV to merge in as lineage context columns")
    ap.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = ap.parse_args()

    inputs = gather_inputs(args)
    if not inputs:
        ap.error("Provide --faa, --fna, or --input-dir")

    print(f"Loading marker HMM panel from {COMMON_HMM} ...", file=sys.stderr)
    hmms, hmm_names = load_hmms()
    load_extended_panel(hmm_names)
    alphabet = pyhmmer.easel.Alphabet.amino()
    print(f"  {len(hmms)} marker models loaded "
          f"(core GVOG7/8/9 subset + extended {len(hmm_names)}-marker panel).", file=sys.stderr)

    tigtog_table = load_tigtog(Path(args.tigtog)) if args.tigtog else None
    if tigtog_table:
        print(f"Loaded TIGTOG predictions for {len(tigtog_table)} sequences.", file=sys.stderr)

    rows = []
    for path in inputs:
        print(f"Processing {path.name} ...", file=sys.stderr)
        try:
            rows.append(process_genome(path, hmms, hmm_names, alphabet, tigtog_table))
        except Exception as exc:
            print(f"  [error] failed on {path}: {exc}", file=sys.stderr)

    if not rows:
        print("No genomes processed successfully.", file=sys.stderr)
        sys.exit(1)

    fieldnames = list(rows[0].keys())
    for row in rows[1:]:
        for k in row:
            if k not in fieldnames:
                fieldnames.append(k)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"\nWrote {len(rows)} genome(s) to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
