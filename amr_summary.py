#!/usr/bin/env python3
"""Summarise AMRFinderPlus results across many genomes.

Reads one AMRFinderPlus TSV file per genome and writes:
  - presence_absence.csv : genes x genomes matrix (1 = gene found)
  - gene_frequency.csv   : how many genomes carry each gene
  - genome_summary.csv   : hits, unique genes and resistance classes per genome
  - qc_report.txt        : files that were empty, malformed or filtered out

Only the Python standard library is needed.
"""

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path

LOG = logging.getLogger("amr_summary")

# AMRFinderPlus renamed some columns between versions, so accept both.
SYMBOL_COLS = ("Element symbol", "Gene symbol")
TYPE_COLS = ("Type", "Element type")
CLASS_COLS = ("Class",)
IDENTITY_COLS = ("% Identity to reference", "% Identity to reference sequence")
COVERAGE_COLS = ("% Coverage of reference", "% Coverage of reference sequence")


def first_present(header, candidates):
    """Return the first column name from candidates found in header, else None."""
    for name in candidates:
        if name in header:
            return name
    return None


def parse_amrfinder_file(path, min_identity=0.0, min_coverage=0.0, element_types=("AMR",)):
    """Parse one AMRFinderPlus TSV file.

    Returns (hits, problems). hits is a list of dicts with keys
    'gene', 'class' and 'type'. problems is a list of text messages
    describing anything unusual (an empty file, missing columns, bad rows).
    """
    hits, problems = [], []
    path = Path(path)

    if path.stat().st_size == 0:
        return hits, ["file is empty (0 bytes)"]

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []

        symbol_col = first_present(header, SYMBOL_COLS)
        if symbol_col is None:
            return hits, ["no gene symbol column found (expected 'Element symbol' or 'Gene symbol')"]

        type_col = first_present(header, TYPE_COLS)
        class_col = first_present(header, CLASS_COLS)
        identity_col = first_present(header, IDENTITY_COLS)
        coverage_col = first_present(header, COVERAGE_COLS)

        if min_identity > 0 and identity_col is None:
            problems.append("identity filter requested but no identity column found; filter skipped")
        if min_coverage > 0 and coverage_col is None:
            problems.append("coverage filter requested but no coverage column found; filter skipped")

        for line_number, row in enumerate(reader, start=2):
            gene = (row.get(symbol_col) or "").strip()
            if not gene:
                problems.append(f"line {line_number}: missing gene symbol, row skipped")
                continue

            hit_type = (row.get(type_col) or "").strip() if type_col else ""
            if element_types and type_col and hit_type not in element_types:
                continue

            try:
                if identity_col and min_identity > 0 and float(row[identity_col]) < min_identity:
                    continue
                if coverage_col and min_coverage > 0 and float(row[coverage_col]) < min_coverage:
                    continue
            except (TypeError, ValueError):
                problems.append(f"line {line_number}: could not read identity/coverage for {gene}, row skipped")
                continue

            hits.append({
                "gene": gene,
                "class": (row.get(class_col) or "").strip() if class_col else "",
                "type": hit_type,
            })

    return hits, problems


def summarise(results):
    """Combine per-genome results.

    results maps genome name -> (hits, problems). A genome counts as failed if
    its file could not be read at all (no hits and a problem was reported).
    Returns a dict with the tables needed for the output files.
    """
    presence = defaultdict(set)      # gene -> genomes carrying it
    gene_class = {}                  # gene -> resistance class
    genome_rows = []
    failed = []

    for genome in sorted(results):
        hits, problems = results[genome]
        unreadable = not hits and any(
            p.startswith(("file is empty", "no gene symbol column")) for p in problems
        )
        if unreadable:
            failed.append(genome)
            genome_rows.append({
                "genome": genome, "n_hits": 0, "n_unique_genes": 0,
                "classes": "", "status": "failed",
            })
            continue

        genes = {h["gene"] for h in hits}
        classes = sorted({h["class"] for h in hits if h["class"]})
        for hit in hits:
            presence[hit["gene"]].add(genome)
            gene_class.setdefault(hit["gene"], hit["class"])

        genome_rows.append({
            "genome": genome,
            "n_hits": len(hits),
            "n_unique_genes": len(genes),
            "classes": ";".join(classes),
            "status": "ok" if hits else "no_hits",
        })

    valid_genomes = [r["genome"] for r in genome_rows if r["status"] != "failed"]
    return {
        "presence": presence,
        "gene_class": gene_class,
        "genome_rows": genome_rows,
        "valid_genomes": valid_genomes,
        "failed": failed,
    }


def write_outputs(summary, results, out_dir):
    """Write the four output files into out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes = summary["valid_genomes"]
    genes = sorted(summary["presence"])

    with (out_dir / "presence_absence.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["gene", "class"] + genomes)
        for gene in genes:
            row = [gene, summary["gene_class"].get(gene, "")]
            row += [1 if g in summary["presence"][gene] else 0 for g in genomes]
            writer.writerow(row)

    with (out_dir / "gene_frequency.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["gene", "class", "n_genomes", "percent_genomes"])
        total = len(genomes)
        ranked = sorted(genes, key=lambda g: (-len(summary["presence"][g]), g))
        for gene in ranked:
            n = len(summary["presence"][gene])
            pct = round(100 * n / total, 1) if total else 0.0
            writer.writerow([gene, summary["gene_class"].get(gene, ""), n, pct])

    with (out_dir / "genome_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["genome", "n_hits", "n_unique_genes", "classes", "status"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary["genome_rows"])

    with (out_dir / "qc_report.txt").open("w") as handle:
        handle.write(f"Files read: {len(results)}\n")
        handle.write(f"Usable genomes: {len(genomes)}\n")
        handle.write(f"Failed files: {len(summary['failed'])}\n")
        no_hits = [r["genome"] for r in summary["genome_rows"] if r["status"] == "no_hits"]
        handle.write(f"Genomes with no AMR hits: {len(no_hits)}\n\n")
        for genome in sorted(results):
            for problem in results[genome][1]:
                handle.write(f"{genome}: {problem}\n")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Summarise AMRFinderPlus results across many genomes."
    )
    parser.add_argument("input_dir", help="folder containing one AMRFinderPlus TSV per genome")
    parser.add_argument("-o", "--output-dir", default="amr_summary_output",
                        help="folder for results (default: amr_summary_output)")
    parser.add_argument("-p", "--pattern", default="*.tsv",
                        help="file name pattern to look for (default: *.tsv)")
    parser.add_argument("--min-identity", type=float, default=0.0,
                        help="ignore hits below this %% identity (default: 0, no filter)")
    parser.add_argument("--min-coverage", type=float, default=0.0,
                        help="ignore hits below this %% coverage (default: 0, no filter)")
    parser.add_argument("--element-types", default="AMR",
                        help="comma-separated element types to keep, e.g. AMR,STRESS "
                             "(default: AMR; use '' to keep everything)")
    parser.add_argument("--strict", action="store_true",
                        help="exit with an error if any file failed to parse")
    parser.add_argument("-v", "--verbose", action="store_true", help="show detailed logging")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        LOG.error("input folder not found: %s", input_dir)
        return 2

    files = sorted(input_dir.glob(args.pattern))
    if not files:
        LOG.error("no files matching '%s' in %s", args.pattern, input_dir)
        return 2

    element_types = tuple(t.strip() for t in args.element_types.split(",") if t.strip())

    results = {}
    for path in files:
        genome = path.stem
        LOG.debug("reading %s", path)
        results[genome] = parse_amrfinder_file(
            path, args.min_identity, args.min_coverage, element_types
        )

    summary = summarise(results)
    write_outputs(summary, results, args.output_dir)

    LOG.info("%d files read, %d usable, %d failed",
             len(results), len(summary["valid_genomes"]), len(summary["failed"]))
    LOG.info("results written to %s", args.output_dir)

    if not summary["valid_genomes"]:
        LOG.error("no usable files; check the input folder and format")
        return 1
    if args.strict and summary["failed"]:
        LOG.error("failed files: %s", ", ".join(summary["failed"]))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
