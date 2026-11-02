# amr-summary

A small command-line tool that summarises [AMRFinderPlus](https://github.com/ncbi/amr) results across many bacterial genomes. It turns a folder of per-genome result files into a gene presence/absence table, gene frequencies and a per-genome summary, and it reports any file that could not be read.

## Context

I wrote this while working on antimicrobial resistance in *Pseudomonas aeruginosa* genomes during my PhD at CSIR-IMTECH, where the same summary had to be rebuilt by hand for many strains. The tool is intentionally small and has no dependencies, so that colleagues and students can run and read it easily.

## What it does

For every AMRFinderPlus TSV file in a folder it:

1. reads the file (both the older `Gene symbol` and newer `Element symbol` column layouts are accepted);
2. keeps AMR hits only by default, and optionally drops hits below a minimum identity or coverage;
3. checks the file: empty files, missing columns and unreadable rows are reported instead of silently ignored;
4. combines all genomes into summary tables.

## Requirements

- Python 3.8 or newer
- No other packages

## Usage

```bash
python amr_summary.py INPUT_DIR -o OUTPUT_DIR
```

Options:

| Option | Meaning | Default |
| --- | --- | --- |
| `-o`, `--output-dir` | where results are written | `amr_summary_output` |
| `-p`, `--pattern` | file name pattern to read | `*.tsv` |
| `--min-identity` | ignore hits below this % identity | 0 (no filter) |
| `--min-coverage` | ignore hits below this % coverage | 0 (no filter) |
| `--element-types` | element types to keep, comma separated (`''` keeps all) | `AMR` |
| `--strict` | exit with an error if any file failed to parse | off |
| `-v`, `--verbose` | detailed logging | off |

The genome name is taken from the file name without its extension (`genome_A.tsv` becomes `genome_A`).

## Input

One AMRFinderPlus result table (tab-separated) per genome, for example produced with:

```bash
amrfinder -n genome_A.fasta -o genome_A.tsv
```

## Output

| File | Content |
| --- | --- |
| `presence_absence.csv` | genes (rows) by genomes (columns), 1 if the gene was found |
| `gene_frequency.csv` | number and percentage of genomes carrying each gene |
| `genome_summary.csv` | hits, unique genes, resistance classes and status (`ok`, `no_hits` or `failed`) per genome |
| `qc_report.txt` | counts of usable and failed files, plus each problem found |

A genome with a valid file but no AMR genes is reported as `no_hits` and stays in the analysis. A genome whose file could not be read is marked `failed` and left out of the matrix.

## Example

The folder `examples/amrfinder_results` holds four small **synthetic** files made for testing (they are not real data):

```bash
python amr_summary.py examples/amrfinder_results -o results
```

The expected output is stored in `examples/expected_output`, and one of the four example files is deliberately empty to show how the QC report works.

## Tests

```bash
python -m unittest discover -s tests -v
```

The tests cover both column layouts, empty and malformed files, the filters, and a full run on the example data compared with the expected output.

## Limitations

- It summarises gene presence only; it does not count copies or report point mutations.
- It reads AMRFinderPlus output only. Other tools (for example CARD RGI) need a different parser.
- Gene names are used as found, so differences in naming between database versions are not harmonised. Run all genomes with the same AMRFinderPlus version and database.

## Use of AI

I used an AI assistant to help draft the code structure and the tests. I reviewed every function, ran the tests and checked the output against the example data by hand.

## Licence

MIT, see `LICENSE`.
