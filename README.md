# cpc2-rs

This repository rebuilds the original CPC2 standalone predictor as a Rust CLI.
The implementation keeps the CPC2 feature extraction and libsvm model, but
removes the Python, Biopython, and external `svm-scale` / `svm-predict`
runtime dependency chain. The trained model and scaling data are embedded in
the binary at compile time.

## Build

```bash
cargo build --release
```

Binary output:

```bash
target/release/cpc2-rs
```

## Run

Basic prediction:

```bash
target/release/cpc2-rs -i data/example.fa -o example_output
```

Check reverse strand and include longest ORF start:

```bash
target/release/cpc2-rs -i data/example.fa -o example_output -r --ORF
```

Also output the putative peptide sequence:

```bash
target/release/cpc2-rs -i data/example.fa -o example_output -r --ORF --peptide
```

Notes:

- `-o/--output` accepts either a final `.txt` path or a basename.
- If the output value does not end with `.txt`, the program appends `.txt`.
- `--peptide` replaces the old `CPC2_output_peptide.py` behavior.

## CLI

```text
Usage:
  cpc2-rs -i <input.fasta> [-o <output>] [-r] [--ORF] [--peptide]

Options:
  -i, --input <FILE>    Input sequence in FASTA format
  -o, --output <FILE>   Output file or basename
  -r, --reverse         Also check the reverse strand
  --ORF                 Add the 1-based longest ORF start column
  --peptide             Add the putative peptide column
  -h, --help            Show help
```

## Output

Default output columns:

```text
#ID    transcript_length    peptide_length    Fickett_score    pI    ORF_integrity    coding_probability    label
```

Optional columns:

- `ORF_Start` when `--ORF` is enabled
- `putative_peptide` when `--peptide` is enabled

## Validation

The Rust implementation is covered by unit tests and example-data regression
tests for both forward-only and reverse-strand modes:

```bash
cargo test
```

## Repository Layout

- `src/` Rust implementation
- `data/` CPC2 model, scaling range, and example FASTA
- `bin/` original Python scripts kept as reference
- `libs/` original libsvm source archive kept as reference
