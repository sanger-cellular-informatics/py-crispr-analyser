# Python CRISPR Analyser

The Python CRISPR Analyser is a set of libraries and scripts that allows you to index [CRISPR](https://en.wikipedia.org/wiki/CRISPR)s from FASTA genome files and use the resulting index to find CRISPRS by using guide RNA sequences and to find off-targets for a given guide RNA sequence.

## About

This tool is a re-write of the original [CRISPR Analyser](https://github.com/htgt/CRISPR-Analyser) written in C++ with the intention that the Python version will be easier to use and extend.
The following workflow is supported:
- **Gathering** CRISPRs from a FASTA file specifying the PAM sequence and the length of the guide RNA,
- **Indexing** the CRISPRs into a binary formatted index for fast search,
- **Searching** for CRISPRs by using a guide RNA sequence,
- **Aligning** the guide RNA sequence to the CRISPRs to find off-targets.

## Installation

To install the Python CRISPR Analyser, first clone the repository:

```bash
git clone git@gitlab.internal.sanger.ac.uk:sci/py-crispr-analyser.git
```

Make sure you have Python 3.10 or later and install [Poetry](https://python-poetry.org/):

```bash
pip install poetry

```

Then install the dependencies:

```bash
poetry install
```

## Usage

To run the **Gather** command use the following command:

```bash
poetry run gather -i <input_fasta> -o <output_file> -p <pam_sequence> -l <guide_length>
```

- The Input File needs to be a FASTA file containing the genenome sequence. For example GRCh38 which can be downloaded from [Ensembl](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/),
- The PAM sequence which can consist of A, C, G, T and N (for any),
- The Guide Length in base pairs e.g. 20,
- The Output File which will be a CSV file (without headers) with columns:
  - Chromasome Name e.g. '18'
  - Position Start - the number of base pairs the start position is offset from the start of the Chromasome (using 5' to 3' orientation),
  - CRISPR Sequence - which includes the PAM,
  - PAM Right? - a 1 or 0 (true or false) indicating if the PAM is on the right side of the CRISPR,
  - Species ID - currently always set to 1.

## Testing

To run unit tests:

```bash
poetry run pytest
```
