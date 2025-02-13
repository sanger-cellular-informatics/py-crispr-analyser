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

Python CRISPR Analyser can be used either as a library or run as a series of scritps.

### Gather

To run the **Gather** command:

```bash
poetry run gather -i <input_fasta> -o <output_file> -p <pam_sequence>
```

The parameters are:
- *-i*, *--ifile* - the Input File needs to be a FASTA file containing the genenome sequence. For example GRCh38 which can be downloaded from [Ensembl](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/) - *Mandatory*,
- *-o*, *--ofile* - the Output File which will be a CSV file (without headers) - *Mandatory*,
- *-p*, *--pam* - the PAM sequence which can consist of A, C, G, T and N (for any) - *Mandatory*,
- *-h*, *--help* - shows the help

For example:

```bash
poetry run gather -i Homo_sapiens.GRCh38.dna.chromosome.18.fa -o chromosome.18.csv -p "NGG"
```

### Index

To run the **Index** command:

```bash
poetry run index -i <input_csv> -a <assembly> -o <offset> -o <output_bin> -s <species> -e <species_id> -g <guide_length> -p <pam_length>
```

The parameters are:
- *-i*, *--ifile* - The Input CSV file, can be declared one or many times (see example below) - *Mandatory*,
- *-o*, *--ofile* - The Output binary file - *Mandatory*,
- *-a*, *--assembly* - The name of the genome assembly - *Mandatory*,
- *-s*, *--species* - The name of the species - *Mandatory*,
- *-f*, *--offset* - the offset after which to start the ID for CRISPRS, defaults to 0,
- *-e*, *--species_id* - The species ID, defaults to 0,
- *-g*, *--guide_length* - The length of the guide RNA, defaults to 20,
- *-p*, *--pam_length* - The length of the PAM, defaults to 3,
- *-h*, *--help* - shows the help

for example:

```bash
poetry run index -i chromosome.1.csv -i chromosome.2.csv -o guides.bin -a GRCh38 -s Human
```

The CSV input file must have the following fields:
- Chromasome Name as a string e.g. '18'
- Position Start - an integer indicating the start position offset from the start of the Chromasome (using 5' to 3' orientation),
- CRISPR Sequence - the string representing the CRISPR including the PAM,
- PAM Right? - a 1 or 0 (true or false) indicating if the PAM is on the right side of the CRISPR,
- Species ID - this is currently always set to 1.

Note that *Species ID* is a legacy field and is not used in the current version of the software.

## Testing

To run the unit tests:

```bash
poetry run pytest
```
