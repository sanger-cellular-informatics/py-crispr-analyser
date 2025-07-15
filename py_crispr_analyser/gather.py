# Copyright (C) 2025 Genome Research Ltd.

import collections
import csv
import getopt
import re
import sys
import time

from .utils import reverse_complement

GUIDE_RNA_LENGTH = 20


def match_pam(
    dna_sequence: str,
    pam_sequence: str,
    pam_on_right: bool,
    legacy_mode: bool = False,
) -> bool:
    """Check if the DNA sequence has a PAM sequence match.

    :param dna_sequence: The string DNA sequence to check.
    :param pam_sequence: The string PAM sequence to match.
    :param pam_on_right: A boolean indicating if PAM sequence is on the right.
    :param legacy_mode: A boolean indicating if non-ACGT chars allowed in PAM
        region of the DNA sequence. Default is False.
    :return: True if the PAM sequence matches, False otherwise.
    """
    start = len(dna_sequence) - len(pam_sequence) if pam_on_right else 0
    for i in range(len(pam_sequence)):
        if legacy_mode is False and dna_sequence[start + i] not in "ACGT":
            return False
        if pam_sequence[i] == "N":
            continue
        if dna_sequence[start + i] != pam_sequence[i]:
            return False
    return True


def gather(
    inputfile: str,
    outputfile: str,
    pam: str,
    verbose: bool = False,
    legacy_mode: bool = False,
) -> None:
    """Run the CRISPR gatherer.

    :param inputfile: The input FASTA file containing DNA sequences.
    :param outputfile: The output CSV file to write results to.
    :param pam: The string PAM sequence to search for e.g. "NGG".
    :param verbose: A boolean indicating if verbose output is enabled.
        Default is False.
    :param legacy_mode: A boolean indicating that species ID column
        is added to CSV file (always equalling 1). Default is False.
    :return: None
    """
    if verbose:
        start = time.time()
    chromosome = ""
    crispr_count = 0
    position = 0
    buffer = collections.deque(maxlen=len(pam) + GUIDE_RNA_LENGTH)

    with open(inputfile, "r") as infile:
        with open(outputfile, "w", newline="") as outfile:
            csvwriter = csv.writer(outfile)
            for line in infile:
                if len(line) == 0:
                    continue
                if line[0] == ">":
                    chromosome = re.search(
                        r">(.*?) dna:chromosome", line
                    ).group(1)
                    if verbose:
                        print(f"Processing chromosome {chromosome}...")
                    position = 0
                    buffer.clear()
                    continue
                else:
                    line = line.strip()
                    for base in line:
                        buffer.append(base)
                        if len(buffer) < len(pam) + GUIDE_RNA_LENGTH:
                            continue
                        position += 1
                        if match_pam(
                            dna_sequence=buffer,
                            pam_sequence=reverse_complement(pam),
                            pam_on_right=False,
                        ):
                            output = (
                                [chromosome, position, "".join(buffer), 0, 1]
                                if legacy_mode
                                else [chromosome, position, "".join(buffer), 0]
                            )
                            csvwriter.writerow(output)
                            crispr_count += 1
                        if match_pam(
                            dna_sequence=buffer,
                            pam_sequence=pam,
                            pam_on_right=True,
                        ):
                            output = (
                                [chromosome, position, "".join(buffer), 1, 1]
                                if legacy_mode
                                else [chromosome, position, "".join(buffer), 1]
                            )
                            csvwriter.writerow(output)
                            crispr_count += 1
    if verbose:
        end = time.time()
        print(f"Gathered {crispr_count} CRISPRs in {end - start} seconds.")


def run(argv=sys.argv[1:]):
    """Run the CRISPR gatherer from the command line."""
    inputfile = ""
    outputfile = ""
    pam = ""

    def usage():
        print(
            """Usage: crispr_analyser_gather [options...]
-h, --help           Print this help message
-i, --ifile <file>   The input FASTA file
-o, --ofile <file>   The output file
-p, --pam <pam seq>  The PAM sequence to search for
"""
        )

    try:
        opts, args = getopt.getopt(
            argv, "hi:o:p:", ["help", "ifile=", "ofile=", "pam="]
        )
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-p", "--pam"):
            pam = arg
        else:
            print("Unhandled Option")
            usage()
            sys.exit(2)
    if inputfile == "" or outputfile == "" or pam == "":
        usage()
        sys.exit(2)

    gather(inputfile, outputfile, pam, verbose=True, legacy_mode=True)
