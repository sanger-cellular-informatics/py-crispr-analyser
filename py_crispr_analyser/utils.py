import numpy as np

COMPLEMENT_MAP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
ENCODING_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
ERROR_STR = np.uint64(0xFFFFFFFFFFFFFFFF)
FILE_VERSION = np.uint16(3)
METADATA_FORMAT = "<QQQB30s30s"

def sequence_to_binary_encoding(sequence: str, pam_right: int) -> np.uint64:
    """Convert a string DNA sequence to bits accounting for pam right or left

    Args:
        sequence: The string DNA sequence to convert
        pam_right: An integer indicating if PAM is on the right
    """
    bits = np.uint64(pam_right)
    for character in sequence:
        if ENCODING_MAP[character] == 4:
            bits = ERROR_STR  # set the error flag if we encounter an 'N'
            break
        else:
            bits <<= 2  # shift left to make room for new char
            bits |= ENCODING_MAP[character]  # add our new char to the end
    return bits

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Args:
        sequence: The string DNA sequence to reverse complement.
    """
    return "".join([COMPLEMENT_MAP[base] for base in reversed(sequence)])
