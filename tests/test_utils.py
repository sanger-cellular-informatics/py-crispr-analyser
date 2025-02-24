import numpy as np
import pytest

import py_crispr_analyser.utils as utils

class TestSequenceToBinaryEncoding:
    def test_when_pam_is_right(self):
        assert utils.sequence_to_binary_encoding(
            sequence="ACGT", pam_right=1
        ) == np.uint64(0b100011011)

    def test_with_pam_right_false(self):
        assert utils.sequence_to_binary_encoding(sequence="ACGT", pam_right=0) == np.uint64(0b11011)

    def test_with_error_flag(self):
        assert utils.sequence_to_binary_encoding(sequence="ACGN", pam_right=0) == utils.ERROR_STR


def test_reverse_complement():
    assert utils.reverse_complement("ATCGN") == "NCGAT"


