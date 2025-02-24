import numpy as np
import pytest
import struct

import py_crispr_analyser.search as search


class TestCheckFileHeader:
    def test_when_version_is_3(self):
        try:
            struct.pack("<BL", 1, 3)
        except ValueError as e:
            pytest.fail(f"Unexpected exception: {e}")

    def test_raises_exception_when_version_is_wrong(self):
        header = struct.pack("<BL", 1, 4)
        with pytest.raises(ValueError) as excinfo:
            search.check_file_header(header)
        assert str(excinfo.value) == "Invalid file version"

    def test_raises_exception_when_header_is_too_short(self):
        header = struct.pack("<B", 1)
        with pytest.raises(ValueError) as excinfo:
            search.check_file_header(header)
        assert str(excinfo.value) == "Invalid file header length"

    def test_raises_exception_when_header_is_too_long(self):
        header = struct.pack("<BLB", 1, 3, 1)
        with pytest.raises(ValueError) as excinfo:
            search.check_file_header(header)
        assert str(excinfo.value) == "Invalid file header length"


class TestGetFileMetadata:
    def test_when_metadata_is_correct(self):
        metadata = struct.pack(
            "<QQQB30s30s",
            1000,
            20,
            0,
            1,
            "Human".encode("utf-8"),
            "GRCh38".encode("utf-8"),
        )
        assert search.get_file_metadata(metadata) == search.Metadata(
            number_of_sequences=1000,
            sequence_length=20,
            offset=0,
            species_id=1,
            species_name="Human",
            assembly="GRCh38",
        )

    def test_raises_exception_when_metadata_is_too_short(self):
        metadata = struct.pack("<Q", 1000)
        with pytest.raises(ValueError) as excinfo:
            search.get_file_metadata(metadata)
        assert str(excinfo.value) == "Invalid metadata length"

    def test_raises_exception_when_metadata_is_too_long(self):
        metadata = struct.pack(
            "<QQQQB30s30s",
            1000,
            20,
            0,
            1000,
            1,
            "Human".encode("utf-8"),
            "GRCh38".encode("utf-8"),
        )
        with pytest.raises(ValueError) as excinfo:
            search.get_file_metadata(metadata)
        assert str(excinfo.value) == "Invalid metadata length"


@pytest.fixture
def guides_array_with_matches():
    return np.array(
        [
            np.uint64(0b10000000001111010000000011110101111011101),
            np.uint64(0b10000011110100000000111101011110111011110),
            np.uint64(0b1000100000010100101111110101001011111111),
            np.uint64(0b11000100000010100101111110101001011111111),
            np.uint64(0b1100110101010001000100010100010001010100),
            np.uint64(0b10111111010100101111111111100111111101),
            np.uint64(0b11000100000010100101111110101001011111111),
        ]
    )


@pytest.fixture
def guides_array_with_no_matches():
    ERROR_STR = np.uint64(0xFFFFFFFFFFFFFFFF)
    return np.array(
        [
            np.uint64(0b111000001010111000001010111000001010111),
            np.uint64(0b1100000101011100000101011100000101011100),
            np.uint64(0b10101110000010101110000010101110000),
            np.uint64(0b111000001010111000001010111000001010111),
            ERROR_STR,
            np.uint64(0b1100110101010001000100010100010001010100),
            np.uint64(0b100010001000101000100010101000100010001),
        ]
    )


@pytest.fixture
def guides_with_matches_file(tmp_path, guides_array_with_matches):
    d = tmp_path / "test"
    d.mkdir()
    in_file = d / "crisprs.bin"
    in_file.write_bytes(
        struct.pack(
            "<BLQQQB30s30sBBBQQQQQQQ",
            # start of file byte
            np.uint8(1),
            # file version
            np.uint(3),
            # metadata
            np.uint64(7),
            np.uint64(20),
            np.uint64(88),
            np.uint8(1),
            b"Human",
            b"GRCh38",
            # padding
            np.uint8(0),
            np.uint8(0),
            np.uint8(0),
            # array of sequences
            *guides_array_with_matches,
        )
    )
    return in_file


@pytest.fixture
def guides_with_no_matches_file(tmp_path, guides_array_with_no_matches):
    d = tmp_path / "test"
    d.mkdir()
    in_file = d / "crisprs.bin"
    in_file.write_bytes(
        struct.pack(
            "<BLQQQB30s30sBBBQQQQQQQ",
            # start of file byte
            np.uint8(1),
            # file version
            np.uint(3),
            # metadata
            np.uint64(7),
            np.uint64(20),
            np.uint64(88),
            np.uint8(1),
            b"Human",
            b"GRCh38",
            # padding
            np.uint8(0),
            np.uint8(0),
            np.uint8(0),
            # array of sequences
            *guides_array_with_no_matches,
        )
    )
    return in_file


class TestGetGuides:
    def test_when_guides_file_is_correct_with_matches(
        self, guides_with_matches_file, guides_array_with_matches
    ):
        in_file = open(guides_with_matches_file, "rb")
        assert np.array_equal(
            search.get_guides(in_file), guides_array_with_matches
        )

    def test_when_guides_file_is_correct_with_no_matches(
        self, guides_with_no_matches_file, guides_array_with_no_matches
    ):
        in_file = open(guides_with_no_matches_file, "rb")
        assert np.array_equal(
            search.get_guides(in_file), guides_array_with_no_matches
        )

    def test_when_guides_size_and_metadata_number_disagree(
        self, tmp_path, capsys
    ):
        """Test that an exception is raised when the number_of_sequences in the
        metadata does not match the number of guides in the following array"""
        d = tmp_path / "test"
        d.mkdir()
        in_file = d / "crisprs.bin"
        in_file.write_bytes(
            struct.pack(
                "<BLQQQB30s30sBBBQQ",
                # start of file byte
                np.uint8(1),
                # file version
                np.uint(3),
                # metadata
                np.uint64(7),  # number of sequences
                np.uint64(20),
                np.uint64(88),
                np.uint8(1),
                b"Human",
                b"GRCh38",
                # padding
                np.uint8(0),
                np.uint8(0),
                np.uint8(0),
                # array of sequences
                np.uint64(0b111000001010111000001010111000001010111),
                np.uint64(0b1100000101011100000101011100000101011100),
            )
        )
        in_file = open(in_file, "rb")
        with pytest.raises(ValueError) as excinfo:
            search.get_guides(in_file)
        assert str(excinfo.value) == "Invalid number of guides"


class TestSearch:
    def test_multiple_sequences_are_found(self, guides_array_with_matches):
        assert search.search(
            guides_array_with_matches, "AAAACTGGAAACTGGTTCTC"
        ) == [1, 3]

    def test_no_sequences_are_found(self, guides_array_with_no_matches):
        assert (
            search.search(guides_array_with_no_matches, "AAAACTGGAAACTGGTTCTC")
            == []
        )


class TestRun:
    def test_multiple_sequences_are_found(
        self, guides_with_matches_file, capsys
    ):
        """Test that the run function prints the correct output when multiple
        sequences are found given an offset"""
        args = [
            "-i",
            guides_with_matches_file,
            "-s",
            "AAAACTGGAAACTGGTTCTC",
        ]
        search.run(args)
        captured = capsys.readouterr()
        assert "Found 2 exact matches" in captured.out
        assert "Found the following matches:\n\t89\n\t91" in captured.out

    def test_no_sequences_are_found(self, guides_with_no_matches_file, capsys):
        args = [
            "-i",
            guides_with_no_matches_file,
            "-s",
            "AAAACTGGAAACTGGTTCTC",
        ]
        search.run(args)
        captured = capsys.readouterr()
        assert "Found 0 exact matches" in captured.out
