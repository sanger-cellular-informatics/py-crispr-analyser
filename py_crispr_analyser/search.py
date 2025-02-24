from dataclasses import dataclass
import getopt
import numpy as np
import struct
import sys
import time
import typing

from .utils import (
    FILE_VERSION,
    METADATA_FORMAT,
    sequence_to_binary_encoding,
    reverse_complement,
)

FILE_HEADER_FORMAT = "<BL"
PADDING_FORMAT = "<BBB"
HEADER_SIZE = struct.calcsize(FILE_HEADER_FORMAT)
METADATA_SIZE = struct.calcsize(METADATA_FORMAT)


@dataclass
class Metadata:
    number_of_sequences: np.uint64
    sequence_length: np.uint64
    offset: np.uint64
    species_id: np.uint8
    species_name: str
    assembly: str


def check_file_header(bytes: bytes) -> None:
    """Check the header of the file from binary data

    Args:
        bytes: The binary data to check
    """
    if len(bytes) != HEADER_SIZE:
        raise ValueError("Invalid file header length")
    values = struct.unpack(FILE_HEADER_FORMAT, bytes)
    version = values[1]
    if version != FILE_VERSION:
        raise ValueError("Invalid file version")


def get_file_metadata(bytes: bytes) -> Metadata:
    """Parse the metadata from binary data

    Args:
        bytes: The binary data to parse
    """
    if len(bytes) != METADATA_SIZE:
        raise ValueError("Invalid metadata length")
    metadata = struct.unpack(METADATA_FORMAT, bytes)
    metadata = Metadata(
        number_of_sequences=metadata[0],
        sequence_length=metadata[1],
        offset=metadata[2],
        species_id=metadata[3],
        species_name=parse_c_string(metadata[4]),
        assembly=parse_c_string(metadata[5]),
    )
    return metadata


def print_metadata(metadata: Metadata) -> None:
    """Print the metadata information

    Args:
        metadata: The Metadata object
    """
    print(f"Assembly is {metadata.assembly} ({metadata.species_name})")
    print(f"File has {metadata.number_of_sequences} sequences")
    print(f"Sequence length is {metadata.sequence_length}")
    print(f"Offset is {metadata.offset}")
    print(f"Species id is {metadata.species_id}")


def parse_c_string(data: bytes) -> str:
    """Parse a null-terminated C string from a byte array

    Args:
        data: The byte array to parse
    """
    null = data.find(b"\x00")
    if null < 0:
        return ""
    cstr = data[:null].decode("utf-8")
    return cstr


8


def get_guides(
    guidesfile_handle: typing.BinaryIO, verbose: bool = False
) -> np.ndarray:
    """Get the array of guides from the binary guides file

    Args:
        guidesfile_handle: The file handle of the guides file
        verbose: A boolean to print verbose output
    """
    if verbose:
        start = time.time()
    guidesfile_handle.seek(HEADER_SIZE)
    # read the number of sequences in the header
    number_of_guides = struct.unpack("<Q", guidesfile_handle.read(8))[0]
    guidesfile_handle.seek(
        HEADER_SIZE + METADATA_SIZE + struct.calcsize(PADDING_FORMAT)
    )
    guides = np.fromfile(guidesfile_handle, dtype=np.uint64, count=-1)
    if number_of_guides != guides.size:
        raise ValueError("Invalid number of guides")
    if verbose:
        print(f"Loading took {time.time() - start:.2f} seconds")
    return guides


def search(
    guides: np.ndarray,
    sequence: str,
    verbose: bool = False,
) -> list[int]:
    """Search for a sequence in an indexed binary file

    Args:
        guides: The numpy uint64 array of guides
        sequence: The query sequence to search for
        verbose: A boolean to print verbose output
    """
    reverse_sequence = reverse_complement(sequence)
    query_sequence = sequence_to_binary_encoding(sequence, 1)
    reverse_query_sequence = sequence_to_binary_encoding(reverse_sequence, 0)
    indices = np.where(
        (guides == query_sequence) | (guides == reverse_query_sequence)
    )
    # the binary index is 0-based,
    # so we add the offset and 1 to make it 1-based as per the db
    # this follows how we numbered the WGE index
    return [x + 1 for x in indices[0].tolist()]


def run(argv=sys.argv[1:]) -> None:
    """Run the search command from the command line."""
    inputfile = ""
    sequence = ""

    def usage() -> None:
        print(
            """Usage: poetry run search [options...]
-h, --help            Print this help message
-i, --ifile <file>    The input binary guides file
-s, --sequence <str>  The guide sequence to search for
"""
        )

    try:
        opts, args = getopt.getopt(
            argv,
            "hi:s:",
            [
                "help",
                "ifile=",
                "sequence=",
            ],
        )
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-s", "--sequence"):
            sequence = arg
    if inputfile == "" or sequence == "":
        usage()
        sys.exit(2)

    with open(inputfile, "rb") as in_file:
        check_file_header(in_file.read(HEADER_SIZE))
        print(f"Version is {FILE_VERSION}")
        metadata = get_file_metadata(in_file.read(METADATA_SIZE))
        print_metadata(metadata)
        guides = get_guides(in_file, verbose=True)
        print(f"Loaded {guides.size} sequences")
        indices = search(
            guides=guides,
            sequence=sequence,
            verbose=True,
        )
        print(f"Found {len(indices)} exact matches")
        print("Found the following matches:")
        for idx in indices:
            print(f"\t{idx + metadata.offset}")


if __name__ == "__main__":
    run()
