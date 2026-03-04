# Copyright (C) 2025-2026 Genome Research Ltd.

import numpy as np
from numba import cuda
import pytest

import py_crispr_analyser.align as align


@pytest.fixture
def query_sequence():
    # AAAACTGGAAACTGGTTCTC pam right
    return np.uint64(0b10000000001111010000000011110101111011101)


@pytest.fixture
def reverse_query_sequence():
    # GAGAACCAGTTTCCAGTTTT pam left
    return np.uint64(0b1000100000010100101111110101001011111111)


@pytest.fixture
def large_guide_list(guide_list):
    """Tile the guide list to ~1 000 000 entries for a more realistic benchmark."""
    reps = (1_000_000 + guide_list.size - 1) // guide_list.size
    return np.tile(guide_list, reps)[:1_000_000]


def bench_find_off_targets(
    benchmark, guide_list, query_sequence, reverse_query_sequence
):
    """Benchmark the pure-Python JIT off-target finder on the fixture dataset."""
    benchmark(
        align.find_off_targets,
        guide_list,
        query_sequence,
        reverse_query_sequence,
        np.uint64(0),
    )


def bench_find_off_targets_large(
    benchmark, large_guide_list, query_sequence, reverse_query_sequence
):
    """Benchmark the pure-Python JIT off-target finder on ~100 000 guides."""
    benchmark(
        align.find_off_targets,
        large_guide_list,
        query_sequence,
        reverse_query_sequence,
        np.uint64(0),
    )


def bench_find_off_targets_cpu(
    benchmark, guide_list, query_sequence, reverse_query_sequence
):
    """Benchmark the parallel-CPU off-target finder on the fixture dataset."""

    def run():
        summary = np.zeros(align.MAX_MISSMATCHES, dtype=np.uint32)
        off_target_ids_idx = np.zeros(1, dtype=np.uint32)
        off_target_ids = np.zeros(align.MAX_OFF_TARGETS, dtype=np.uint32)
        align.find_off_targets_cpu(
            guide_list,
            query_sequence,
            reverse_query_sequence,
            summary,
            off_target_ids_idx,
            off_target_ids,
            np.uint64(0),
        )

    benchmark(run)


def bench_find_off_targets_cpu_large(
    benchmark, large_guide_list, query_sequence, reverse_query_sequence
):
    """Benchmark the parallel-CPU off-target finder on ~100 000 guides."""
    # Size the off_target_ids array to the guide list to avoid overflow
    n = large_guide_list.size

    def run():
        summary = np.zeros(align.MAX_MISSMATCHES, dtype=np.uint32)
        off_target_ids_idx = np.zeros(1, dtype=np.uint32)
        off_target_ids = np.zeros(n, dtype=np.uint32)
        align.find_off_targets_cpu(
            large_guide_list,
            query_sequence,
            reverse_query_sequence,
            summary,
            off_target_ids_idx,
            off_target_ids,
            np.uint64(0),
        )

    benchmark(run)


@pytest.mark.skipif(
    cuda.is_available() is False, reason="CUDA is not available"
)
def bench_find_off_targets_kernel(
    benchmark, guide_list, query_sequence, reverse_query_sequence
):
    """Benchmark the CUDA kernel off-target finder on the fixture dataset."""
    threads_per_block = 256
    blocks_per_grid = (
        guide_list.size + threads_per_block - 1
    ) // threads_per_block
    device_guides = cuda.to_device(guide_list)

    def run():
        summary = np.zeros(align.MAX_MISSMATCHES, dtype=np.uint32)
        off_target_ids_idx = np.zeros(1, dtype=np.uint32)
        off_target_ids = np.zeros(align.MAX_OFF_TARGETS, dtype=np.uint32)
        device_summary = cuda.to_device(summary)
        device_off_target_ids_idx = cuda.to_device(off_target_ids_idx)
        device_off_target_ids = cuda.to_device(off_target_ids)
        align.find_off_targets_kernel[blocks_per_grid, threads_per_block](
            device_guides,
            query_sequence,
            reverse_query_sequence,
            device_summary,
            device_off_target_ids_idx,
            device_off_target_ids,
            np.uint64(0),
        )
        cuda.synchronize()

    benchmark(run)


@pytest.mark.skipif(
    cuda.is_available() is False, reason="CUDA is not available"
)
def bench_find_off_targets_kernel_large(
    benchmark, large_guide_list, query_sequence, reverse_query_sequence
):
    """Benchmark the CUDA kernel off-target finder on ~100 000 guides."""
    threads_per_block = 256
    blocks_per_grid = (
        large_guide_list.size + threads_per_block - 1
    ) // threads_per_block
    device_guides = cuda.to_device(large_guide_list)

    def run():
        summary = np.zeros(align.MAX_MISSMATCHES, dtype=np.uint32)
        off_target_ids_idx = np.zeros(1, dtype=np.uint32)
        off_target_ids = np.zeros(align.MAX_OFF_TARGETS, dtype=np.uint32)
        device_summary = cuda.to_device(summary)
        device_off_target_ids_idx = cuda.to_device(off_target_ids_idx)
        device_off_target_ids = cuda.to_device(off_target_ids)
        align.find_off_targets_kernel[blocks_per_grid, threads_per_block](
            device_guides,
            query_sequence,
            reverse_query_sequence,
            device_summary,
            device_off_target_ids_idx,
            device_off_target_ids,
            np.uint64(0),
        )
        cuda.synchronize()

    benchmark(run)
