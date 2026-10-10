#!/usr/bin/env python3

"""
Utilities for orienting and slicing observed FASTQ reads during
MHASS controlled-read preparation.

Important:
- Sequence bases are never corrected using the truth template.
- Reverse-oriented reads are reverse-complemented.
- FASTQ qualities are reversed, not complemented.
- Sequence and quality coordinates must always remain synchronized.
"""

from __future__ import annotations

from Bio.Seq import Seq


def orient_sequence_and_qualities(
    sequence: str,
    qualities: list[int],
    orientation: str,
) -> tuple[str, list[int]]:
    """
    Return sequence and qualities in canonical truth-template orientation.

    orientation must be:
        "as_read"
        "revcomp"
    """

    sequence = sequence.upper()
    qualities = list(qualities)

    if len(sequence) != len(qualities):
        raise ValueError(
            "Sequence/quality length mismatch before orientation: "
            f"sequence={len(sequence)}, qualities={len(qualities)}"
        )

    if orientation == "as_read":
        canonical_sequence = sequence
        canonical_qualities = qualities.copy()

    elif orientation == "revcomp":
        canonical_sequence = str(
            Seq(sequence).reverse_complement()
        )
        canonical_qualities = list(reversed(qualities))

    else:
        raise ValueError(
            f"Unsupported orientation: {orientation}"
        )

    if len(canonical_sequence) != len(canonical_qualities):
        raise RuntimeError(
            "Sequence/quality length mismatch after orientation."
        )

    return canonical_sequence, canonical_qualities


def slice_sequence_and_qualities(
    sequence: str,
    qualities: list[int],
    start: int,
    end: int,
) -> tuple[str, list[int]]:
    """
    Slice sequence and FASTQ qualities using the same 0-based,
    half-open interval [start, end).
    """

    if len(sequence) != len(qualities):
        raise ValueError(
            "Sequence/quality length mismatch before slicing: "
            f"sequence={len(sequence)}, qualities={len(qualities)}"
        )

    if start < 0:
        raise ValueError(
            f"Invalid negative start coordinate: {start}"
        )

    if end < start:
        raise ValueError(
            f"Invalid interval: start={start}, end={end}"
        )

    if end > len(sequence):
        raise ValueError(
            f"End coordinate exceeds sequence length: "
            f"end={end}, length={len(sequence)}"
        )

    sliced_sequence = sequence[start:end]
    sliced_qualities = qualities[start:end]

    if len(sliced_sequence) != len(sliced_qualities):
        raise RuntimeError(
            "Sequence/quality length mismatch after slicing."
        )

    return sliced_sequence, sliced_qualities
