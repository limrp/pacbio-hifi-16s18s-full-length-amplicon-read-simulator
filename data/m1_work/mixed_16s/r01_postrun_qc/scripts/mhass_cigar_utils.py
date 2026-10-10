#!/usr/bin/env python3

"""
Utilities for parsing edlib extended CIGAR strings and projecting
reference boundaries onto observed-query coordinates.

Coordinate convention:
- reference = expected MHASS construct R
- query     = canonicalized observed CCS read Q
- coordinates are 0-based, half-open where intervals are used

edlib extended CIGAR operations:
    =  match
    X  substitution
    I  insertion in query relative to reference
    D  deletion from query relative to reference
"""

from __future__ import annotations

import re


_CIGAR_TOKEN = re.compile(r"(\d+)([=XID])")


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """Parse an edlib extended CIGAR into (length, operation) tuples."""

    if not cigar:
        raise ValueError("CIGAR is empty.")

    tokens = _CIGAR_TOKEN.findall(cigar)

    reconstructed = "".join(
        f"{length}{operation}"
        for length, operation in tokens
    )

    if reconstructed != cigar:
        raise ValueError(f"Unsupported or malformed CIGAR: {cigar}")

    parsed: list[tuple[int, str]] = []

    for length_text, operation in tokens:
        length = int(length_text)

        if length <= 0:
            raise ValueError(
                f"CIGAR operation has non-positive length: "
                f"{length}{operation}"
            )

        parsed.append((length, operation))

    return parsed


def _reference_region(
    reference_position: int,
    target_start: int,
    target_end: int,
) -> str:
    """Classify a reference base position into left/target/right."""

    if reference_position < target_start:
        return "left_wrapper"

    if reference_position < target_end:
        return "target"

    return "right_wrapper"


def walk_cigar(
    cigar: str,
    *,
    reference_length: int,
    query_length: int,
    target_start: int,
    target_end: int,
) -> dict[str, int | bool]:
    """
    Walk an edlib extended CIGAR and project biological boundaries.

    The reference is the known complete MHASS construct R.
    The query is the canonicalized observed CCS sequence Q.

    Insertions exactly at target_start or target_end are recorded
    separately because their assignment to wrapper vs biological target
    is intrinsically ambiguous from alignment alone.
    """

    if not (0 <= target_start <= target_end <= reference_length):
        raise ValueError(
            "Invalid target interval: "
            f"[{target_start}, {target_end}) for "
            f"reference length {reference_length}"
        )

    counts = {
        "left_wrapper_substitutions": 0,
        "left_wrapper_insertions": 0,
        "left_wrapper_deletions": 0,
        "target_substitutions": 0,
        "target_insertions": 0,
        "target_deletions": 0,
        "right_wrapper_substitutions": 0,
        "right_wrapper_insertions": 0,
        "right_wrapper_deletions": 0,
        "left_boundary_insertions": 0,
        "right_boundary_insertions": 0,
    }

    boundary_hits = {
        target_start: [],
        target_end: [],
    }

    reference_position = 0
    query_position = 0

    def record_boundary_vertex() -> None:
        if reference_position in boundary_hits:
            boundary_hits[reference_position].append(query_position)

    # Alignment begins at vertex (0, 0).
    record_boundary_vertex()

    for run_length, operation in parse_cigar(cigar):
        for _ in range(run_length):

            if operation == "=":
                reference_position += 1
                query_position += 1

            elif operation == "X":
                region = _reference_region(
                    reference_position,
                    target_start,
                    target_end,
                )

                counts[f"{region}_substitutions"] += 1

                reference_position += 1
                query_position += 1

            elif operation == "D":
                region = _reference_region(
                    reference_position,
                    target_start,
                    target_end,
                )

                counts[f"{region}_deletions"] += 1

                reference_position += 1

            elif operation == "I":
                if reference_position == target_start:
                    counts["left_boundary_insertions"] += 1

                elif reference_position == target_end:
                    counts["right_boundary_insertions"] += 1

                else:
                    if reference_position < target_start:
                        region = "left_wrapper"

                    elif reference_position < target_end:
                        region = "target"

                    else:
                        region = "right_wrapper"

                    counts[f"{region}_insertions"] += 1

                query_position += 1

            else:
                raise RuntimeError(
                    f"Unexpected CIGAR operation: {operation}"
                )

            # A CIGAR path is a series of alignment vertices.
            # Multiple query coordinates at the same reference boundary
            # indicate an insertion exactly at that boundary.
            record_boundary_vertex()

    if reference_position != reference_length:
        raise ValueError(
            "CIGAR/reference length inconsistency: "
            f"walked {reference_position}, "
            f"expected {reference_length}"
        )

    if query_position != query_length:
        raise ValueError(
            "CIGAR/query length inconsistency: "
            f"walked {query_position}, "
            f"expected {query_length}"
        )

    start_hits = boundary_hits[target_start]
    end_hits = boundary_hits[target_end]

    if not start_hits:
        raise RuntimeError(
            f"Target-start boundary {target_start} was never reached."
        )

    if not end_hits:
        raise RuntimeError(
            f"Target-end boundary {target_end} was never reached."
        )

    alignment_start_min = min(start_hits)
    alignment_start_max = max(start_hits)

    alignment_end_min = min(end_hits)
    alignment_end_max = max(end_hits)

    left_boundary_ambiguous = (
        alignment_start_min != alignment_start_max
    )

    right_boundary_ambiguous = (
        alignment_end_min != alignment_end_max
    )

    left_wrapper_length = target_start
    right_wrapper_length = reference_length - target_end

    fixed_target_start = left_wrapper_length
    fixed_target_end = query_length - right_wrapper_length

    result: dict[str, int | bool] = {
        **counts,
        "fixed_target_start": fixed_target_start,
        "fixed_target_end": fixed_target_end,
        "alignment_target_start_min": alignment_start_min,
        "alignment_target_start_max": alignment_start_max,
        "alignment_target_end_min": alignment_end_min,
        "alignment_target_end_max": alignment_end_max,
        "left_boundary_ambiguous": left_boundary_ambiguous,
        "right_boundary_ambiguous": right_boundary_ambiguous,
        "start_shift_min": (
            alignment_start_min - fixed_target_start
        ),
        "start_shift_max": (
            alignment_start_max - fixed_target_start
        ),
        "end_shift_min": (
            alignment_end_min - fixed_target_end
        ),
        "end_shift_max": (
            alignment_end_max - fixed_target_end
        ),
    }

    return result
