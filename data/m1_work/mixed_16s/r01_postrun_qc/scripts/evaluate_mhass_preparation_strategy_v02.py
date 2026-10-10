#!/usr/bin/env python3

"""
Evaluate MHASS read provenance, orientation, and biological-target
boundary projection for controlled simulated reads.

This diagnostic does NOT:
- trim reads
- modify reads
- correct sequencing errors
- write FASTQ files

It reconstructs the expected complete MHASS construct for each realized
read, determines canonical orientation using global edit distance, and
uses the selected global alignment path to evaluate wrapper errors and
project biological-target boundaries onto observed-read coordinates.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import edlib
from Bio import SeqIO
from Bio.Seq import Seq
from mhass_cigar_utils import walk_cigar


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate MHASS provenance and orientation using global edit-distance "
            "alignment against the known expected construct."
        )
    )

    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)

    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Process only the first N FASTQ reads. Default: all reads.",
    )

    return parser.parse_args()


def load_truth(path: Path) -> dict[str, str]:
    truth: dict[str, str] = {}

    for record in SeqIO.parse(path, "fasta"):
        if record.id in truth:
            raise ValueError(f"Duplicate truth FASTA ID: {record.id}")

        truth[record.id] = str(record.seq).upper()

    if not truth:
        raise ValueError("Truth FASTA contains no sequences.")

    return truth


def load_barcodes(path: Path) -> dict[str, dict[str, str]]:
    barcodes: dict[str, dict[str, str]] = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "SampleID",
            "BarcodeID",
            "ForwardBarcode",
            "RevCompReverse",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Barcode table is missing required columns: "
                f"{sorted(missing)}. "
                f"Found columns: {reader.fieldnames}"
            )

        for row in reader:
            sample = row["SampleID"]

            if sample in barcodes:
                raise ValueError(f"Duplicate barcode entry for sample: {sample}")

            barcodes[sample] = {
                "barcode_id": row["BarcodeID"],
                "forward": row["ForwardBarcode"].upper(),
                "revcomp_reverse": row["RevCompReverse"].upper(),
            }

    return barcodes


def load_mapping(path: Path) -> dict[str, tuple[str, str, str]]:
    """
    Return:
        template_group -> (truth_template_id, sample_id, barcode_id)

    Multiple intended-molecule rows may point to the same TemplateFile.
    They must all agree on ASVID, SampleID, and BarcodeID.
    """

    mapping: dict[str, tuple[str, str, str]] = {}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"ASVID", "SampleID", "BarcodeID", "TemplateFile"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Mapping table is missing required columns: "
                f"{sorted(missing)}. "
                f"Found columns: {reader.fieldnames}"
            )

        for row in reader:
            group = row["TemplateFile"]

            if group.endswith(".fasta"):
                group = group[:-6]

            value = (
                row["ASVID"],
                row["SampleID"],
                row["BarcodeID"],
            )

            if group in mapping and mapping[group] != value:
                raise ValueError(
                    "Conflicting provenance for template group "
                    f"{group}: {mapping[group]} vs {value}"
                )

            mapping[group] = value

    return mapping


def build_expected_construct(
    truth_sequence: str,
    forward_barcode: str,
    revcomp_reverse_barcode: str,
) -> str:
    return (
        "A"
        + forward_barcode
        + truth_sequence
        + revcomp_reverse_barcode
        + "A"
    )


def edit_distance(query: str, reference: str) -> int:
    result = edlib.align(
        query,
        reference,
        mode="NW",
        task="distance",
    )

    distance = result["editDistance"]

    if distance < 0:
        raise RuntimeError("edlib failed to produce a global edit distance.")

    return distance

def alignment_path(
    query: str,
    reference: str,
) -> tuple[int, str]:
    result = edlib.align(
        query,
        reference,
        mode="NW",
        task="path",
    )

    distance = result["editDistance"]
    cigar = result["cigar"]

    if distance < 0:
        raise RuntimeError(
            "edlib failed to produce a global alignment."
        )

    if not cigar:
        raise RuntimeError(
            "edlib did not return a CIGAR alignment path."
        )

    return distance, cigar
    
def main() -> int:
    args = parse_args()

    truth = load_truth(args.truth)
    barcodes = load_barcodes(args.barcodes)
    mapping = load_mapping(args.mapping)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "read_index",
        "read_id",
        "template_group",
        "truth_template_id",
        "sample_id",
        "barcode_id",
        "observed_length",
        "expected_length",
        "length_delta",
        "d_as_read",
        "d_revcomp",
        "selected_orientation",
        "best_edit_distance",
        "orientation_margin",
        "best_edit_rate",
        
        "alignment_cigar",

        "left_wrapper_substitutions",
        "left_wrapper_insertions",
        "left_wrapper_deletions",

        "target_substitutions",
        "target_insertions",
        "target_deletions",

        "right_wrapper_substitutions",
        "right_wrapper_insertions",
        "right_wrapper_deletions",

        "left_boundary_insertions",
        "right_boundary_insertions",

        "fixed_target_start",
        "fixed_target_end",

        "alignment_target_start_min",
        "alignment_target_start_max",
        "alignment_target_end_min",
        "alignment_target_end_max",

        "start_shift_min",
        "start_shift_max",
        "end_shift_min",
        "end_shift_max",

        "left_boundary_ambiguous",
        "right_boundary_ambiguous",
        "boundary_status",
    ]

    processed = 0
    forward_count = 0
    reverse_count = 0
    ambiguous_count = 0

    with args.output.open("w", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for read_index, record in enumerate(
            SeqIO.parse(args.fastq, "fastq"),
            start=1,
        ):
            if args.max_reads is not None and read_index > args.max_reads:
                break

            read_id = record.id
            template_group = read_id.split("/", 1)[0]

            if template_group not in mapping:
                raise ValueError(
                    f"Unknown FASTQ template group: {template_group}"
                )

            truth_template_id, sample_id, barcode_id = mapping[template_group]

            if truth_template_id not in truth:
                raise ValueError(
                    f"Truth template missing from FASTA: {truth_template_id}"
                )

            if sample_id not in barcodes:
                raise ValueError(
                    f"Sample missing from barcode table: {sample_id}"
                )

            barcode_info = barcodes[sample_id]

            if barcode_info["barcode_id"] != barcode_id:
                raise ValueError(
                    f"Barcode mismatch for {template_group}: "
                    f"mapping={barcode_id}, "
                    f"sample table={barcode_info['barcode_id']}"
                )

            truth_sequence = truth[truth_template_id]

            expected = build_expected_construct(
                truth_sequence=truth_sequence,
                forward_barcode=barcode_info["forward"],
                revcomp_reverse_barcode=barcode_info["revcomp_reverse"],
            )
            
            # Reference coordinates of the biological interval in expected construct R.
            target_start = 1 + len(barcode_info["forward"])
            target_end = target_start + len(truth_sequence)

            expected_right_wrapper_length = (
                len(barcode_info["revcomp_reverse"]) + 1
            )

            if len(expected) - target_end != expected_right_wrapper_length:
                raise RuntimeError(
                    f"Expected-construct structure inconsistency for {read_id}: "
                    f"right wrapper length="
                    f"{len(expected) - target_end}, "
                    f"expected={expected_right_wrapper_length}"
                )
    
            q_as_read = str(record.seq).upper()
            q_revcomp = str(Seq(q_as_read).reverse_complement())

            d_as_read = edit_distance(q_as_read, expected)
            d_revcomp = edit_distance(q_revcomp, expected)

            if d_as_read < d_revcomp:
                orientation = "as_read"
                best_distance = d_as_read
                canonical_query = q_as_read
                forward_count += 1

            elif d_revcomp < d_as_read:
                orientation = "revcomp"
                best_distance = d_revcomp
                canonical_query = q_revcomp
                reverse_count += 1

            else:
                orientation = "ambiguous"
                best_distance = d_as_read
                canonical_query = None
                ambiguous_count += 1

            if canonical_query is None:
                raise RuntimeError(
                    f"Cannot perform boundary projection for "
                    f"ambiguous read: {read_id}"
                )
            
            orientation_margin = abs(d_as_read - d_revcomp)
            best_edit_rate = best_distance / len(expected)
            
            path_distance, cigar = alignment_path(
                canonical_query,
                expected,
            )

            if path_distance != best_distance:
                raise RuntimeError(
                    f"Distance/path inconsistency for {read_id}: "
                    f"orientation distance={best_distance}, "
                    f"path distance={path_distance}"
                )
    
            boundary_qc = walk_cigar(
                cigar,
                reference_length=len(expected),
                query_length=len(canonical_query),
                target_start=target_start,
                target_end=target_end,
            )

            if (
                boundary_qc["left_boundary_ambiguous"]
                or boundary_qc["right_boundary_ambiguous"]
            ):
                boundary_status = "ambiguous"
            else:
                boundary_status = "ok"

            writer.writerow(
                {
                    "read_index": read_index,
                    "read_id": read_id,
                    "template_group": template_group,
                    "truth_template_id": truth_template_id,
                    "sample_id": sample_id,
                    "barcode_id": barcode_id,
                    "observed_length": len(q_as_read),
                    "expected_length": len(expected),
                    "length_delta": len(q_as_read) - len(expected),
                    "d_as_read": d_as_read,
                    "d_revcomp": d_revcomp,
                    "selected_orientation": orientation,
                    "best_edit_distance": best_distance,
                    "orientation_margin": orientation_margin,
                    "best_edit_rate": f"{best_edit_rate:.8f}",
                    
                    "alignment_cigar": cigar,

                    "left_wrapper_substitutions":
                        boundary_qc["left_wrapper_substitutions"],
                    "left_wrapper_insertions":
                        boundary_qc["left_wrapper_insertions"],
                    "left_wrapper_deletions":
                        boundary_qc["left_wrapper_deletions"],

                    "target_substitutions":
                        boundary_qc["target_substitutions"],
                    "target_insertions":
                        boundary_qc["target_insertions"],
                    "target_deletions":
                        boundary_qc["target_deletions"],

                    "right_wrapper_substitutions":
                        boundary_qc["right_wrapper_substitutions"],
                    "right_wrapper_insertions":
                        boundary_qc["right_wrapper_insertions"],
                    "right_wrapper_deletions":
                        boundary_qc["right_wrapper_deletions"],

                    "left_boundary_insertions":
                        boundary_qc["left_boundary_insertions"],
                    "right_boundary_insertions":
                        boundary_qc["right_boundary_insertions"],

                    "fixed_target_start":
                        boundary_qc["fixed_target_start"],
                    "fixed_target_end":
                        boundary_qc["fixed_target_end"],

                    "alignment_target_start_min":
                        boundary_qc["alignment_target_start_min"],
                    "alignment_target_start_max":
                        boundary_qc["alignment_target_start_max"],
                    "alignment_target_end_min":
                        boundary_qc["alignment_target_end_min"],
                    "alignment_target_end_max":
                        boundary_qc["alignment_target_end_max"],

                    "start_shift_min":
                        boundary_qc["start_shift_min"],
                    "start_shift_max":
                        boundary_qc["start_shift_max"],
                    "end_shift_min":
                        boundary_qc["end_shift_min"],
                    "end_shift_max":
                        boundary_qc["end_shift_max"],

                    "left_boundary_ambiguous":
                        boundary_qc["left_boundary_ambiguous"],
                    "right_boundary_ambiguous":
                        boundary_qc["right_boundary_ambiguous"],

                    "boundary_status":
                        boundary_status,
                }
            )

            processed += 1

    print(f"processed_reads\t{processed}")
    print(f"as_read\t{forward_count}")
    print(f"revcomp\t{reverse_count}")
    print(f"ambiguous\t{ambiguous_count}")
    print(f"output\t{args.output}")

    if processed == 0:
        print("ERROR: no FASTQ records were processed.", file=sys.stderr)
        return 1

    if ambiguous_count > 0:
        print(
            "WARNING: ambiguous orientation calls were observed.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
