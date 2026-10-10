#!/usr/bin/env python3

"""
Independently validate a complete prepared MHASS dataset.

PURPOSE
-------
This script validates the output of prepare_mhass_sample_reads.py
without reusing the production sequence/quality manipulation helpers.

For every realized read, it compares four sources:

    1. Raw MHASS FASTQ
    2. Validated v0.2 diagnostic TSV
    3. Production truth_reads.tsv
    4. Actual prepared per-sample FASTQ

The validator checks:

- Exact read-ID set equality.
- Production orientation vs validated v0.2 orientation.
- Production target coordinates vs validated v0.2 coordinates.
- Boundary ambiguity flags and boundary policy.
- Sample assignment.
- output_fastq assignment.
- Raw and prepared sequence lengths.
- FASTQ sequence/quality synchronization.
- Independently reconstructed prepared sequence.
- Independently reconstructed prepared quality scores.
- Expected marker/dataset/replicate/seed metadata.

IMPORTANT
---------
The validator does NOT call:

    mhass_fastq_utils.py
    prepare_mhass_sample_reads.py

for sequence orientation or slicing.

That independence helps detect bugs in the production implementation.
"""

# Allow modern type annotations consistently.
from __future__ import annotations

# argparse handles command-line arguments.
import argparse

# csv reads our tab-separated diagnostic and provenance files.
import csv

# Counter makes it convenient to summarize categories such as
# sample, orientation, and boundary status.
from collections import Counter

# Path provides explicit and safe filesystem path handling.
from pathlib import Path

# SeqIO parses FASTQ records and their quality scores.
from Bio import SeqIO

# Seq performs reverse complementation independently from the
# production FASTQ helper.
from Bio.Seq import Seq

# SeqRecord is used only for type annotations and clarity.
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    """
    Define all command-line inputs required for full validation.
    """

    # Create the command-line parser.
    parser = argparse.ArgumentParser(
        description=(
            "Independently validate a complete prepared MHASS dataset "
            "against raw reads and validated v0.2 diagnostics."
        )
    )

    # Raw combined MHASS FASTQ containing all realized R01 reads.
    parser.add_argument(
        "--raw-fastq",
        required=True,
        type=Path,
        help="Original MHASS combined_reads.fastq.",
    )

    # Complete validated v0.2 diagnostic TSV for R01.
    parser.add_argument(
        "--diagnostic",
        required=True,
        type=Path,
        help="Complete validated v0.2 diagnostic TSV.",
    )

    # Production provenance table generated during preparation.
    parser.add_argument(
        "--truth-reads",
        required=True,
        type=Path,
        help="Production truth_reads.tsv.",
    )

    # Directory containing S01.fastq, S02.fastq, S03.fastq, etc.
    parser.add_argument(
        "--prepared-dir",
        required=True,
        type=Path,
        help="Directory containing prepared per-sample FASTQ files.",
    )

    # Expected experiment metadata.
    #
    # We validate these fields so a technically correct read cannot
    # accidentally be labelled as the wrong marker or replicate.
    parser.add_argument(
        "--expected-marker",
        required=True,
        help="Expected marker, for example 16S.",
    )

    parser.add_argument(
        "--expected-dataset-id",
        required=True,
        help="Expected dataset identifier, for example toyctrl16.",
    )

    parser.add_argument(
        "--expected-simulation-replicate",
        required=True,
        help="Expected simulation replicate, for example R01.",
    )

    parser.add_argument(
        "--expected-master-seed",
        required=True,
        type=int,
        help="Expected master simulation seed.",
    )

    # Parse and return the values supplied on the command line.
    return parser.parse_args()


def parse_boolean(value: str) -> bool:
    """
    Convert textual TSV boolean values to real Python booleans.

    Accepted representations:

        True / true
        False / false

    Anything else is considered invalid metadata.
    """

    # Remove surrounding whitespace and normalize case.
    normalized = value.strip().lower()

    # Recognize true.
    if normalized == "true":
        return True

    # Recognize false.
    if normalized == "false":
        return False

    # Refuse to guess if the value has an unexpected representation.
    raise ValueError(
        f"Invalid boolean value: {value!r}"
    )


def load_tsv_by_read_id(
    path: Path,
    required_columns: set[str],
) -> dict[str, dict[str, str]]:
    """
    Load a TSV as:

        read_id -> complete TSV row

    The function also validates:

    - all required columns exist;
    - every read_id is non-empty;
    - no read_id occurs twice.
    """

    # Dictionary in which each read ID will map to one complete row.
    rows: dict[str, dict[str, str]] = {}

    # Open the file without newline transformation.
    with path.open(newline="") as handle:

        # Read fields by their header names.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Determine which columns are available.
        observed_columns = set(
            reader.fieldnames or []
        )

        # Determine whether any required columns are missing.
        missing_columns = (
            required_columns - observed_columns
        )

        # Fail immediately if the TSV schema is incomplete.
        if missing_columns:
            raise ValueError(
                f"{path} is missing required columns: "
                f"{sorted(missing_columns)}"
            )

        # Process every data row.
        for row in reader:

            # read_id is our unique primary key.
            read_id = row["read_id"]

            # Empty identifiers are never acceptable.
            if not read_id:
                raise ValueError(
                    f"Empty read_id in TSV: {path}"
                )

            # Duplicate IDs would violate our required 1:1 relationship.
            if read_id in rows:
                raise ValueError(
                    f"Duplicate read_id in {path}: {read_id}"
                )

            # Store the complete row.
            rows[read_id] = row

    # Return the indexed TSV.
    return rows


def load_fastq_by_read_id(
    path: Path,
) -> dict[str, SeqRecord]:
    """
    Load one FASTQ as:

        read_id -> SeqRecord

    Duplicate read IDs are forbidden.
    """

    # Result dictionary.
    records: dict[str, SeqRecord] = {}

    # Parse records using a real FASTQ parser.
    for record in SeqIO.parse(
        path,
        "fastq",
    ):

        # Duplicate IDs would make read-level validation ambiguous.
        if record.id in records:
            raise ValueError(
                f"Duplicate FASTQ read ID in {path}: "
                f"{record.id}"
            )

        # Store the complete sequence + quality record.
        records[record.id] = record

    # Return all records.
    return records


def load_prepared_fastqs(
    prepared_dir: Path,
) -> tuple[
    dict[str, SeqRecord],
    dict[str, str],
    Counter,
]:
    """
    Load every *.fastq file in the prepared output directory.

    Returns three objects:

        records
            read_id -> SeqRecord

        locations
            read_id -> filename containing that read

        per_file_counts
            filename -> number of FASTQ records

    This lets us validate not only a read's contents, but also that it
    was written to the correct sample FASTQ file.
    """

    # Find all FASTQ files immediately inside the prepared directory.
    fastq_paths = sorted(
        prepared_dir.glob("*.fastq")
    )

    # A prepared dataset without FASTQ files is invalid.
    if not fastq_paths:
        raise ValueError(
            f"No *.fastq files found in {prepared_dir}"
        )

    # Global read lookup across every sample.
    records: dict[str, SeqRecord] = {}

    # Record which physical FASTQ file contains each read.
    locations: dict[str, str] = {}

    # Count records independently for each output FASTQ.
    per_file_counts: Counter = Counter()

    # Process one prepared FASTQ file at a time.
    for path in fastq_paths:

        # Parse every record in the current file.
        for record in SeqIO.parse(
            path,
            "fastq",
        ):

            # A read may occur in exactly one prepared FASTQ.
            #
            # This catches both duplicates within one file and the same
            # read accidentally appearing in multiple sample FASTQs.
            if record.id in records:
                raise ValueError(
                    "Duplicate prepared read ID across FASTQs: "
                    f"{record.id}"
                )

            # Save the complete record.
            records[record.id] = record

            # Save only the filename, such as S01.fastq.
            locations[record.id] = path.name

            # Count this record for its physical FASTQ file.
            per_file_counts[path.name] += 1

    # Return all three validation structures.
    return records, locations, per_file_counts


def report_set_difference(
    label: str,
    values: set[str],
    max_examples: int = 5,
) -> None:
    """
    Print the size of a set difference.

    If the difference is non-empty, print a few example read IDs so a
    failure can be diagnosed without flooding the terminal with
    thousands of lines.
    """

    # Always print the count.
    print(
        f"{label}\t{len(values)}"
    )

    # Show only a small deterministic sample if there are failures.
    for read_id in sorted(values)[:max_examples]:
        print(
            f"{label}_example\t{read_id}"
        )


def main() -> int:
    """
    Perform full independent R01 validation.
    """

    # Read command-line arguments.
    args = parse_args()

    # Define the production truth_reads.tsv columns needed below.
    truth_required_columns = {
        "read_id",
        "marker",
        "dataset_id",
        "sample_id",
        "simulation_replicate",
        "master_seed",
        "truth_template_id",
        "barcode_id",
        "mhass_template_file",
        "input_orientation",
        "boundary_status",
        "left_boundary_ambiguous",
        "right_boundary_ambiguous",
        "target_start",
        "target_end",
        "raw_length",
        "prepared_length",
        "output_fastq",
    }

    # Define the validated v0.2 diagnostic columns needed below.
    diagnostic_required_columns = {
        "read_id",
        "selected_orientation",
        "alignment_target_start_min",
        "alignment_target_end_max",
        "left_boundary_ambiguous",
        "right_boundary_ambiguous",
        "observed_length",
    }

    # Load production metadata.
    truth_rows = load_tsv_by_read_id(
        args.truth_reads,
        truth_required_columns,
    )

    # Load complete validated v0.2 diagnostics.
    diagnostic_rows = load_tsv_by_read_id(
        args.diagnostic,
        diagnostic_required_columns,
    )

    # Load every raw MHASS read.
    raw_records = load_fastq_by_read_id(
        args.raw_fastq
    )

    # Load all prepared sample FASTQs.
    (
        prepared_records,
        prepared_locations,
        per_file_counts,
    ) = load_prepared_fastqs(
        args.prepared_dir
    )

    # ------------------------------------------------------------------
    # DATASET-LEVEL COUNTS
    # ------------------------------------------------------------------

    # Print basic source sizes first.
    print(
        "raw_fastq_records\t"
        f"{len(raw_records)}"
    )

    print(
        "diagnostic_rows\t"
        f"{len(diagnostic_rows)}"
    )

    print(
        "truth_rows\t"
        f"{len(truth_rows)}"
    )

    print(
        "prepared_fastq_records\t"
        f"{len(prepared_records)}"
    )

    print(
        "prepared_fastq_files\t"
        f"{len(per_file_counts)}"
    )

    # Print one count per prepared FASTQ.
    for filename in sorted(per_file_counts):
        print(
            f"prepared_file_{filename}\t"
            f"{per_file_counts[filename]}"
        )

    # ------------------------------------------------------------------
    # EXACT READ-ID SET COMPARISONS
    # ------------------------------------------------------------------

    # Build sets once so later comparisons are easy to read.
    raw_ids = set(raw_records)
    diagnostic_ids = set(diagnostic_rows)
    truth_ids = set(truth_rows)
    prepared_ids = set(prepared_records)

    # Compare diagnostic IDs against raw FASTQ IDs in both directions.
    raw_missing_from_diagnostic = (
        raw_ids - diagnostic_ids
    )

    diagnostic_missing_from_raw = (
        diagnostic_ids - raw_ids
    )

    # Compare production truth rows against raw IDs in both directions.
    raw_missing_from_truth = (
        raw_ids - truth_ids
    )

    truth_missing_from_raw = (
        truth_ids - raw_ids
    )

    # Compare actual prepared FASTQ IDs against raw IDs.
    raw_missing_from_prepared = (
        raw_ids - prepared_ids
    )

    prepared_missing_from_raw = (
        prepared_ids - raw_ids
    )

    # Report every set difference.
    report_set_difference(
        "raw_missing_from_diagnostic",
        raw_missing_from_diagnostic,
    )

    report_set_difference(
        "diagnostic_missing_from_raw",
        diagnostic_missing_from_raw,
    )

    report_set_difference(
        "raw_missing_from_truth",
        raw_missing_from_truth,
    )

    report_set_difference(
        "truth_missing_from_raw",
        truth_missing_from_raw,
    )

    report_set_difference(
        "raw_missing_from_prepared",
        raw_missing_from_prepared,
    )

    report_set_difference(
        "prepared_missing_from_raw",
        prepared_missing_from_raw,
    )

    # ------------------------------------------------------------------
    # FASTQ-FILE SET VALIDATION
    # ------------------------------------------------------------------

    # Files declared by truth_reads.tsv.
    declared_output_files = {
        row["output_fastq"]
        for row in truth_rows.values()
    }

    # Files actually present and containing records.
    actual_output_files = set(
        per_file_counts
    )

    # A declared file must physically exist.
    declared_files_missing = (
        declared_output_files
        - actual_output_files
    )

    # A physical prepared FASTQ should also be declared in provenance.
    undeclared_actual_files = (
        actual_output_files
        - declared_output_files
    )

    # Report those differences.
    report_set_difference(
        "declared_files_missing",
        declared_files_missing,
    )

    report_set_difference(
        "undeclared_actual_files",
        undeclared_actual_files,
    )

    # ------------------------------------------------------------------
    # READ-LEVEL VALIDATION COUNTERS
    # ------------------------------------------------------------------

    # Counter stores every possible mismatch class.
    errors: Counter = Counter()

    # These counters provide useful summaries even when validation passes.
    orientation_counts: Counter = Counter()
    boundary_status_counts: Counter = Counter()
    sample_counts: Counter = Counter()

    # Only reads present in all four data sources can undergo complete
    # read-level reconstruction.
    common_ids = (
        raw_ids
        & diagnostic_ids
        & truth_ids
        & prepared_ids
    )

    # Validate each common read deterministically by sorted ID.
    for read_id in sorted(common_ids):

        # Retrieve each representation of this read.
        raw_record = raw_records[read_id]
        diagnostic = diagnostic_rows[read_id]
        truth_row = truth_rows[read_id]
        prepared_record = prepared_records[read_id]

        # Record useful output summaries.
        orientation_counts[
            truth_row["input_orientation"]
        ] += 1

        boundary_status_counts[
            truth_row["boundary_status"]
        ] += 1

        sample_counts[
            truth_row["sample_id"]
        ] += 1

        # --------------------------------------------------------------
        # EXPERIMENT METADATA
        # --------------------------------------------------------------

        # Marker must be the one requested for this run.
        if (
            truth_row["marker"]
            != args.expected_marker
        ):
            errors["marker_mismatches"] += 1

        # Dataset ID must identify the intended simulated dataset.
        if (
            truth_row["dataset_id"]
            != args.expected_dataset_id
        ):
            errors["dataset_id_mismatches"] += 1

        # Replicate must be R01 for the current validation.
        if (
            truth_row["simulation_replicate"]
            != args.expected_simulation_replicate
        ):
            errors[
                "simulation_replicate_mismatches"
            ] += 1

        # Seed must correspond to this simulation replicate.
        if (
            truth_row["master_seed"]
            != str(args.expected_master_seed)
        ):
            errors["master_seed_mismatches"] += 1

        # --------------------------------------------------------------
        # ORIENTATION
        # --------------------------------------------------------------

        # Production orientation must exactly reproduce v0.2.
        production_orientation = truth_row[
            "input_orientation"
        ]

        diagnostic_orientation = diagnostic[
            "selected_orientation"
        ]

        if (
            production_orientation
            != diagnostic_orientation
        ):
            errors["orientation_mismatches"] += 1

        # --------------------------------------------------------------
        # BOUNDARY COORDINATES
        # --------------------------------------------------------------

        # v0.2 minimum target start is the production start policy.
        expected_start = int(
            diagnostic[
                "alignment_target_start_min"
            ]
        )

        # v0.2 maximum target end is the production end policy.
        expected_end = int(
            diagnostic[
                "alignment_target_end_max"
            ]
        )

        # Coordinates actually recorded by production.
        production_start = int(
            truth_row["target_start"]
        )

        production_end = int(
            truth_row["target_end"]
        )

        # Check start independently.
        if production_start != expected_start:
            errors["target_start_mismatches"] += 1

        # Check end independently.
        if production_end != expected_end:
            errors["target_end_mismatches"] += 1

        # --------------------------------------------------------------
        # BOUNDARY AMBIGUITY
        # --------------------------------------------------------------

        # Convert diagnostic flags to Python booleans.
        diagnostic_left_ambiguous = parse_boolean(
            diagnostic[
                "left_boundary_ambiguous"
            ]
        )

        diagnostic_right_ambiguous = parse_boolean(
            diagnostic[
                "right_boundary_ambiguous"
            ]
        )

        # Convert production flags independently.
        production_left_ambiguous = parse_boolean(
            truth_row[
                "left_boundary_ambiguous"
            ]
        )

        production_right_ambiguous = parse_boolean(
            truth_row[
                "right_boundary_ambiguous"
            ]
        )

        # Compare left flag.
        if (
            production_left_ambiguous
            != diagnostic_left_ambiguous
        ):
            errors[
                "left_boundary_flag_mismatches"
            ] += 1

        # Compare right flag.
        if (
            production_right_ambiguous
            != diagnostic_right_ambiguous
        ):
            errors[
                "right_boundary_flag_mismatches"
            ] += 1

        # Derive the expected policy label independently.
        if (
            diagnostic_left_ambiguous
            or diagnostic_right_ambiguous
        ):
            expected_boundary_status = (
                "ambiguous_retained"
            )
        else:
            expected_boundary_status = "ok"

        # Production must record the policy correctly.
        if (
            truth_row["boundary_status"]
            != expected_boundary_status
        ):
            errors[
                "boundary_status_mismatches"
            ] += 1

        # --------------------------------------------------------------
        # RAW FASTQ LENGTH AND QUALITY STRUCTURE
        # --------------------------------------------------------------

        # Recover the raw observed sequence.
        raw_sequence = str(
            raw_record.seq
        ).upper()

        # Recover its original numeric Phred qualities.
        raw_qualities = list(
            raw_record.letter_annotations[
                "phred_quality"
            ]
        )

        # Raw FASTQ must have exactly one quality value per base.
        if len(raw_sequence) != len(raw_qualities):
            errors[
                "raw_sequence_quality_length_mismatches"
            ] += 1

        # Production metadata must record the true raw length.
        if (
            int(truth_row["raw_length"])
            != len(raw_sequence)
        ):
            errors[
                "truth_raw_length_mismatches"
            ] += 1

        # Diagnostic observed_length must independently agree too.
        if (
            int(diagnostic["observed_length"])
            != len(raw_sequence)
        ):
            errors[
                "diagnostic_raw_length_mismatches"
            ] += 1

        # --------------------------------------------------------------
        # INDEPENDENT ORIENTATION RECONSTRUCTION
        # --------------------------------------------------------------

        # IMPORTANT:
        #
        # We do not call mhass_fastq_utils here.
        #
        # Instead we independently reproduce orientation with standard
        # Biopython + Python operations.
        if diagnostic_orientation == "as_read":

            # Forward-oriented reads keep their sequence unchanged.
            canonical_sequence = raw_sequence

            # Their quality order also remains unchanged.
            canonical_qualities = raw_qualities

        elif diagnostic_orientation == "revcomp":

            # Reverse-complement the raw sequence.
            canonical_sequence = str(
                Seq(
                    raw_sequence
                ).reverse_complement()
            )

            # Reverse, but do NOT complement, quality scores.
            canonical_qualities = list(
                reversed(raw_qualities)
            )

        else:
            # v0.2 should never contain any other orientation category.
            raise RuntimeError(
                "Unexpected diagnostic orientation for "
                f"{read_id}: "
                f"{diagnostic_orientation}"
            )

        # --------------------------------------------------------------
        # INDEPENDENT PREPARED READ RECONSTRUCTION
        # --------------------------------------------------------------

        # Coordinates must describe a valid half-open interval.
        if (
            expected_start < 0
            or expected_end < expected_start
            or expected_end > len(canonical_sequence)
        ):
            errors["invalid_diagnostic_coordinates"] += 1

            # We cannot safely slice this read, so skip sequence/quality
            # reconstruction for it.
            continue

        # Independently slice the canonical observed sequence.
        expected_sequence = canonical_sequence[
            expected_start:expected_end
        ]

        # Independently slice matching canonical quality scores.
        expected_qualities = canonical_qualities[
            expected_start:expected_end
        ]

        # Recover the actual production sequence.
        actual_sequence = str(
            prepared_record.seq
        ).upper()

        # Recover the actual production quality scores.
        actual_qualities = list(
            prepared_record.letter_annotations[
                "phred_quality"
            ]
        )

        # Prepared FASTQ itself must have synchronized bases/qualities.
        if (
            len(actual_sequence)
            != len(actual_qualities)
        ):
            errors[
                "prepared_sequence_quality_length_mismatches"
            ] += 1

        # Compare every prepared sequence base.
        if actual_sequence != expected_sequence:
            errors["sequence_mismatches"] += 1

        # Compare every Phred quality value, position by position.
        if actual_qualities != expected_qualities:
            errors["quality_mismatches"] += 1

        # --------------------------------------------------------------
        # PREPARED LENGTH
        # --------------------------------------------------------------

        # The coordinate interval itself predicts output length.
        expected_prepared_length = (
            expected_end - expected_start
        )

        # truth_reads.tsv must record that same length.
        if (
            int(truth_row["prepared_length"])
            != expected_prepared_length
        ):
            errors[
                "truth_prepared_length_mismatches"
            ] += 1

        # Actual FASTQ length must also equal the predicted length.
        if (
            len(actual_sequence)
            != expected_prepared_length
        ):
            errors[
                "actual_prepared_length_mismatches"
            ] += 1

        # --------------------------------------------------------------
        # SAMPLE / OUTPUT-FILE ASSIGNMENT
        # --------------------------------------------------------------

        # Example:
        #
        # sample_id = S02
        #
        # therefore the expected output filename is:
        #
        # S02.fastq
        expected_output_filename = (
            f"{truth_row['sample_id']}.fastq"
        )

        # truth_reads.tsv must name the expected file.
        if (
            truth_row["output_fastq"]
            != expected_output_filename
        ):
            errors[
                "truth_output_filename_mismatches"
            ] += 1

        # The read must physically occur in that same file.
        actual_output_filename = (
            prepared_locations[read_id]
        )

        if (
            actual_output_filename
            != expected_output_filename
        ):
            errors[
                "physical_output_file_mismatches"
            ] += 1

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    # Number of reads that could be completely compared.
    print(
        "reads_fully_compared\t"
        f"{len(common_ids)}"
    )

    # Print useful production summaries.
    for orientation in sorted(orientation_counts):
        print(
            f"orientation_{orientation}\t"
            f"{orientation_counts[orientation]}"
        )

    for status in sorted(boundary_status_counts):
        print(
            f"boundary_status_{status}\t"
            f"{boundary_status_counts[status]}"
        )

    for sample in sorted(sample_counts):
        print(
            f"sample_{sample}\t"
            f"{sample_counts[sample]}"
        )

    # Define every expected mismatch category explicitly.
    #
    # Printing zeroes is important because absence of terminal output
    # should not be confused with a test that was never performed.
    error_categories = [
        "marker_mismatches",
        "dataset_id_mismatches",
        "simulation_replicate_mismatches",
        "master_seed_mismatches",
        "orientation_mismatches",
        "target_start_mismatches",
        "target_end_mismatches",
        "left_boundary_flag_mismatches",
        "right_boundary_flag_mismatches",
        "boundary_status_mismatches",
        "raw_sequence_quality_length_mismatches",
        "truth_raw_length_mismatches",
        "diagnostic_raw_length_mismatches",
        "invalid_diagnostic_coordinates",
        "prepared_sequence_quality_length_mismatches",
        "sequence_mismatches",
        "quality_mismatches",
        "truth_prepared_length_mismatches",
        "actual_prepared_length_mismatches",
        "truth_output_filename_mismatches",
        "physical_output_file_mismatches",
    ]

    # Print every error count, including zeros.
    for category in error_categories:
        print(
            f"{category}\t"
            f"{errors[category]}"
        )

    # ------------------------------------------------------------------
    # FINAL PASS / FAIL DECISION
    # ------------------------------------------------------------------

    # Count all read-ID/file-set discrepancies.
    set_failures = (
        len(raw_missing_from_diagnostic)
        + len(diagnostic_missing_from_raw)
        + len(raw_missing_from_truth)
        + len(truth_missing_from_raw)
        + len(raw_missing_from_prepared)
        + len(prepared_missing_from_raw)
        + len(declared_files_missing)
        + len(undeclared_actual_files)
    )

    # Count all read-level errors.
    read_level_failures = sum(
        errors.values()
    )

    # Total validation failures.
    total_failures = (
        set_failures
        + read_level_failures
    )

    # Print the final failure count.
    print(
        "total_validation_failures\t"
        f"{total_failures}"
    )

    # Any disagreement results in a non-zero program exit status.
    if total_failures > 0:
        print("validation_status\tFAIL")
        return 1

    # Otherwise every dataset-level and read-level check passed.
    print("validation_status\tPASS")
    return 0


# Run main() only when the script is executed directly.
#
# Importing this file from a future unit test will therefore not
# automatically launch validation.
if __name__ == "__main__":
    raise SystemExit(main())
