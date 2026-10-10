#!/usr/bin/env python3

"""
Audit the complete M1 controlled MHASS test simulation.

PURPOSE
-------
This script validates the complete logical chain:

    frozen biological truth
        ->
    frozen intended counts
        ->
    MHASS intended-molecule mapping
        ->
    realized raw CCS reads
        ->
    realized-count QC tables
        ->
    realized-read truth table
        ->
    prepared per-sample FASTQs

This is intentionally a SIMULATION-LEVEL audit.

It does NOT repeat the detailed sequence, quality-score, orientation,
or boundary reconstruction already performed by
validate_prepared_dataset.py.

Instead, it asks whether all representations of the experiment agree
on:

- templates,
- samples,
- intended counts,
- mapping provenance,
- realized counts,
- realized read identities,
- realized read provenance,
- prepared sample assignment,
- marker/dataset/replicate/seed metadata.

A successful run must finish with:

    total_audit_failures    0
    simulation_audit_status PASS
"""

# Allow modern Python type annotations consistently.
from __future__ import annotations

# argparse parses explicit command-line inputs.
import argparse

# csv reads our tab-separated metadata and QC files.
import csv

# Counter efficiently counts template/sample combinations.
# defaultdict gives us convenient containers for error examples.
from collections import Counter, defaultdict

# Path provides safer filesystem-path handling.
from pathlib import Path

# SeqIO parses FASTA and FASTQ files correctly.
from Bio import SeqIO

# Seq is used for one small barcode-integrity check:
# ReverseBarcode -> reverse complement -> RevCompReverse.
from Bio.Seq import Seq


def parse_args() -> argparse.Namespace:
    """
    Define every file and expected design value required by the audit.

    Expected values are passed explicitly rather than hidden in the
    script so the audit command itself documents the dataset being
    validated.
    """

    # Create the command-line parser.
    parser = argparse.ArgumentParser(
        description=(
            "Audit a complete controlled MHASS simulation from frozen "
            "design through realized and prepared reads."
        )
    )

    # Frozen biological truth FASTA.
    parser.add_argument(
        "--truth",
        required=True,
        type=Path,
        help="Frozen truth_sequences.fasta.",
    )

    # Frozen intended template-by-sample counts.
    parser.add_argument(
        "--fixed-counts",
        required=True,
        type=Path,
        help="Frozen fixed_counts.tsv.",
    )

    # MHASS intended-molecule mapping.
    parser.add_argument(
        "--mapping",
        required=True,
        type=Path,
        help="MHASS sequence_file_mapping.tsv.",
    )

    # MHASS sample/barcode metadata.
    parser.add_argument(
        "--barcodes",
        required=True,
        type=Path,
        help="MHASS sample_barcode_map.tsv.",
    )

    # Raw realized CCS reads.
    parser.add_argument(
        "--raw-fastq",
        required=True,
        type=Path,
        help="Raw MHASS combined_reads.fastq.",
    )

    # QC table containing realized counts by template and sample.
    parser.add_argument(
        "--realized-template-counts",
        required=True,
        type=Path,
        help="realized_template_counts.tsv.",
    )

    # QC table comparing intended and realized counts for every
    # template/sample combination.
    parser.add_argument(
        "--intended-realized",
        required=True,
        type=Path,
        help="intended_vs_realized_template_sample.tsv.",
    )

    # Production realized-read truth/provenance table.
    parser.add_argument(
        "--truth-reads",
        required=True,
        type=Path,
        help="Prepared dataset truth_reads.tsv.",
    )

    # Directory containing S01.fastq, S02.fastq, and S03.fastq.
    parser.add_argument(
        "--prepared-dir",
        required=True,
        type=Path,
        help="Directory containing prepared per-sample FASTQs.",
    )

    # Explicit expected sample names.
    #
    # nargs="+" means the command can supply:
    #
    #   --expected-samples S01 S02 S03
    parser.add_argument(
        "--expected-samples",
        required=True,
        nargs="+",
        help="Expected sample IDs.",
    )

    # Explicit number of biological templates expected in this design.
    parser.add_argument(
        "--expected-template-count",
        required=True,
        type=int,
        help="Expected number of truth templates.",
    )

    # Explicit expected intended total across all samples.
    parser.add_argument(
        "--expected-intended-total",
        required=True,
        type=int,
        help="Expected total intended molecules.",
    )

    # Every sample in this design was configured to 5120 intended reads.
    parser.add_argument(
        "--expected-intended-per-sample",
        required=True,
        type=int,
        help="Expected intended count for each sample.",
    )

    # Explicit expected total realized CCS reads.
    parser.add_argument(
        "--expected-realized-total",
        required=True,
        type=int,
        help="Expected realized raw read count.",
    )

    # Expected realized-read metadata values.
    parser.add_argument(
        "--expected-marker",
        required=True,
        help="Expected marker, for example 16S.",
    )

    parser.add_argument(
        "--expected-dataset-id",
        required=True,
        help="Expected dataset identifier.",
    )

    parser.add_argument(
        "--expected-simulation-replicate",
        required=True,
        help="Expected simulation replicate.",
    )

    parser.add_argument(
        "--expected-master-seed",
        required=True,
        type=int,
        help="Expected MHASS master seed.",
    )

    # Parse the command line and return the resulting namespace.
    return parser.parse_args()


def add_error(
    errors: Counter,
    examples: dict[str, list[str]],
    category: str,
    detail: str,
    max_examples: int = 5,
) -> None:
    """
    Record one validation failure.

    We count every failure but store only a few example messages so a
    broken dataset remains diagnosable without flooding the terminal.
    """

    # Increment this category's failure count.
    errors[category] += 1

    # Preserve only the first few examples.
    if len(examples[category]) < max_examples:
        examples[category].append(detail)


def require_columns(
    reader: csv.DictReader,
    required: set[str],
    path: Path,
) -> None:
    """
    Fail immediately if a TSV is missing required columns.

    Schema failure is different from a scientific mismatch: if the file
    cannot even be interpreted safely, continuing would be misleading.
    """

    # Convert the actual header to a set.
    observed = set(reader.fieldnames or [])

    # Find any required fields that are absent.
    missing = required - observed

    # Stop if the schema is incomplete.
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {sorted(missing)}. "
            f"Observed columns: {reader.fieldnames}"
        )


def parse_nonnegative_int(
    value: str,
    context: str,
) -> int:
    """
    Parse a TSV value as a non-negative integer.

    Counts such as intended_reads and realized_reads must never be
    negative or fractional.
    """

    # Convert textual value to integer.
    number = int(value)

    # Negative counts are invalid.
    if number < 0:
        raise ValueError(
            f"Negative count in {context}: {number}"
        )

    # Return the validated integer.
    return number


def load_truth_ids(
    path: Path,
) -> set[str]:
    """
    Load truth FASTA IDs and reject duplicates.

    The truth FASTA defines the allowed biological template universe.
    """

    # Store IDs that have already appeared.
    truth_ids: set[str] = set()

    # Read every FASTA record.
    for record in SeqIO.parse(path, "fasta"):

        # Duplicate truth identifiers make provenance ambiguous.
        if record.id in truth_ids:
            raise ValueError(
                f"Duplicate truth FASTA ID: {record.id}"
            )

        # Record the unique template ID.
        truth_ids.add(record.id)

    # Empty truth is never valid.
    if not truth_ids:
        raise ValueError(
            "Truth FASTA contains no records."
        )

    # Return the complete template-ID set.
    return truth_ids


def load_fixed_counts(
    path: Path,
    expected_samples: list[str],
) -> tuple[
    dict[tuple[str, str], int],
    set[str],
    Counter,
]:
    """
    Load frozen intended counts.

    Actual schema:

        ASVID    S01    S02    S03

    Returns:

        pair_counts:
            (template_id, sample_id) -> intended count

        template_ids:
            all ASVID values present

        sample_totals:
            sample_id -> total intended reads
    """

    # Store one intended count for each template/sample combination.
    pair_counts: dict[tuple[str, str], int] = {}

    # Track templates present in the design.
    template_ids: set[str] = set()

    # Accumulate total intended molecules per sample.
    sample_totals: Counter = Counter()

    # Open the frozen TSV.
    with path.open(newline="") as handle:

        # Read by column name.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Require ASVID plus exactly the expected sample fields.
        require_columns(
            reader,
            {"ASVID", *expected_samples},
            path,
        )

        # Reject unexpected extra sample/design columns.
        actual_columns = list(reader.fieldnames or [])

        expected_columns = [
            "ASVID",
            *expected_samples,
        ]

        if actual_columns != expected_columns:
            raise ValueError(
                "fixed_counts.tsv columns/order differ from expected: "
                f"expected={expected_columns}, observed={actual_columns}"
            )

        # Process one template row at a time.
        for row in reader:

            # Frozen truth template identifier.
            template_id = row["ASVID"]

            # Duplicate rows would give two definitions of intended counts.
            if template_id in template_ids:
                raise ValueError(
                    f"Duplicate ASVID in fixed counts: {template_id}"
                )

            # Record the template.
            template_ids.add(template_id)

            # Read each sample's intended count.
            for sample_id in expected_samples:

                # Parse and validate the count.
                count = parse_nonnegative_int(
                    row[sample_id],
                    (
                        f"{path}: template={template_id}, "
                        f"sample={sample_id}"
                    ),
                )

                # Store exact template/sample intent.
                pair_counts[
                    (template_id, sample_id)
                ] = count

                # Add to the sample total.
                sample_totals[sample_id] += count

    # Return all frozen design structures.
    return pair_counts, template_ids, sample_totals


def load_barcodes(
    path: Path,
) -> dict[str, dict[str, str]]:
    """
    Load and validate the MHASS sample barcode table.

    Actual schema:

        SampleID
        BarcodeID
        ForwardBarcode
        ReverseBarcode
        RevCompReverse
    """

    # Sample -> barcode metadata.
    barcodes: dict[str, dict[str, str]] = {}

    # Also track BarcodeIDs to ensure they are unique between samples.
    used_barcode_ids: set[str] = set()

    # Open the table.
    with path.open(newline="") as handle:

        # Read by header names.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Validate the exact fields needed.
        require_columns(
            reader,
            {
                "SampleID",
                "BarcodeID",
                "ForwardBarcode",
                "ReverseBarcode",
                "RevCompReverse",
            },
            path,
        )

        # Process each sample.
        for row in reader:

            # Extract fields.
            sample_id = row["SampleID"]
            barcode_id = row["BarcodeID"]

            forward = row["ForwardBarcode"].upper()
            reverse = row["ReverseBarcode"].upper()
            revcomp_reverse = row["RevCompReverse"].upper()

            # Sample must appear once.
            if sample_id in barcodes:
                raise ValueError(
                    f"Duplicate SampleID in barcode table: {sample_id}"
                )

            # Barcode ID must also uniquely identify one sample.
            if barcode_id in used_barcode_ids:
                raise ValueError(
                    f"Duplicate BarcodeID in barcode table: {barcode_id}"
                )

            # Record this BarcodeID as used.
            used_barcode_ids.add(barcode_id)

            # Ensure barcode sequences are non-empty.
            if not forward or not reverse or not revcomp_reverse:
                raise ValueError(
                    f"Empty barcode sequence for sample {sample_id}"
                )

            # Independently confirm that RevCompReverse really is the
            # reverse complement of ReverseBarcode.
            calculated_revcomp = str(
                Seq(reverse).reverse_complement()
            )

            if calculated_revcomp != revcomp_reverse:
                raise ValueError(
                    "RevCompReverse mismatch for sample "
                    f"{sample_id}: calculated={calculated_revcomp}, "
                    f"table={revcomp_reverse}"
                )

            # Store normalized metadata.
            barcodes[sample_id] = {
                "barcode_id": barcode_id,
                "forward": forward,
                "reverse": reverse,
                "revcomp_reverse": revcomp_reverse,
            }

    # Return complete sample barcode metadata.
    return barcodes


def load_mapping(
    path: Path,
    errors: Counter,
    examples: dict[str, list[str]],
) -> tuple[
    Counter,
    dict[str, tuple[str, str, str, str]],
    int,
]:
    """
    Load MHASS sequence_file_mapping.tsv.

    Actual schema:

        ASVID
        SampleID
        BarcodeID
        TemplateFile

    Returns:

        pair_counts:
            number of intended mapping rows per template/sample

        group_provenance:
            FASTQ template-group prefix ->
            (ASVID, SampleID, BarcodeID, TemplateFile)

        row_count:
            total intended mapping rows
    """

    # Count mapping rows by template/sample.
    pair_counts: Counter = Counter()

    # Map each grouped TemplateFile to one consistent provenance tuple.
    group_provenance: dict[
        str,
        tuple[str, str, str, str],
    ] = {}

    # Count every intended mapping row.
    row_count = 0

    # Open the mapping file.
    with path.open(newline="") as handle:

        # Read by header name.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Validate schema.
        require_columns(
            reader,
            {
                "ASVID",
                "SampleID",
                "BarcodeID",
                "TemplateFile",
            },
            path,
        )

        # Process all 15,360 intended rows.
        for row in reader:

            # Increment total row count.
            row_count += 1

            # Extract provenance.
            template_id = row["ASVID"]
            sample_id = row["SampleID"]
            barcode_id = row["BarcodeID"]
            template_file = row["TemplateFile"]

            # Count this intended template/sample combination.
            pair_counts[
                (template_id, sample_id)
            ] += 1

            # FASTQ IDs use the TemplateFile basename without ".fasta".
            template_group = template_file

            if template_group.endswith(".fasta"):
                template_group = template_group[:-6]

            # This provenance must be identical every time the grouped
            # TemplateFile appears in the intended-molecule mapping.
            provenance = (
                template_id,
                sample_id,
                barcode_id,
                template_file,
            )

            # Detect conflicts before replacing an existing mapping.
            if (
                template_group in group_provenance
                and group_provenance[template_group] != provenance
            ):
                add_error(
                    errors,
                    examples,
                    "mapping_group_provenance_conflicts",
                    (
                        f"{template_group}: "
                        f"{group_provenance[template_group]} "
                        f"vs {provenance}"
                    ),
                )

            else:
                # Store or reaffirm consistent provenance.
                group_provenance[
                    template_group
                ] = provenance

    # Return all intended-mapping information.
    return pair_counts, group_provenance, row_count


def derive_raw_realized_counts(
    path: Path,
    group_provenance: dict[
        str,
        tuple[str, str, str, str],
    ],
    errors: Counter,
    examples: dict[str, list[str]],
) -> tuple[
    set[str],
    Counter,
    Counter,
]:
    """
    Parse raw realized FASTQ and independently derive realized counts.

    Returns:

        raw_ids:
            all unique realized read IDs

        pair_counts:
            realized counts by (template_id, sample_id)

        sample_counts:
            realized counts by sample
    """

    # Set guarantees uniqueness of realized read IDs.
    raw_ids: set[str] = set()

    # Independently derived realized template/sample counts.
    pair_counts: Counter = Counter()

    # Independently derived realized sample totals.
    sample_counts: Counter = Counter()

    # Parse each realized FASTQ read.
    for record in SeqIO.parse(path, "fastq"):

        # Read ID is the primary realized molecule identifier.
        read_id = record.id

        # Duplicated raw IDs would invalidate 1:1 provenance.
        if read_id in raw_ids:
            add_error(
                errors,
                examples,
                "duplicate_raw_read_ids",
                read_id,
            )

            # Do not count the duplicate twice.
            continue

        # Record the unique realized ID.
        raw_ids.add(read_id)

        # Example:
        #
        # template10_A1_np17/S/3/ccs
        #
        # becomes:
        #
        # template10_A1_np17
        template_group = read_id.split("/", 1)[0]

        # Every realized read must resolve to intended provenance.
        if template_group not in group_provenance:
            add_error(
                errors,
                examples,
                "raw_reads_missing_mapping_provenance",
                read_id,
            )

            # Cannot derive template/sample safely for this read.
            continue

        # Retrieve the mapping provenance.
        (
            template_id,
            sample_id,
            _barcode_id,
            _template_file,
        ) = group_provenance[template_group]

        # Count this realized template/sample combination.
        pair_counts[
            (template_id, sample_id)
        ] += 1

        # Count this realized sample.
        sample_counts[sample_id] += 1

    # Return independently derived realized truth.
    return raw_ids, pair_counts, sample_counts


def load_realized_template_counts(
    path: Path,
    expected_samples: list[str],
) -> tuple[
    dict[tuple[str, str], int],
    set[str],
]:
    """
    Load realized_template_counts.tsv.

    Actual schema:

        template_id    S01    S02    S03
    """

    # Table values keyed by template/sample.
    pair_counts: dict[tuple[str, str], int] = {}

    # Templates represented by table rows.
    template_ids: set[str] = set()

    # Open the QC table.
    with path.open(newline="") as handle:

        # Read by header names.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Require template_id plus the three expected sample columns.
        require_columns(
            reader,
            {"template_id", *expected_samples},
            path,
        )

        # Require exactly the expected column order.
        expected_columns = [
            "template_id",
            *expected_samples,
        ]

        actual_columns = list(
            reader.fieldnames or []
        )

        if actual_columns != expected_columns:
            raise ValueError(
                "realized_template_counts.tsv columns/order differ "
                f"from expected: expected={expected_columns}, "
                f"observed={actual_columns}"
            )

        # Process every template row.
        for row in reader:

            # Retrieve template ID.
            template_id = row["template_id"]

            # Reject duplicate template rows.
            if template_id in template_ids:
                raise ValueError(
                    "Duplicate template_id in realized counts: "
                    f"{template_id}"
                )

            # Record template.
            template_ids.add(template_id)

            # Parse one count per sample.
            for sample_id in expected_samples:

                count = parse_nonnegative_int(
                    row[sample_id],
                    (
                        f"{path}: template={template_id}, "
                        f"sample={sample_id}"
                    ),
                )

                pair_counts[
                    (template_id, sample_id)
                ] = count

    # Return table values.
    return pair_counts, template_ids


def load_intended_realized(
    path: Path,
) -> dict[
    tuple[str, str],
    dict[str, str],
]:
    """
    Load intended_vs_realized_template_sample.tsv.

    Actual schema:

        template_id
        sample_id
        intended_reads
        realized_reads
        lost_reads
        realization_fraction
    """

    # One row per template/sample combination.
    rows: dict[
        tuple[str, str],
        dict[str, str],
    ] = {}

    # Open TSV.
    with path.open(newline="") as handle:

        # Read by header names.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Require all six fields.
        require_columns(
            reader,
            {
                "template_id",
                "sample_id",
                "intended_reads",
                "realized_reads",
                "lost_reads",
                "realization_fraction",
            },
            path,
        )

        # Process all rows.
        for row in reader:

            # Pair is the table's natural unique key.
            key = (
                row["template_id"],
                row["sample_id"],
            )

            # Every template/sample pair must appear only once.
            if key in rows:
                raise ValueError(
                    "Duplicate template/sample row in "
                    f"{path}: {key}"
                )

            # Preserve the complete textual row.
            rows[key] = row

    # Return indexed table.
    return rows


def fraction_matches(
    text: str,
    intended: int,
    realized: int,
    tolerance: float = 1e-6,
) -> bool:
    """
    Check realization_fraction safely.

    For intended > 0:

        expected = realized / intended

    For intended == 0:

        realized must also be zero, and we accept either:
        - blank / NA / NaN / N/A / "."
        - numeric zero

    This avoids assuming one specific missing-value convention.
    """

    # Normalize surrounding whitespace.
    normalized = text.strip()

    # Handle intended-zero combinations explicitly.
    if intended == 0:

        # A zero-intended pair must never realize reads.
        if realized != 0:
            return False

        # Common textual missing-value representations are acceptable.
        if normalized.lower() in {
            "",
            "na",
            "nan",
            "n/a",
            ".",
        }:
            return True

        # Otherwise allow explicit numeric zero.
        try:
            return abs(float(normalized)) <= tolerance

        except ValueError:
            return False

    # For positive intended counts, calculate expected fraction.
    expected = realized / intended

    # realization_fraction must be numeric here.
    try:
        observed = float(normalized)

    except ValueError:
        return False

    # Compare with tolerance to allow rounded decimal output.
    return abs(observed - expected) <= tolerance


def load_truth_reads(
    path: Path,
) -> dict[str, dict[str, str]]:
    """
    Load production truth_reads.tsv keyed by read_id.
    """

    # Realized-read provenance rows.
    rows: dict[str, dict[str, str]] = {}

    # Open table.
    with path.open(newline="") as handle:

        # Read by named fields.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Require fields needed by this simulation-level audit.
        require_columns(
            reader,
            {
                "read_id",
                "marker",
                "dataset_id",
                "sample_id",
                "simulation_replicate",
                "master_seed",
                "truth_template_id",
                "barcode_id",
                "mhass_template_file",
                "output_fastq",
            },
            path,
        )

        # Read every realized-read truth row.
        for row in reader:

            # Primary key.
            read_id = row["read_id"]

            # Duplicate truth rows violate 1:1 provenance.
            if read_id in rows:
                raise ValueError(
                    f"Duplicate read_id in truth_reads.tsv: {read_id}"
                )

            # Store complete row.
            rows[read_id] = row

    # Return indexed truth table.
    return rows


def load_prepared_fastq_locations(
    prepared_dir: Path,
) -> tuple[
    set[str],
    dict[str, str],
    Counter,
]:
    """
    Load prepared FASTQ IDs without revalidating sequence contents.

    Returns:

        prepared_ids:
            all unique prepared read IDs

        locations:
            read_id -> physical FASTQ filename

        file_counts:
            FASTQ filename -> number of records
    """

    # Find all prepared FASTQs.
    fastq_paths = sorted(
        prepared_dir.glob("*.fastq")
    )

    # At least one FASTQ must exist.
    if not fastq_paths:
        raise ValueError(
            f"No prepared FASTQ files found in {prepared_dir}"
        )

    # Unique read IDs across all files.
    prepared_ids: set[str] = set()

    # Physical filename containing each read.
    locations: dict[str, str] = {}

    # Number of reads in each physical file.
    file_counts: Counter = Counter()

    # Parse each sample FASTQ.
    for path in fastq_paths:

        for record in SeqIO.parse(
            path,
            "fastq",
        ):

            # Duplicate IDs across prepared files are forbidden.
            if record.id in prepared_ids:
                raise ValueError(
                    "Duplicate prepared read ID across FASTQs: "
                    f"{record.id}"
                )

            # Register the unique read.
            prepared_ids.add(record.id)

            # Save the physical filename.
            locations[record.id] = path.name

            # Increment physical file count.
            file_counts[path.name] += 1

    # Return prepared-dataset structures.
    return prepared_ids, locations, file_counts


def compare_id_sets(
    left_name: str,
    left_ids: set[str],
    right_name: str,
    right_ids: set[str],
    errors: Counter,
    examples: dict[str, list[str]],
) -> None:
    """
    Compare two read-ID sets in both directions.

    Equal counts alone are insufficient: the identities themselves
    must be identical.
    """

    # IDs present on the left but absent on the right.
    left_only = left_ids - right_ids

    # IDs present on the right but absent on the left.
    right_only = right_ids - left_ids

    # Record every missing ID as a failure.
    for read_id in left_only:
        add_error(
            errors,
            examples,
            f"{left_name}_missing_from_{right_name}",
            read_id,
        )

    for read_id in right_only:
        add_error(
            errors,
            examples,
            f"{right_name}_missing_from_{left_name}",
            read_id,
        )


def main() -> int:
    """
    Execute the complete simulation-level audit.
    """

    # Read explicit command-line configuration.
    args = parse_args()

    # Normalize expected samples into a stable ordered list.
    expected_samples = list(
        args.expected_samples
    )

    # Duplicate sample names on the command line would be misleading.
    if len(expected_samples) != len(set(expected_samples)):
        raise ValueError(
            f"Duplicate expected sample IDs: {expected_samples}"
        )

    # Error counter.
    errors: Counter = Counter()

    # A few diagnostic examples for each error category.
    examples: dict[str, list[str]] = defaultdict(list)

    # ==============================================================
    # 1. FROZEN BIOLOGICAL TRUTH
    # ==============================================================

    # Load the allowed template universe.
    truth_ids = load_truth_ids(
        args.truth
    )

    # Validate expected template count.
    if len(truth_ids) != args.expected_template_count:
        add_error(
            errors,
            examples,
            "truth_template_count_mismatches",
            (
                f"observed={len(truth_ids)}, "
                f"expected={args.expected_template_count}"
            ),
        )

    # ==============================================================
    # 2. FROZEN INTENDED COUNTS
    # ==============================================================

    (
        fixed_pair_counts,
        fixed_template_ids,
        fixed_sample_totals,
    ) = load_fixed_counts(
        args.fixed_counts,
        expected_samples,
    )

    # Frozen counts should describe exactly the truth FASTA templates.
    for template_id in (
        truth_ids - fixed_template_ids
    ):
        add_error(
            errors,
            examples,
            "truth_templates_missing_from_fixed_counts",
            template_id,
        )

    for template_id in (
        fixed_template_ids - truth_ids
    ):
        add_error(
            errors,
            examples,
            "fixed_count_templates_missing_from_truth",
            template_id,
        )

    # Calculate complete intended molecule total.
    fixed_intended_total = sum(
        fixed_pair_counts.values()
    )

    # Validate known intended total: 15,360.
    if (
        fixed_intended_total
        != args.expected_intended_total
    ):
        add_error(
            errors,
            examples,
            "fixed_intended_total_mismatches",
            (
                f"observed={fixed_intended_total}, "
                f"expected={args.expected_intended_total}"
            ),
        )

    # Validate 5,120 intended molecules per sample.
    for sample_id in expected_samples:

        observed = fixed_sample_totals[
            sample_id
        ]

        if (
            observed
            != args.expected_intended_per_sample
        ):
            add_error(
                errors,
                examples,
                "fixed_sample_total_mismatches",
                (
                    f"{sample_id}: observed={observed}, "
                    f"expected={args.expected_intended_per_sample}"
                ),
            )

    # ==============================================================
    # 3. BARCODE TABLE
    # ==============================================================

    # Load and internally validate barcode metadata.
    barcodes = load_barcodes(
        args.barcodes
    )

    # Expected sample set must exactly match barcode-table sample set.
    expected_sample_set = set(
        expected_samples
    )

    barcode_sample_set = set(
        barcodes
    )

    for sample_id in (
        expected_sample_set
        - barcode_sample_set
    ):
        add_error(
            errors,
            examples,
            "expected_samples_missing_from_barcodes",
            sample_id,
        )

    for sample_id in (
        barcode_sample_set
        - expected_sample_set
    ):
        add_error(
            errors,
            examples,
            "unexpected_barcode_samples",
            sample_id,
        )

    # ==============================================================
    # 4. MHASS INTENDED-MOLECULE MAPPING
    # ==============================================================

    (
        mapping_pair_counts,
        group_provenance,
        mapping_row_count,
    ) = load_mapping(
        args.mapping,
        errors,
        examples,
    )

    # Total mapping rows must equal complete intended design total.
    if mapping_row_count != fixed_intended_total:
        add_error(
            errors,
            examples,
            "mapping_total_vs_fixed_total_mismatches",
            (
                f"mapping={mapping_row_count}, "
                f"fixed={fixed_intended_total}"
            ),
        )

    # Also compare against explicit expected 15,360 total.
    if (
        mapping_row_count
        != args.expected_intended_total
    ):
        add_error(
            errors,
            examples,
            "mapping_total_vs_expected_mismatches",
            (
                f"observed={mapping_row_count}, "
                f"expected={args.expected_intended_total}"
            ),
        )

    # Mapping counts must exactly reproduce every frozen pair count.
    all_design_pairs = set(
        fixed_pair_counts
    ) | set(
        mapping_pair_counts
    )

    for pair in all_design_pairs:

        fixed_count = fixed_pair_counts.get(
            pair,
            0,
        )

        mapping_count = mapping_pair_counts.get(
            pair,
            0,
        )

        if fixed_count != mapping_count:
            add_error(
                errors,
                examples,
                "mapping_pair_count_mismatches",
                (
                    f"{pair}: fixed={fixed_count}, "
                    f"mapping={mapping_count}"
                ),
            )

    # Every mapping provenance group must use valid template/sample/barcode.
    for (
        template_group,
        provenance,
    ) in group_provenance.items():

        (
            template_id,
            sample_id,
            barcode_id,
            _template_file,
        ) = provenance

        # Template must exist in biological truth.
        if template_id not in truth_ids:
            add_error(
                errors,
                examples,
                "mapping_templates_missing_from_truth",
                (
                    f"{template_group}: "
                    f"{template_id}"
                ),
            )

        # Sample must be expected.
        if sample_id not in expected_sample_set:
            add_error(
                errors,
                examples,
                "mapping_unexpected_samples",
                (
                    f"{template_group}: "
                    f"{sample_id}"
                ),
            )

            # Cannot safely check barcode if sample itself is unknown.
            continue

        # Sample must have barcode metadata.
        if sample_id not in barcodes:
            add_error(
                errors,
                examples,
                "mapping_samples_missing_barcode_metadata",
                (
                    f"{template_group}: "
                    f"{sample_id}"
                ),
            )

            continue

        # Mapping BarcodeID must equal the sample barcode-table BarcodeID.
        expected_barcode_id = barcodes[
            sample_id
        ]["barcode_id"]

        if barcode_id != expected_barcode_id:
            add_error(
                errors,
                examples,
                "mapping_barcode_id_mismatches",
                (
                    f"{template_group}: "
                    f"mapping={barcode_id}, "
                    f"expected={expected_barcode_id}"
                ),
            )

    # ==============================================================
    # 5. RAW REALIZED CCS READS
    # ==============================================================

    (
        raw_ids,
        raw_pair_counts,
        raw_sample_counts,
    ) = derive_raw_realized_counts(
        args.raw_fastq,
        group_provenance,
        errors,
        examples,
    )

    # Validate complete realized total: 10,228.
    if (
        len(raw_ids)
        != args.expected_realized_total
    ):
        add_error(
            errors,
            examples,
            "raw_realized_total_mismatches",
            (
                f"observed={len(raw_ids)}, "
                f"expected={args.expected_realized_total}"
            ),
        )

    # Realized reads may never exceed intended reads for a pair.
    for pair, realized_count in raw_pair_counts.items():

        intended_count = fixed_pair_counts.get(
            pair,
            0,
        )

        if realized_count > intended_count:
            add_error(
                errors,
                examples,
                "realized_exceeds_intended",
                (
                    f"{pair}: intended={intended_count}, "
                    f"realized={realized_count}"
                ),
            )

    # Explicitly audit intended-zero combinations.
    zero_intended_pairs = {
        pair
        for pair, count
        in fixed_pair_counts.items()
        if count == 0
    }

    for pair in zero_intended_pairs:

        if raw_pair_counts.get(pair, 0) != 0:
            add_error(
                errors,
                examples,
                "zero_intended_pairs_with_realized_reads",
                (
                    f"{pair}: realized="
                    f"{raw_pair_counts.get(pair, 0)}"
                ),
            )

    # ==============================================================
    # 6. REALIZED_TEMPLATE_COUNTS.TSV
    # ==============================================================

    (
        realized_table_counts,
        realized_table_templates,
    ) = load_realized_template_counts(
        args.realized_template_counts,
        expected_samples,
    )

    # Realized count table should represent exactly truth templates.
    for template_id in (
        truth_ids
        - realized_table_templates
    ):
        add_error(
            errors,
            examples,
            "truth_templates_missing_from_realized_table",
            template_id,
        )

    for template_id in (
        realized_table_templates
        - truth_ids
    ):
        add_error(
            errors,
            examples,
            "realized_table_templates_missing_from_truth",
            template_id,
        )

    # Compare every table cell against counts independently derived
    # from the raw FASTQ + mapping provenance.
    all_realized_pairs = set(
        fixed_pair_counts
    ) | set(
        realized_table_counts
    )

    for pair in all_realized_pairs:

        independently_derived = (
            raw_pair_counts.get(
                pair,
                0,
            )
        )

        table_value = (
            realized_table_counts.get(
                pair,
                0,
            )
        )

        if table_value != independently_derived:
            add_error(
                errors,
                examples,
                "realized_table_count_mismatches",
                (
                    f"{pair}: table={table_value}, "
                    f"derived={independently_derived}"
                ),
            )

    # ==============================================================
    # 7. INTENDED_VS_REALIZED TABLE
    # ==============================================================

    intended_realized_rows = load_intended_realized(
        args.intended_realized
    )

    # The expected pair universe is all 20 templates x 3 samples.
    expected_pairs = set(
        fixed_pair_counts
    )

    # Check for missing pair rows.
    for pair in (
        expected_pairs
        - set(intended_realized_rows)
    ):
        add_error(
            errors,
            examples,
            "pairs_missing_from_intended_realized_table",
            str(pair),
        )

    # Check for unexpected pair rows.
    for pair in (
        set(intended_realized_rows)
        - expected_pairs
    ):
        add_error(
            errors,
            examples,
            "unexpected_pairs_in_intended_realized_table",
            str(pair),
        )

    # Validate every expected pair row mathematically.
    for pair in expected_pairs:

        # Missing rows were already recorded above.
        if pair not in intended_realized_rows:
            continue

        # Retrieve textual table row.
        row = intended_realized_rows[
            pair
        ]

        # Frozen intended truth.
        expected_intended = fixed_pair_counts[
            pair
        ]

        # Independently derived realized truth.
        expected_realized = raw_pair_counts.get(
            pair,
            0,
        )

        # Parse table's intended field.
        observed_intended = parse_nonnegative_int(
            row["intended_reads"],
            f"{pair}: intended_reads",
        )

        # Parse table's realized field.
        observed_realized = parse_nonnegative_int(
            row["realized_reads"],
            f"{pair}: realized_reads",
        )

        # Parse table's lost count.
        observed_lost = parse_nonnegative_int(
            row["lost_reads"],
            f"{pair}: lost_reads",
        )

        # Expected lost count.
        expected_lost = (
            expected_intended
            - expected_realized
        )

        # Compare intended field.
        if observed_intended != expected_intended:
            add_error(
                errors,
                examples,
                "intended_realized_intended_mismatches",
                (
                    f"{pair}: table={observed_intended}, "
                    f"expected={expected_intended}"
                ),
            )

        # Compare realized field.
        if observed_realized != expected_realized:
            add_error(
                errors,
                examples,
                "intended_realized_realized_mismatches",
                (
                    f"{pair}: table={observed_realized}, "
                    f"expected={expected_realized}"
                ),
            )

        # Compare lost field.
        if observed_lost != expected_lost:
            add_error(
                errors,
                examples,
                "intended_realized_lost_mismatches",
                (
                    f"{pair}: table={observed_lost}, "
                    f"expected={expected_lost}"
                ),
            )

        # Check realization_fraction.
        if not fraction_matches(
            row["realization_fraction"],
            expected_intended,
            expected_realized,
        ):
            add_error(
                errors,
                examples,
                "realization_fraction_mismatches",
                (
                    f"{pair}: table="
                    f"{row['realization_fraction']!r}, "
                    f"expected="
                    f"{expected_realized}/{expected_intended}"
                ),
            )

    # ==============================================================
    # 8. PRODUCTION TRUTH_READS.TSV
    # ==============================================================

    truth_rows = load_truth_reads(
        args.truth_reads
    )

    # IDs must match raw realized read IDs exactly.
    truth_read_ids = set(
        truth_rows
    )

    compare_id_sets(
        "raw",
        raw_ids,
        "truth",
        truth_read_ids,
        errors,
        examples,
    )

    # Independently count truth rows by template/sample.
    truth_pair_counts: Counter = Counter()

    # Validate every realized-read truth row against mapping provenance.
    for read_id, row in truth_rows.items():

        # Recover FASTQ/mapping template group from read ID.
        template_group = read_id.split(
            "/",
            1,
        )[0]

        # Every truth row must trace back to intended mapping provenance.
        if template_group not in group_provenance:
            add_error(
                errors,
                examples,
                "truth_reads_missing_mapping_provenance",
                read_id,
            )

            continue

        # Retrieve authoritative intended provenance.
        (
            expected_template_id,
            expected_sample_id,
            expected_barcode_id,
            expected_template_file,
        ) = group_provenance[
            template_group
        ]

        # Compare biological template.
        if (
            row["truth_template_id"]
            != expected_template_id
        ):
            add_error(
                errors,
                examples,
                "truth_read_template_mismatches",
                (
                    f"{read_id}: table="
                    f"{row['truth_template_id']}, "
                    f"mapping={expected_template_id}"
                ),
            )

        # Compare sample.
        if (
            row["sample_id"]
            != expected_sample_id
        ):
            add_error(
                errors,
                examples,
                "truth_read_sample_mismatches",
                (
                    f"{read_id}: table="
                    f"{row['sample_id']}, "
                    f"mapping={expected_sample_id}"
                ),
            )

        # Compare barcode ID.
        if (
            row["barcode_id"]
            != expected_barcode_id
        ):
            add_error(
                errors,
                examples,
                "truth_read_barcode_mismatches",
                (
                    f"{read_id}: table="
                    f"{row['barcode_id']}, "
                    f"mapping={expected_barcode_id}"
                ),
            )

        # Compare MHASS grouped template filename.
        if (
            row["mhass_template_file"]
            != expected_template_file
        ):
            add_error(
                errors,
                examples,
                "truth_read_template_file_mismatches",
                (
                    f"{read_id}: table="
                    f"{row['mhass_template_file']}, "
                    f"mapping={expected_template_file}"
                ),
            )

        # Validate marker metadata.
        if row["marker"] != args.expected_marker:
            add_error(
                errors,
                examples,
                "truth_read_marker_mismatches",
                (
                    f"{read_id}: {row['marker']}"
                ),
            )

        # Validate dataset ID.
        if (
            row["dataset_id"]
            != args.expected_dataset_id
        ):
            add_error(
                errors,
                examples,
                "truth_read_dataset_id_mismatches",
                (
                    f"{read_id}: {row['dataset_id']}"
                ),
            )

        # Validate replicate.
        if (
            row["simulation_replicate"]
            != args.expected_simulation_replicate
        ):
            add_error(
                errors,
                examples,
                "truth_read_replicate_mismatches",
                (
                    f"{read_id}: "
                    f"{row['simulation_replicate']}"
                ),
            )

        # Validate master seed.
        if (
            row["master_seed"]
            != str(args.expected_master_seed)
        ):
            add_error(
                errors,
                examples,
                "truth_read_master_seed_mismatches",
                (
                    f"{read_id}: "
                    f"{row['master_seed']}"
                ),
            )

        # output_fastq must follow sample identity.
        expected_output_fastq = (
            f"{row['sample_id']}.fastq"
        )

        if (
            row["output_fastq"]
            != expected_output_fastq
        ):
            add_error(
                errors,
                examples,
                "truth_read_output_filename_mismatches",
                (
                    f"{read_id}: table="
                    f"{row['output_fastq']}, "
                    f"expected={expected_output_fastq}"
                ),
            )

        # Count using mapping-authoritative provenance, not potentially
        # incorrect values from the truth row itself.
        truth_pair_counts[
            (
                expected_template_id,
                expected_sample_id,
            )
        ] += 1

    # Truth rows must independently reproduce raw realized pair counts.
    for pair in expected_pairs:

        if (
            truth_pair_counts.get(pair, 0)
            != raw_pair_counts.get(pair, 0)
        ):
            add_error(
                errors,
                examples,
                "truth_read_pair_count_mismatches",
                (
                    f"{pair}: truth="
                    f"{truth_pair_counts.get(pair, 0)}, "
                    f"raw="
                    f"{raw_pair_counts.get(pair, 0)}"
                ),
            )

    # ==============================================================
    # 9. PREPARED FASTQS
    # ==============================================================

    (
        prepared_ids,
        prepared_locations,
        prepared_file_counts,
    ) = load_prepared_fastq_locations(
        args.prepared_dir
    )

    # Prepared read IDs must equal raw realized IDs exactly.
    compare_id_sets(
        "raw",
        raw_ids,
        "prepared",
        prepared_ids,
        errors,
        examples,
    )

    # Expected physical FASTQ filenames.
    expected_fastq_files = {
        f"{sample_id}.fastq"
        for sample_id in expected_samples
    }

    # Actual physical FASTQ filenames.
    actual_fastq_files = set(
        prepared_file_counts
    )

    # Detect missing expected files.
    for filename in (
        expected_fastq_files
        - actual_fastq_files
    ):
        add_error(
            errors,
            examples,
            "prepared_expected_fastq_files_missing",
            filename,
        )

    # Detect unexpected extra FASTQ files.
    for filename in (
        actual_fastq_files
        - expected_fastq_files
    ):
        add_error(
            errors,
            examples,
            "prepared_unexpected_fastq_files",
            filename,
        )

    # Verify every prepared read is physically in the sample file
    # assigned by truth_reads.tsv.
    for read_id in (
        prepared_ids
        & truth_read_ids
    ):

        # Truth row's biological sample.
        sample_id = truth_rows[
            read_id
        ]["sample_id"]

        # Expected physical file.
        expected_filename = (
            f"{sample_id}.fastq"
        )

        # Actual physical file.
        actual_filename = (
            prepared_locations[
                read_id
            ]
        )

        if actual_filename != expected_filename:
            add_error(
                errors,
                examples,
                "prepared_physical_sample_assignment_mismatches",
                (
                    f"{read_id}: actual={actual_filename}, "
                    f"expected={expected_filename}"
                ),
            )

    # ==============================================================
    # 10. CROSS-LAYER SAMPLE TOTALS
    # ==============================================================

    # Truth-read sample totals.
    truth_sample_counts: Counter = Counter(
        row["sample_id"]
        for row in truth_rows.values()
    )

    # Compare raw, truth, and prepared counts per sample.
    for sample_id in expected_samples:

        # Independently derived raw count.
        raw_count = raw_sample_counts[
            sample_id
        ]

        # Production truth-table count.
        truth_count = truth_sample_counts[
            sample_id
        ]

        # Physical prepared FASTQ count.
        prepared_count = (
            prepared_file_counts[
                f"{sample_id}.fastq"
            ]
        )

        # Raw vs truth.
        if raw_count != truth_count:
            add_error(
                errors,
                examples,
                "raw_vs_truth_sample_count_mismatches",
                (
                    f"{sample_id}: raw={raw_count}, "
                    f"truth={truth_count}"
                ),
            )

        # Raw vs prepared.
        if raw_count != prepared_count:
            add_error(
                errors,
                examples,
                "raw_vs_prepared_sample_count_mismatches",
                (
                    f"{sample_id}: raw={raw_count}, "
                    f"prepared={prepared_count}"
                ),
            )

    # ==============================================================
    # 11. REPORT HIGH-LEVEL AUDIT COUNTS
    # ==============================================================

    # Frozen biological template count.
    print(
        "truth_templates\t"
        f"{len(truth_ids)}"
    )

    # Number of expected samples.
    print(
        "expected_samples\t"
        f"{len(expected_samples)}"
    )

    # Explicit sample names.
    print(
        "sample_ids\t"
        + ",".join(expected_samples)
    )

    # Frozen intended total.
    print(
        "fixed_intended_total\t"
        f"{fixed_intended_total}"
    )

    # One intended total per sample.
    for sample_id in expected_samples:
        print(
            f"fixed_intended_{sample_id}\t"
            f"{fixed_sample_totals[sample_id]}"
        )

    # Mapping size.
    print(
        "mapping_rows\t"
        f"{mapping_row_count}"
    )

    # Number of unique grouped template FASTAs in MHASS mapping.
    print(
        "mapping_unique_template_groups\t"
        f"{len(group_provenance)}"
    )

    # Raw realized total.
    print(
        "raw_realized_reads\t"
        f"{len(raw_ids)}"
    )

    # One realized raw count per sample.
    for sample_id in expected_samples:
        print(
            f"raw_realized_{sample_id}\t"
            f"{raw_sample_counts[sample_id]}"
        )

    # Number of intended-zero template/sample cells.
    print(
        "zero_intended_pairs\t"
        f"{len(zero_intended_pairs)}"
    )

    # Realized QC table size.
    print(
        "realized_table_templates\t"
        f"{len(realized_table_templates)}"
    )

    # Intended-vs-realized table row count.
    print(
        "intended_realized_rows\t"
        f"{len(intended_realized_rows)}"
    )

    # Production truth-table size.
    print(
        "truth_read_rows\t"
        f"{len(truth_rows)}"
    )

    # Prepared physical read total.
    print(
        "prepared_reads\t"
        f"{len(prepared_ids)}"
    )

    # One physical FASTQ count per sample file.
    for filename in sorted(
        prepared_file_counts
    ):
        print(
            f"prepared_{filename}\t"
            f"{prepared_file_counts[filename]}"
        )

    # ==============================================================
    # 12. PRINT ALL AUDIT FAILURE COUNTS
    # ==============================================================

    # Define expected failure categories explicitly.
    #
    # Printing zeroes proves the check was performed rather than merely
    # producing no output.
    error_categories = [
        "truth_template_count_mismatches",
        "truth_templates_missing_from_fixed_counts",
        "fixed_count_templates_missing_from_truth",
        "fixed_intended_total_mismatches",
        "fixed_sample_total_mismatches",
        "expected_samples_missing_from_barcodes",
        "unexpected_barcode_samples",
        "mapping_group_provenance_conflicts",
        "mapping_total_vs_fixed_total_mismatches",
        "mapping_total_vs_expected_mismatches",
        "mapping_pair_count_mismatches",
        "mapping_templates_missing_from_truth",
        "mapping_unexpected_samples",
        "mapping_samples_missing_barcode_metadata",
        "mapping_barcode_id_mismatches",
        "duplicate_raw_read_ids",
        "raw_reads_missing_mapping_provenance",
        "raw_realized_total_mismatches",
        "realized_exceeds_intended",
        "zero_intended_pairs_with_realized_reads",
        "truth_templates_missing_from_realized_table",
        "realized_table_templates_missing_from_truth",
        "realized_table_count_mismatches",
        "pairs_missing_from_intended_realized_table",
        "unexpected_pairs_in_intended_realized_table",
        "intended_realized_intended_mismatches",
        "intended_realized_realized_mismatches",
        "intended_realized_lost_mismatches",
        "realization_fraction_mismatches",
        "raw_missing_from_truth",
        "truth_missing_from_raw",
        "truth_reads_missing_mapping_provenance",
        "truth_read_template_mismatches",
        "truth_read_sample_mismatches",
        "truth_read_barcode_mismatches",
        "truth_read_template_file_mismatches",
        "truth_read_marker_mismatches",
        "truth_read_dataset_id_mismatches",
        "truth_read_replicate_mismatches",
        "truth_read_master_seed_mismatches",
        "truth_read_output_filename_mismatches",
        "truth_read_pair_count_mismatches",
        "raw_missing_from_prepared",
        "prepared_missing_from_raw",
        "prepared_expected_fastq_files_missing",
        "prepared_unexpected_fastq_files",
        "prepared_physical_sample_assignment_mismatches",
        "raw_vs_truth_sample_count_mismatches",
        "raw_vs_prepared_sample_count_mismatches",
    ]

    # Print every failure category, including categories with zero errors.
    for category in error_categories:
        print(
            f"{category}\t"
            f"{errors[category]}"
        )

    # ==============================================================
    # 13. PRINT LIMITED FAILURE EXAMPLES
    # ==============================================================

    # Only failure categories with recorded examples are printed here.
    #
    # This section is usually empty for a successful audit.
    for category in sorted(examples):

        for detail in examples[category]:
            print(
                f"{category}_example\t"
                f"{detail}"
            )

    # ==============================================================
    # 14. FINAL PASS / FAIL
    # ==============================================================

    # Sum every audit discrepancy.
    total_failures = sum(
        errors.values()
    )

    # Print explicit aggregate count.
    print(
        "total_audit_failures\t"
        f"{total_failures}"
    )

    # Any discrepancy results in failure and Unix exit status 1.
    if total_failures > 0:
        print(
            "simulation_audit_status\tFAIL"
        )

        return 1

    # No discrepancies means the complete simulation chain is internally
    # consistent.
    print(
        "simulation_audit_status\tPASS"
    )

    # Unix success.
    return 0


# Run the audit only when the file is executed directly.
#
# Importing it later for tests will not automatically launch the audit.
if __name__ == "__main__":
    raise SystemExit(main())
