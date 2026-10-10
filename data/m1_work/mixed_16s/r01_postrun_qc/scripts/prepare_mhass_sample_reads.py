#!/usr/bin/env python3

"""
Prepare analysis-ready FASTQ reads from controlled MHASS simulations.

PURPOSE
-------
MHASS surrounds each biological truth sequence with artificial sequence:

    "A"
    + ForwardBarcode
    + biological truth sequence
    + RevCompReverse barcode
    + "A"

The realized CCS read may contain substitutions, insertions, deletions,
or terminal sequence loss relative to that expected construct.

Therefore, simply removing a fixed number of bases such as:

    read[35:-35]

is not reliable.

Instead, this script:

1. Recovers each read's known provenance from the MHASS mapping table.
2. Reconstructs the expected MHASS molecule for that read.
3. Determines whether the observed FASTQ read is already in the same
   orientation as the truth template or must be reverse-complemented.
4. Globally aligns the canonical observed read to the expected molecule.
5. Projects the known biological reference interval onto observed-read
   coordinates using the validated CIGAR walker.
6. Removes only the artificial MHASS wrapper.
7. Preserves the actually observed biological sequence, including
   sequencing errors.
8. Preserves the corresponding FASTQ quality scores.
9. Writes one prepared FASTQ per sample.
10. Writes one provenance/truth row per prepared read.

IMPORTANT
---------
The biological truth sequence is used only to determine provenance,
orientation, and alignment-derived boundaries.

The truth sequence is NEVER substituted for the observed read sequence.
"""

# Allow modern Python type annotations to be evaluated consistently.
from __future__ import annotations

# argparse parses command-line arguments such as --truth and --fastq.
import argparse

# csv provides safe reading/writing of tab-separated metadata files.
import csv

# Counter lets us count output reads per sample and other categories.
from collections import Counter

# Path provides safer filesystem-path handling than raw strings.
from pathlib import Path

# edlib performs the same global edit-distance/path alignment that we
# already validated in evaluate_mhass_preparation_strategy_v02.py.
import edlib

# SeqIO reads and writes FASTA/FASTQ records.
from Bio import SeqIO

# Seq is used to construct output sequence objects.
from Bio.Seq import Seq

# SeqRecord is used to create each prepared FASTQ record while retaining
# the original read ID and description.
from Bio.SeqRecord import SeqRecord

# Reuse the CIGAR-coordinate walker that already passed our artificial
# tests and the full-R01 accounting invariants.
from mhass_cigar_utils import walk_cigar

# Reuse the FASTQ orientation/slicing helpers that passed all 9 tests.
from mhass_fastq_utils import (
    orient_sequence_and_qualities,
    slice_sequence_and_qualities,
)


def parse_args() -> argparse.Namespace:
    """
    Define and parse all command-line arguments.

    Keeping analysis metadata such as marker, replicate, and seed as
    explicit command-line arguments makes the output provenance
    reproducible instead of hiding those values inside the script.
    """

    # Create the command-line argument parser.
    parser = argparse.ArgumentParser(
        description=(
            "Prepare analysis-ready MHASS reads by reconstructing provenance, "
            "canonicalizing orientation, projecting biological boundaries, "
            "and removing only the artificial MHASS wrapper."
        )
    )

    # FASTA containing the biological truth templates.
    parser.add_argument(
        "--truth",
        required=True,
        type=Path,
        help="FASTA containing biological truth templates.",
    )

    # MHASS mapping:
    # TemplateFile -> ASVID / SampleID / BarcodeID.
    parser.add_argument(
        "--mapping",
        required=True,
        type=Path,
        help="MHASS sequence_file_mapping.tsv.",
    )

    # MHASS barcode table:
    # SampleID -> BarcodeID / ForwardBarcode / RevCompReverse.
    parser.add_argument(
        "--barcodes",
        required=True,
        type=Path,
        help="MHASS sample_barcode_map.tsv.",
    )

    # Raw realized MHASS CCS FASTQ.
    parser.add_argument(
        "--fastq",
        required=True,
        type=Path,
        help="Raw MHASS combined_reads.fastq.",
    )

    # New directory in which prepared FASTQs and truth_reads.tsv
    # will be written.
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="New output directory. It must not already exist.",
    )

    # Metadata copied into truth_reads.tsv.
    parser.add_argument(
        "--marker",
        required=True,
        help="Marker name, for example 16S.",
    )

    parser.add_argument(
        "--dataset-id",
        required=True,
        help="Dataset identifier, for example toyctrl16.",
    )

    parser.add_argument(
        "--simulation-replicate",
        required=True,
        help="Simulation replicate, for example R01.",
    )

    parser.add_argument(
        "--master-seed",
        required=True,
        type=int,
        help="MHASS master random seed, for example 16001.",
    )

    # This argument is only for staged validation.
    #
    # Example:
    # --max-reads 10
    #
    # processes only the first ten FASTQ records.
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Process only the first N reads. Default: process all reads.",
    )

    # Parse the actual command-line values and return them.
    return parser.parse_args()


def load_truth(path: Path) -> dict[str, str]:
    """
    Load the truth FASTA into:

        truth_template_id -> uppercase sequence

    Duplicate FASTA IDs are forbidden because provenance must map each
    identifier to exactly one biological truth sequence.
    """

    # Create an initially empty dictionary.
    truth: dict[str, str] = {}

    # Read every FASTA record.
    for record in SeqIO.parse(path, "fasta"):

        # Duplicate template IDs would make ground truth ambiguous.
        if record.id in truth:
            raise ValueError(
                f"Duplicate truth FASTA ID: {record.id}"
            )

        # Store the sequence in uppercase so comparisons are consistent.
        truth[record.id] = str(record.seq).upper()

    # An empty truth FASTA cannot support preparation.
    if not truth:
        raise ValueError(
            "Truth FASTA contains no sequences."
        )

    # Return the complete lookup dictionary.
    return truth


def load_barcodes(
    path: Path,
) -> dict[str, dict[str, str]]:
    """
    Load the MHASS barcode table.

    Returned structure:

        barcodes[sample_id] = {
            "barcode_id": ...,
            "forward": ...,
            "revcomp_reverse": ...,
        }
    """

    # Create an empty sample -> barcode-information dictionary.
    barcodes: dict[str, dict[str, str]] = {}

    # Open the TSV without newline translation.
    with path.open(newline="") as handle:

        # DictReader lets us refer to columns by their names rather
        # than fragile numeric positions.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # These are the columns required by our preparation algorithm.
        required = {
            "SampleID",
            "BarcodeID",
            "ForwardBarcode",
            "RevCompReverse",
        }

        # Determine whether any required column is missing.
        missing = required - set(reader.fieldnames or [])

        # Fail rather than guessing if the table structure differs.
        if missing:
            raise ValueError(
                "Barcode table is missing required columns: "
                f"{sorted(missing)}. "
                f"Found columns: {reader.fieldnames}"
            )

        # Process one barcode-table row at a time.
        for row in reader:

            # SampleID is our primary key.
            sample_id = row["SampleID"]

            # A sample must map to exactly one barcode definition.
            if sample_id in barcodes:
                raise ValueError(
                    f"Duplicate barcode entry for sample: {sample_id}"
                )

            # Avoid accidental path creation such as "S01/foo".
            # Our current samples S01/S02/S03 pass this check.
            if "/" in sample_id or "\\" in sample_id:
                raise ValueError(
                    f"Unsafe SampleID for output filename: {sample_id}"
                )

            # Store normalized uppercase barcode sequences.
            barcodes[sample_id] = {
                "barcode_id": row["BarcodeID"],
                "forward": row["ForwardBarcode"].upper(),
                "revcomp_reverse": row["RevCompReverse"].upper(),
            }

    # Return all sample/barcode metadata.
    return barcodes


def load_mapping(
    path: Path,
) -> dict[str, tuple[str, str, str, str]]:
    """
    Load MHASS sequence_file_mapping.tsv.

    Returned structure:

        template_group -> (
            truth_template_id,
            sample_id,
            barcode_id,
            original_template_file,
        )

    FASTQ headers contain a prefix such as:

        template10_A1_np17/S/3/ccs

    therefore the template group is:

        template10_A1_np17

    MHASS TemplateFile contains the corresponding:

        template10_A1_np17.fasta
    """

    # Dictionary indexed by the FASTQ template-group prefix.
    mapping: dict[str, tuple[str, str, str, str]] = {}

    # Open the mapping TSV.
    with path.open(newline="") as handle:

        # Read columns by header name.
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        # Required provenance fields.
        required = {
            "ASVID",
            "SampleID",
            "BarcodeID",
            "TemplateFile",
        }

        # Determine whether the file has every required column.
        missing = required - set(reader.fieldnames or [])

        # Stop if the file schema is not what we expect.
        if missing:
            raise ValueError(
                "Mapping table is missing required columns: "
                f"{sorted(missing)}. "
                f"Found columns: {reader.fieldnames}"
            )

        # Process every intended-molecule mapping row.
        for row in reader:

            # Preserve the original filename for truth_reads.tsv.
            original_template_file = row["TemplateFile"]

            # The FASTQ read prefix does not contain ".fasta",
            # so remove that suffix for dictionary lookup.
            template_group = original_template_file

            if template_group.endswith(".fasta"):
                template_group = template_group[:-6]

            # Provenance that all repeated rows for this template group
            # must agree on.
            value = (
                row["ASVID"],
                row["SampleID"],
                row["BarcodeID"],
                original_template_file,
            )

            # MHASS may repeat the same TemplateFile for multiple intended
            # molecule rows. Repetition is allowed only when provenance
            # is identical.
            if (
                template_group in mapping
                and mapping[template_group] != value
            ):
                raise ValueError(
                    "Conflicting provenance for template group "
                    f"{template_group}: "
                    f"{mapping[template_group]} vs {value}"
                )

            # Store or reaffirm the mapping.
            mapping[template_group] = value

    # Return template-group provenance.
    return mapping


def build_expected_construct(
    truth_sequence: str,
    forward_barcode: str,
    revcomp_reverse_barcode: str,
) -> str:
    """
    Reconstruct the exact molecule that MHASS intended to simulate.

    Structure:

        A
        + ForwardBarcode
        + biological truth sequence
        + RevCompReverse
        + A
    """

    # Concatenate the exact expected components in MHASS order.
    return (
        "A"
        + forward_barcode
        + truth_sequence
        + revcomp_reverse_barcode
        + "A"
    )


def edit_distance(
    query: str,
    reference: str,
) -> int:
    """
    Compute global edit distance between observed query and expected R.

    mode="NW" means global Needleman-Wunsch-style alignment:
    the full query and full reference participate in the comparison.
    """

    # Ask edlib only for the minimum edit distance.
    result = edlib.align(
        query,
        reference,
        mode="NW",
        task="distance",
    )

    # Extract the numeric minimum edit distance.
    distance = result["editDistance"]

    # Negative distance indicates alignment failure.
    if distance < 0:
        raise RuntimeError(
            "edlib failed to produce a global edit distance."
        )

    # Return the valid edit distance.
    return distance


def alignment_path(
    query: str,
    reference: str,
) -> tuple[int, str]:
    """
    Run the same global alignment but request the CIGAR path.

    Returns:

        (edit_distance, cigar)
    """

    # task="path" asks edlib for both distance and alignment CIGAR.
    result = edlib.align(
        query,
        reference,
        mode="NW",
        task="path",
    )

    # Extract the independently reported path edit distance.
    distance = result["editDistance"]

    # Extract the extended CIGAR string.
    cigar = result["cigar"]

    # Stop on alignment failure.
    if distance < 0:
        raise RuntimeError(
            "edlib failed to produce a global alignment."
        )

    # Boundary projection cannot work without a CIGAR path.
    if not cigar:
        raise RuntimeError(
            "edlib did not return a CIGAR alignment path."
        )

    # Return both values for consistency checking.
    return distance, cigar


def main() -> int:
    """
    Execute controlled MHASS read preparation.
    """

    # Read command-line arguments.
    args = parse_args()

    # Refuse accidental overwriting.
    #
    # Every test run should therefore use a new output directory.
    if args.outdir.exists():
        raise FileExistsError(
            f"Output directory already exists: {args.outdir}"
        )

    # Create the new output directory.
    args.outdir.mkdir(
        parents=True,
        exist_ok=False,
    )

    # Load all reference/provenance information once.
    truth = load_truth(args.truth)
    barcodes = load_barcodes(args.barcodes)
    mapping = load_mapping(args.mapping)

    # Define the per-read provenance table path.
    truth_reads_path = args.outdir / "truth_reads.tsv"

    # These columns document exactly how each output read was prepared.
    truth_fieldnames = [
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
    ]

    # Keep one open FASTQ handle per sample.
    #
    # We open each sample file only when that sample first appears.
    sample_handles: dict[str, object] = {}

    # Count prepared reads per sample.
    sample_counts: Counter[str] = Counter()

    # Count how many reads used each orientation.
    orientation_counts: Counter[str] = Counter()

    # Count how many boundary-ambiguous reads were retained.
    boundary_counts: Counter[str] = Counter()

    # Record read IDs so duplicated FASTQ identifiers cause a hard failure.
    seen_read_ids: set[str] = set()

    # Total number of successfully prepared reads.
    processed = 0

    # The truth TSV is managed by a context manager so it closes safely.
    with truth_reads_path.open(
        "w",
        newline="",
    ) as truth_handle:

        # Create the tab-separated writer.
        truth_writer = csv.DictWriter(
            truth_handle,
            fieldnames=truth_fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        # Write the truth_reads.tsv header.
        truth_writer.writeheader()

        try:
            # Read every raw MHASS FASTQ record.
            for read_index, record in enumerate(
                SeqIO.parse(args.fastq, "fastq"),
                start=1,
            ):

                # Optional smoke-test limit.
                if (
                    args.max_reads is not None
                    and read_index > args.max_reads
                ):
                    break

                # Biopython's record.id is the first token in the
                # original FASTQ header.
                read_id = record.id

                # Duplicate IDs would break our required one-read/one-truth-row
                # relationship, so fail immediately.
                if read_id in seen_read_ids:
                    raise ValueError(
                        f"Duplicate FASTQ read ID: {read_id}"
                    )

                # Mark the read ID as observed.
                seen_read_ids.add(read_id)

                # MHASS read IDs have the form:
                #
                # template10_A1_np17/S/3/ccs
                #
                # Everything before the first slash identifies the
                # grouped input-template FASTA.
                template_group = read_id.split("/", 1)[0]

                # Provenance must be known.
                if template_group not in mapping:
                    raise ValueError(
                        "Unknown FASTQ template group: "
                        f"{template_group}"
                    )

                # Recover validated provenance.
                (
                    truth_template_id,
                    sample_id,
                    barcode_id,
                    original_template_file,
                ) = mapping[template_group]

                # The truth template must exist.
                if truth_template_id not in truth:
                    raise ValueError(
                        "Truth template missing from FASTA: "
                        f"{truth_template_id}"
                    )

                # The mapped sample must exist in the barcode table.
                if sample_id not in barcodes:
                    raise ValueError(
                        "Sample missing from barcode table: "
                        f"{sample_id}"
                    )

                # Retrieve this sample's barcode metadata.
                barcode_info = barcodes[sample_id]

                # Mapping-table BarcodeID and sample-barcode BarcodeID
                # must agree exactly.
                if (
                    barcode_info["barcode_id"]
                    != barcode_id
                ):
                    raise ValueError(
                        f"Barcode mismatch for {template_group}: "
                        f"mapping={barcode_id}, "
                        f"sample table="
                        f"{barcode_info['barcode_id']}"
                    )

                # Retrieve the biological truth sequence.
                truth_sequence = truth[truth_template_id]

                # Reconstruct the exact expected MHASS construct R.
                expected = build_expected_construct(
                    truth_sequence=truth_sequence,
                    forward_barcode=barcode_info["forward"],
                    revcomp_reverse_barcode=(
                        barcode_info["revcomp_reverse"]
                    ),
                )

                # The biological target begins after:
                #
                # 1 terminal artificial A
                # +
                # the forward barcode.
                target_start_reference = (
                    1 + len(barcode_info["forward"])
                )

                # The biological target ends after the complete truth sequence.
                target_end_reference = (
                    target_start_reference
                    + len(truth_sequence)
                )

                # Validate that the expected right wrapper has exactly:
                #
                # RevCompReverse barcode + terminal A
                expected_right_wrapper_length = (
                    len(barcode_info["revcomp_reverse"])
                    + 1
                )

                # If this relationship ever fails, our construct-building
                # and coordinate definitions have diverged.
                if (
                    len(expected) - target_end_reference
                    != expected_right_wrapper_length
                ):
                    raise RuntimeError(
                        "Expected-construct structure inconsistency "
                        f"for {read_id}"
                    )

                # Store the raw observed sequence in uppercase.
                raw_sequence = str(record.seq).upper()

                # Retrieve the parsed Phred quality score for each raw base.
                #
                # The quality list must remain synchronized with sequence
                # coordinates throughout every transformation.
                raw_qualities = record.letter_annotations.get(
                    "phred_quality"
                )

                # FASTQ input must actually contain quality scores.
                if raw_qualities is None:
                    raise ValueError(
                        f"FASTQ qualities missing for read: {read_id}"
                    )

                # Convert to a normal list so orientation/slicing helpers
                # receive exactly the tested datatype.
                raw_qualities = list(raw_qualities)

                # A FASTQ record must have exactly one quality value per base.
                if len(raw_sequence) != len(raw_qualities):
                    raise ValueError(
                        "Raw sequence/quality length mismatch for "
                        f"{read_id}: "
                        f"sequence={len(raw_sequence)}, "
                        f"qualities={len(raw_qualities)}"
                    )

                # Construct the reverse-complement candidate solely for
                # orientation comparison.
                reverse_sequence = str(
                    Seq(raw_sequence).reverse_complement()
                )

                # Compare both possible orientations to expected R.
                d_as_read = edit_distance(
                    raw_sequence,
                    expected,
                )

                d_revcomp = edit_distance(
                    reverse_sequence,
                    expected,
                )

                # Select whichever orientation has the smaller global
                # edit distance.
                if d_as_read < d_revcomp:
                    orientation = "as_read"
                    best_distance = d_as_read

                elif d_revcomp < d_as_read:
                    orientation = "revcomp"
                    best_distance = d_revcomp

                # A tie means orientation cannot be determined safely.
                # Production code must fail rather than guess.
                else:
                    raise RuntimeError(
                        "Ambiguous orientation for read: "
                        f"{read_id}"
                    )

                # Count the selected orientation.
                orientation_counts[orientation] += 1

                # Apply the selected orientation to BOTH sequence and
                # quality scores using our independently tested helper.
                (
                    canonical_sequence,
                    canonical_qualities,
                ) = orient_sequence_and_qualities(
                    raw_sequence,
                    raw_qualities,
                    orientation,
                )

                # Obtain the selected global CIGAR path.
                path_distance, cigar = alignment_path(
                    canonical_sequence,
                    expected,
                )

                # The path alignment must reproduce the same minimum
                # distance used to select orientation.
                if path_distance != best_distance:
                    raise RuntimeError(
                        "Distance/path inconsistency for "
                        f"{read_id}: "
                        f"orientation distance={best_distance}, "
                        f"path distance={path_distance}"
                    )

                # Walk the validated CIGAR path and project biological
                # reference boundaries onto observed-query coordinates.
                boundary_qc = walk_cigar(
                    cigar,
                    reference_length=len(expected),
                    query_length=len(canonical_sequence),
                    target_start=target_start_reference,
                    target_end=target_end_reference,
                )

                # Our validated policy keeps insertions that fall exactly
                # at either artificial/biological boundary:
                #
                # left  -> choose minimum possible start
                # right -> choose maximum possible end
                chosen_start = int(
                    boundary_qc[
                        "alignment_target_start_min"
                    ]
                )

                chosen_end = int(
                    boundary_qc[
                        "alignment_target_end_max"
                    ]
                )

                # Retrieve explicit ambiguity flags from the CIGAR walker.
                left_boundary_ambiguous = bool(
                    boundary_qc[
                        "left_boundary_ambiguous"
                    ]
                )

                right_boundary_ambiguous = bool(
                    boundary_qc[
                        "right_boundary_ambiguous"
                    ]
                )

                # If either boundary contains an insertion ambiguity,
                # record that we intentionally retained it.
                if (
                    left_boundary_ambiguous
                    or right_boundary_ambiguous
                ):
                    boundary_status = "ambiguous_retained"

                # Otherwise both boundaries were uniquely projected.
                else:
                    boundary_status = "ok"

                # Count each boundary-status category.
                boundary_counts[boundary_status] += 1

                # Slice sequence and quality scores using exactly the
                # same half-open interval [chosen_start, chosen_end).
                (
                    prepared_sequence,
                    prepared_qualities,
                ) = slice_sequence_and_qualities(
                    canonical_sequence,
                    canonical_qualities,
                    chosen_start,
                    chosen_end,
                )

                # Independently predict how long the observed biological
                # interval should be from CIGAR operations.
                #
                # Start from truth-reference length...
                predicted_prepared_length = len(
                    truth_sequence
                )

                # ...add observed insertions inside the biological target...
                predicted_prepared_length += int(
                    boundary_qc["target_insertions"]
                )

                # ...subtract reference bases deleted inside the target...
                predicted_prepared_length -= int(
                    boundary_qc["target_deletions"]
                )

                # ...and retain insertions that occur exactly at the
                # left artificial/biological boundary.
                predicted_prepared_length += int(
                    boundary_qc[
                        "left_boundary_insertions"
                    ]
                )

                # Also retain insertions exactly at the right boundary.
                predicted_prepared_length += int(
                    boundary_qc[
                        "right_boundary_insertions"
                    ]
                )

                # Production preparation must satisfy the same extraction
                # invariant already validated across all 10,228 R01 reads.
                if (
                    len(prepared_sequence)
                    != predicted_prepared_length
                ):
                    raise RuntimeError(
                        "Prepared-length invariant failed for "
                        f"{read_id}: "
                        f"observed={len(prepared_sequence)}, "
                        f"predicted={predicted_prepared_length}"
                    )

                # The FASTQ invariant must hold after every transformation:
                #
                # one sequence base <-> one quality score.
                if (
                    len(prepared_sequence)
                    != len(prepared_qualities)
                ):
                    raise RuntimeError(
                        "Prepared sequence/quality mismatch for "
                        f"{read_id}"
                    )

                # The output filename is deterministic from SampleID.
                output_fastq_name = (
                    f"{sample_id}.fastq"
                )

                # Construct its complete path inside the new output directory.
                output_fastq_path = (
                    args.outdir
                    / output_fastq_name
                )

                # Open that sample's FASTQ only the first time the sample
                # is encountered.
                if sample_id not in sample_handles:
                    sample_handles[sample_id] = (
                        output_fastq_path.open("w")
                    )

                # Create the prepared FASTQ SeqRecord.
                #
                # IMPORTANT:
                # - id remains the original read ID.
                # - description remains the original FASTQ description.
                # - sequence is the observed prepared sequence, not truth.
                prepared_record = SeqRecord(
                    Seq(prepared_sequence),
                    id=record.id,
                    name=record.name,
                    description=record.description,
                )

                # Attach the prepared qualities to the prepared sequence.
                prepared_record.letter_annotations[
                    "phred_quality"
                ] = prepared_qualities

                # Write exactly one prepared FASTQ record.
                records_written = SeqIO.write(
                    prepared_record,
                    sample_handles[sample_id],
                    "fastq",
                )

                # SeqIO.write() should report one written record here.
                if records_written != 1:
                    raise RuntimeError(
                        "Unexpected FASTQ write count for "
                        f"{read_id}: {records_written}"
                    )

                # Write the corresponding provenance/truth row.
                truth_writer.writerow(
                    {
                        "read_id": read_id,
                        "marker": args.marker,
                        "dataset_id": args.dataset_id,
                        "sample_id": sample_id,
                        "simulation_replicate": (
                            args.simulation_replicate
                        ),
                        "master_seed": args.master_seed,
                        "truth_template_id": (
                            truth_template_id
                        ),
                        "barcode_id": barcode_id,
                        "mhass_template_file": (
                            original_template_file
                        ),
                        "input_orientation": orientation,
                        "boundary_status": boundary_status,
                        "left_boundary_ambiguous": (
                            left_boundary_ambiguous
                        ),
                        "right_boundary_ambiguous": (
                            right_boundary_ambiguous
                        ),
                        "target_start": chosen_start,
                        "target_end": chosen_end,
                        "raw_length": len(raw_sequence),
                        "prepared_length": len(
                            prepared_sequence
                        ),
                        "output_fastq": (
                            output_fastq_name
                        ),
                    }
                )

                # Update per-sample output counts.
                sample_counts[sample_id] += 1

                # Update total successfully prepared reads.
                processed += 1

        finally:
            # Whether preparation succeeds or fails, close every sample
            # FASTQ handle that was opened.
            for handle in sample_handles.values():
                handle.close()

    # An empty FASTQ or an unexpectedly zero-sized smoke test is invalid.
    if processed == 0:
        raise RuntimeError(
            "No FASTQ reads were prepared."
        )

    # Print a compact reproducibility/QC summary.
    print(f"processed_reads\t{processed}")

    # Print per-orientation counts in deterministic order.
    print(
        "as_read\t"
        f"{orientation_counts['as_read']}"
    )

    print(
        "revcomp\t"
        f"{orientation_counts['revcomp']}"
    )

    # Print boundary-policy counts.
    print(
        "boundary_ok\t"
        f"{boundary_counts['ok']}"
    )

    print(
        "boundary_ambiguous_retained\t"
        f"{boundary_counts['ambiguous_retained']}"
    )

    # Print sample counts alphabetically for deterministic output.
    for sample_id in sorted(sample_counts):
        print(
            f"sample_{sample_id}\t"
            f"{sample_counts[sample_id]}"
        )

    # Print the output paths for convenient logging.
    print(f"truth_reads\t{truth_reads_path}")
    print(f"outdir\t{args.outdir}")

    # Conventional successful program exit status.
    return 0


# Only run main() when this file is executed directly.
#
# Importing functions from this module in tests will therefore NOT
# automatically start FASTQ preparation.
if __name__ == "__main__":
    raise SystemExit(main())
