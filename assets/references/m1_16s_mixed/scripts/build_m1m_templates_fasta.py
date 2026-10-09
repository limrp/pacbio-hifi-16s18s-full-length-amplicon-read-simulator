#!/usr/bin/env python3

import csv
import hashlib
from pathlib import Path


BASE = Path("data/m1_work/mixed_16s")

SOURCE_FASTA = (
    BASE
    / "m1_mixed_16s_selected_full_length.fasta"
)

MANIFEST = (
    BASE
    / "final_truth_build"
    / "m1m_truth_manifest.tsv"
)

OUTPUT_FASTA = (
    BASE
    / "final_truth_build"
    / "m1m_16s_templates.fasta"
)

OUTPUT_HASHES = (
    BASE
    / "final_truth_build"
    / "m1m_16s_templates_sha256.tsv"
)


def read_fasta(path):
    records = []
    header = None
    chunks = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        (
                            header.split()[0],
                            "".join(chunks).upper(),
                        )
                    )

                header = line[1:]
                chunks = []

            else:
                if header is None:
                    raise ValueError(
                        f"Sequence before FASTA header in {path}"
                    )

                chunks.append(line)

    if header is not None:
        records.append(
            (
                header.split()[0],
                "".join(chunks).upper(),
            )
        )

    return records


def read_manifest(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        rows = list(reader)

    result = {}

    for row in rows:
        source_id = row["source_sequence_id"]

        if source_id in result:
            raise ValueError(
                f"Duplicate source_sequence_id in manifest: {source_id}"
            )

        result[source_id] = row

    return result


def sha256_sequence(sequence):
    return hashlib.sha256(
        sequence.encode("ascii")
    ).hexdigest()


def wrap(sequence, width=80):
    for start in range(0, len(sequence), width):
        yield sequence[start:start + width]


def main():

    fasta_records = read_fasta(SOURCE_FASTA)
    manifest = read_manifest(MANIFEST)

    fasta_ids = [source_id for source_id, _ in fasta_records]

    if len(fasta_records) != 20:
        raise ValueError(
            f"Expected 20 FASTA records; found {len(fasta_records)}"
        )

    if len(set(fasta_ids)) != len(fasta_ids):
        raise ValueError(
            "Duplicate source IDs detected in FASTA"
        )

    if set(fasta_ids) != set(manifest):
        raise ValueError(
            "FASTA source IDs and manifest source IDs do not match"
        )

    seen_template_ids = set()
    hash_rows = []

    with OUTPUT_FASTA.open("w") as fasta_out:

        for source_id, sequence in fasta_records:

            row = manifest[source_id]

            template_id = row["template_id"]
            expected_length = int(row["sequence_length"])
            expected_sha256 = row["sequence_sha256"]

            observed_length = len(sequence)
            observed_sha256 = sha256_sequence(sequence)

            if observed_length != expected_length:
                raise ValueError(
                    f"Length mismatch for {source_id}: "
                    f"{observed_length} != {expected_length}"
                )

            if observed_sha256 != expected_sha256:
                raise ValueError(
                    f"SHA-256 mismatch for {source_id}"
                )

            if template_id in seen_template_ids:
                raise ValueError(
                    f"Duplicate stable template ID: {template_id}"
                )

            seen_template_ids.add(template_id)

            header = (
                f"{template_id} "
                f"source_sequence_id={source_id} "
                f"organism_id={row['organism_id']} "
                f"domain={row['domain']}"
            )

            fasta_out.write(f">{header}\n")

            for sequence_line in wrap(sequence):
                fasta_out.write(sequence_line + "\n")

            hash_rows.append(
                {
                    "template_id": template_id,
                    "source_sequence_id": source_id,
                    "sequence_length": observed_length,
                    "sequence_sha256": observed_sha256,
                }
            )

    with OUTPUT_HASHES.open("w", newline="") as handle:

        fieldnames = [
            "template_id",
            "source_sequence_id",
            "sequence_length",
            "sequence_sha256",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(hash_rows)

    print(f"Wrote {len(fasta_records)} stable templates.")
    print(f"Reference FASTA: {OUTPUT_FASTA}")
    print(f"Sequence hashes: {OUTPUT_HASHES}")


if __name__ == "__main__":
    main()
