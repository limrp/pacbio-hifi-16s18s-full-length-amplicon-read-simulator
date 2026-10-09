#!/usr/bin/env python3

import csv
import hashlib
from pathlib import Path


BASE = Path("data/m1_work/mixed_16s")

FASTA = BASE / "m1_mixed_16s_selected_full_length.fasta"

SOURCE_META = (
    BASE
    / "final_truth_build"
    / "final_source_metadata.tsv"
)

UC = (
    BASE
    / "final_panel_qc"
    / "clusters_97.uc"
)

OUTPUT = (
    BASE
    / "final_truth_build"
    / "m1m_truth_manifest.tsv"
)


def read_fasta(path):
    records = []
    current_id = None
    chunks = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records.append(
                        (current_id, "".join(chunks).upper())
                    )

                current_id = line[1:].split()[0]
                chunks = []

            else:
                if current_id is None:
                    raise ValueError(
                        f"Sequence found before FASTA header in {path}"
                    )

                chunks.append(line)

    if current_id is not None:
        records.append(
            (current_id, "".join(chunks).upper())
        )

    return records


def read_source_metadata(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        rows = list(reader)

    required = {
        "source_sequence_id",
        "domain",
        "organism_name",
        "source_type",
        "source_accession",
    }

    if set(reader.fieldnames or []) != required:
        raise ValueError(
            f"Unexpected metadata columns: {reader.fieldnames}"
        )

    result = {}

    for row in rows:
        source_id = row["source_sequence_id"]

        if source_id in result:
            raise ValueError(
                f"Duplicate metadata ID: {source_id}"
            )

        result[source_id] = row

    return result


def read_uc(path):
    """
    Return:

        sequence_id ->
            raw cluster ID
            centroid source sequence ID
    """

    membership = {}
    centroids = {}

    rows = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            record_type = fields[0]

            if record_type not in {"S", "H"}:
                continue

            cluster_id = fields[1]
            sequence_id = fields[8]

            rows.append(
                (record_type, cluster_id, sequence_id, fields)
            )

            if record_type == "S":
                centroids[cluster_id] = sequence_id

    for record_type, cluster_id, sequence_id, fields in rows:

        if cluster_id not in centroids:
            raise ValueError(
                f"No centroid found for cluster {cluster_id}"
            )

        membership[sequence_id] = {
            "vsearch_cluster_id": cluster_id,
            "centroid_source_id": centroids[cluster_id],
        }

    return membership


def sha256_sequence(sequence):
    return hashlib.sha256(
        sequence.encode("ascii")
    ).hexdigest()


def main():

    fasta_records = read_fasta(FASTA)
    metadata = read_source_metadata(SOURCE_META)
    clusters = read_uc(UC)

    fasta_ids = [record_id for record_id, _ in fasta_records]

    if len(fasta_ids) != 20:
        raise ValueError(
            f"Expected 20 FASTA records, found {len(fasta_ids)}"
        )

    if len(set(fasta_ids)) != len(fasta_ids):
        raise ValueError("Duplicate FASTA IDs detected")

    if set(fasta_ids) != set(metadata):
        raise ValueError(
            "FASTA IDs and metadata IDs do not match exactly"
        )

    if set(fasta_ids) != set(clusters):
        raise ValueError(
            "FASTA IDs and UC clustering IDs do not match exactly"
        )

    # Assign stable organism IDs in first-appearance order.
    organism_ids = {}

    # Assign stable OTU IDs in first-appearance order.
    otu_ids = {}

    rows = []

    for index, (source_id, sequence) in enumerate(
        fasta_records,
        start=1,
    ):

        meta = metadata[source_id]

        organism_name = meta["organism_name"]

        if organism_name not in organism_ids:
            organism_ids[organism_name] = (
                f"M1M_ORG{len(organism_ids) + 1:03d}"
            )

        raw_cluster = clusters[source_id]["vsearch_cluster_id"]

        if raw_cluster not in otu_ids:
            otu_ids[raw_cluster] = (
                f"M1M_OTU97_{len(otu_ids) + 1:03d}"
            )

        rows.append(
            {
                "template_id": f"M1M_16S_T{index:03d}",
                "marker": "16S",
                "template_sequence_type": "full_length_16S",
                "domain": meta["domain"],
                "organism_id": organism_ids[organism_name],
                "organism_name": organism_name,
                "source_type": meta["source_type"],
                "source_accession": meta["source_accession"],
                "source_sequence_id": source_id,
                "sequence_length": len(sequence),
                "sequence_sha256": sha256_sequence(sequence),
                "exact_variant_id": f"M1M_16S_VAR{index:03d}",
                "otu97_id": otu_ids[raw_cluster],
                "vsearch_cluster_id": raw_cluster,
                "otu97_centroid_source_id": (
                    clusters[source_id]["centroid_source_id"]
                ),
                "selected_final": "yes",
            }
        )

    fieldnames = list(rows[0])

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} templates.")
    print(f"Organisms: {len(organism_ids)}")
    print(f"97% OTUs: {len(otu_ids)}")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
