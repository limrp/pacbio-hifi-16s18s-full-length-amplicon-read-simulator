#!/usr/bin/env python3

import csv
from pathlib import Path


BASE = Path(
    "data/m1_work/mixed_16s/final_truth_build"
)

PRIMER_DIR = BASE / "primer_compatibility"

CORE_MANIFEST = BASE / "m1m_truth_manifest.tsv"
TAXONOMY_QUERIES = BASE / "m1m_taxonomy_queries.tsv"

TAXONOMY_SUMMARY = (
    BASE
    / "m1m_ncbi_taxonomy_2026-09-30"
    / "ncbi_dataset"
    / "data"
    / "taxonomy_summary.tsv"
)

LOCUS_PROVENANCE = (
    BASE
    / "m1m_genome_locus_provenance.tsv"
)

PROJECT_QC = (
    PRIMER_DIR
    / "project_primer_qc.tsv"
)

ARCHAEA_FORWARD_QC = (
    PRIMER_DIR
    / "archaea_primer_qc_window50.tsv"
)

ARCHAEA_REVERSE_QC = (
    PRIMER_DIR
    / "archaea_primer_qc_window200.tsv"
)

PROJECT_EXACT_BED = (
    PRIMER_DIR
    / "project_exact_amplicons.bed"
)

ARCHAEA_EXACT_BED = (
    PRIMER_DIR
    / "archaea_exact_amplicons.bed"
)

OUTPUT = (
    BASE
    / "m1m_truth_manifest_enriched.tsv"
)


PROJECT_FORWARD = "AGRGTTYGATYMTGGCTCAG"
PROJECT_REVERSE = "RGYTACCTTGTTACGACTT"

ARCHAEA_FORWARD = "TCCGGTTGATCCYGCCGG"
ARCHAEA_REVERSE = "CRGTGWGTRCAAGGRGCA"


def read_tsv(path):
    """Read one TSV and return (fieldnames, rows)."""

    if not path.exists():
        raise FileNotFoundError(
            f"Required file does not exist: {path}"
        )

    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        rows = list(reader)

    if not reader.fieldnames:
        raise ValueError(
            f"No header found in: {path}"
        )

    return reader.fieldnames, rows


def require_columns(fieldnames, required, label):
    """Fail if an expected column is missing."""

    missing = sorted(
        set(required) - set(fieldnames)
    )

    if missing:
        raise ValueError(
            f"{label} is missing columns: {missing}"
        )


def unique_map(rows, key, label):
    """
    Convert rows into:
        value_of_key -> row

    Fail on duplicate keys.
    """

    result = {}

    for row in rows:
        value = row[key]

        if not value:
            raise ValueError(
                f"Empty {key} in {label}"
            )

        if value in result:
            raise ValueError(
                f"Duplicate {key}={value!r} "
                f"in {label}"
            )

        result[value] = row

    return result


def read_primer_side(
    path,
    primer_side,
    expected_primer,
    expected_template_ids,
    expected_lengths,
):
    """
    Read one selected primer side from one QC file.

    This lets us deliberately use:

      project pair:
        window100 for forward + reverse

      archaeal pair:
        window50  for forward
        window200 for reverse
    """

    fieldnames, rows = read_tsv(path)

    require_columns(
        fieldnames,
        {
            "sequence_id",
            "sequence_length",
            "primer",
            "primer_iupac",
            "distance_from_expected_end",
            "best_mismatches",
            "three_prime_mismatches",
        },
        str(path),
    )

    selected = [
        row
        for row in rows
        if row["primer"] == primer_side
    ]

    result = unique_map(
        selected,
        "sequence_id",
        f"{path}:{primer_side}",
    )

    if set(result) != expected_template_ids:
        missing = sorted(
            expected_template_ids - set(result)
        )

        extra = sorted(
            set(result) - expected_template_ids
        )

        raise ValueError(
            f"Primer-QC template mismatch "
            f"for {path}:{primer_side}\n"
            f"Missing: {missing}\n"
            f"Extra: {extra}"
        )

    for template_id, row in result.items():

        if row["primer_iupac"] != expected_primer:
            raise ValueError(
                f"Unexpected {primer_side} primer "
                f"for {template_id}: "
                f"{row['primer_iupac']}"
            )

        observed_length = int(
            row["sequence_length"]
        )

        if observed_length != expected_lengths[template_id]:
            raise ValueError(
                f"Length mismatch in primer QC "
                f"for {template_id}: "
                f"{observed_length} != "
                f"{expected_lengths[template_id]}"
            )

    return result


def read_bed_ids(path):
    """
    Read the first BED field as the template ID.

    Require at most one exact amplicon per template
    for this controlled reference panel.
    """

    if not path.exists():
        raise FileNotFoundError(path)

    ids = []

    with path.open() as handle:

        for raw_line in handle:

            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            ids.append(
                line.split("\t")[0]
            )

    if len(ids) != len(set(ids)):
        raise ValueError(
            f"Duplicate template IDs in BED: {path}"
        )

    return set(ids)


def exact_pair_ids(
    template_ids,
    forward_rows,
    reverse_rows,
):
    """Templates whose two selected primer sites are exact."""

    return {
        template_id
        for template_id in template_ids
        if (
            int(
                forward_rows[template_id][
                    "best_mismatches"
                ]
            ) == 0
            and
            int(
                reverse_rows[template_id][
                    "best_mismatches"
                ]
            ) == 0
        )
    }


def main():

    if OUTPUT.exists():
        raise FileExistsError(
            f"STOP: output already exists: {OUTPUT}"
        )

    # ==============================================================
    # 1. Core validated truth
    # ==============================================================

    core_fields, core_rows = read_tsv(
        CORE_MANIFEST
    )

    require_columns(
        core_fields,
        {
            "template_id",
            "marker",
            "template_sequence_type",
            "domain",
            "organism_id",
            "organism_name",
            "source_type",
            "source_accession",
            "source_sequence_id",
            "sequence_length",
            "sequence_sha256",
            "exact_variant_id",
            "otu97_id",
            "vsearch_cluster_id",
            "otu97_centroid_source_id",
            "selected_final",
        },
        "core manifest",
    )

    core = unique_map(
        core_rows,
        "template_id",
        "core manifest",
    )

    if len(core) != 20:
        raise ValueError(
            f"Expected 20 templates; found {len(core)}"
        )

    template_ids = set(core)

    organism_ids = {
        row["organism_id"]
        for row in core_rows
    }

    otu_ids = {
        row["otu97_id"]
        for row in core_rows
    }

    if len(organism_ids) != 16:
        raise ValueError(
            f"Expected 16 organisms; "
            f"found {len(organism_ids)}"
        )

    if len(otu_ids) != 16:
        raise ValueError(
            f"Expected 16 OTUs; "
            f"found {len(otu_ids)}"
        )

    expected_lengths = {
        template_id: int(
            row["sequence_length"]
        )
        for template_id, row in core.items()
    }

    # ==============================================================
    # 2. Taxonomy-query mapping
    # ==============================================================

    taxonomy_query_fields, taxonomy_query_rows = (
        read_tsv(TAXONOMY_QUERIES)
    )

    require_columns(
        taxonomy_query_fields,
        {
            "organism_id",
            "organism_name",
            "taxonomy_query_name",
        },
        "taxonomy query table",
    )

    taxonomy_queries = unique_map(
        taxonomy_query_rows,
        "organism_id",
        "taxonomy query table",
    )

    if set(taxonomy_queries) != organism_ids:
        raise ValueError(
            "Taxonomy-query organism IDs do not "
            "match the core manifest exactly."
        )

    # Verify that our benchmark organism labels have
    # not silently changed between files.
    for row in core_rows:

        query_row = taxonomy_queries[
            row["organism_id"]
        ]

        if (
            query_row["organism_name"]
            != row["organism_name"]
        ):
            raise ValueError(
                "Organism-name mismatch for "
                f"{row['organism_id']}"
            )

    # ==============================================================
    # 3. NCBI taxonomy snapshot
    # ==============================================================

    tax_fields, tax_rows = read_tsv(
        TAXONOMY_SUMMARY
    )

    required_tax_fields = {
        "Taxid",
        "Tax name",
        "Rank",
        "Domain/Realm name",
        "Phylum name",
        "Class name",
        "Order name",
        "Family name",
        "Genus name",
        "Species name",
    }

    require_columns(
        tax_fields,
        required_tax_fields,
        "NCBI taxonomy summary",
    )

    if len(tax_rows) != 16:
        raise ValueError(
            f"Expected 16 NCBI taxonomy rows; "
            f"found {len(tax_rows)}"
        )

    # Important:
    #
    # The downloaded "Query" column contains numeric
    # resolved TaxIDs in this snapshot, so we DO NOT
    # join using Query.
    #
    # We already inspected the output and confirmed
    # that all taxonomy_query_name values exactly match
    # the NCBI "Tax name" values.
    tax_by_name = unique_map(
        tax_rows,
        "Tax name",
        "NCBI taxonomy summary",
    )

    expected_tax_names = {
        row["taxonomy_query_name"]
        for row in taxonomy_query_rows
    }

    if set(tax_by_name) != expected_tax_names:

        missing = sorted(
            expected_tax_names - set(tax_by_name)
        )

        extra = sorted(
            set(tax_by_name) - expected_tax_names
        )

        raise ValueError(
            "NCBI taxon names do not exactly match "
            "our taxonomy queries.\n"
            f"Missing: {missing}\n"
            f"Extra: {extra}"
        )

    # Confirm that the NCBI domain agrees with the
    # domain already present in our core truth.
    for row in core_rows:

        query_name = taxonomy_queries[
            row["organism_id"]
        ]["taxonomy_query_name"]

        tax_row = tax_by_name[query_name]

        if (
            tax_row["Domain/Realm name"]
            != row["domain"]
        ):
            raise ValueError(
                f"Domain mismatch for "
                f"{row['template_id']}: "
                f"core={row['domain']!r}, "
                f"NCBI="
                f"{tax_row['Domain/Realm name']!r}"
            )

    # ==============================================================
    # 4. Physical-locus provenance
    # ==============================================================

    locus_fields, locus_rows = read_tsv(
        LOCUS_PROVENANCE
    )

    require_columns(
        locus_fields,
        {
            "template_id",
            "source_sequence_id",
            "parent_accession",
            "gene",
            "locus_tag",
            "old_locus_tag",
            "start_1based",
            "end_1based",
            "strand",
            "provenance_source",
        },
        "locus provenance",
    )

    loci = unique_map(
        locus_rows,
        "template_id",
        "locus provenance",
    )

    genome_template_ids = {
        row["template_id"]
        for row in core_rows
        if row["source_type"]
        == "genome_annotated_locus"
    }

    if set(loci) != genome_template_ids:
        raise ValueError(
            "Locus-provenance IDs do not exactly "
            "match genome-derived templates."
        )

    if len(loci) != 9:
        raise ValueError(
            f"Expected 9 genome-derived loci; "
            f"found {len(loci)}"
        )

    for template_id, locus in loci.items():

        core_row = core[template_id]

        if (
            locus["source_sequence_id"]
            != core_row["source_sequence_id"]
        ):
            raise ValueError(
                f"Source-sequence mismatch "
                f"for {template_id}"
            )

        if (
            locus["parent_accession"]
            != core_row["source_accession"]
        ):
            raise ValueError(
                f"Parent-accession mismatch "
                f"for {template_id}"
            )

    # ==============================================================
    # 5. Primer-QC evidence
    # ==============================================================

    project_forward = read_primer_side(
        PROJECT_QC,
        "forward",
        PROJECT_FORWARD,
        template_ids,
        expected_lengths,
    )

    project_reverse = read_primer_side(
        PROJECT_QC,
        "reverse",
        PROJECT_REVERSE,
        template_ids,
        expected_lengths,
    )

    # Deliberately asymmetric archaeal-specific pair:
    #
    # forward -> 50-nt search window
    # reverse -> 200-nt search window
    #
    # This avoids the false A1401R result produced
    # by the earlier 100-nt window.
    archaea_forward = read_primer_side(
        ARCHAEA_FORWARD_QC,
        "forward",
        ARCHAEA_FORWARD,
        template_ids,
        expected_lengths,
    )

    archaea_reverse = read_primer_side(
        ARCHAEA_REVERSE_QC,
        "reverse",
        ARCHAEA_REVERSE,
        template_ids,
        expected_lengths,
    )

    project_exact = exact_pair_ids(
        template_ids,
        project_forward,
        project_reverse,
    )

    archaea_exact = exact_pair_ids(
        template_ids,
        archaea_forward,
        archaea_reverse,
    )

    # Independently validate those exact-pair sets
    # against SeqKit's exact amplicon search.
    project_seqkit_exact = read_bed_ids(
        PROJECT_EXACT_BED
    )

    archaea_seqkit_exact = read_bed_ids(
        ARCHAEA_EXACT_BED
    )

    if project_exact != project_seqkit_exact:
        raise ValueError(
            "Python and SeqKit disagree for the "
            "project exact-primer-pair set."
        )

    if archaea_exact != archaea_seqkit_exact:
        raise ValueError(
            "Python and SeqKit disagree for the "
            "archaeal exact-primer-pair set."
        )

    if len(project_exact) != 8:
        raise ValueError(
            f"Expected 8 project exact pairs; "
            f"found {len(project_exact)}"
        )

    if len(archaea_exact) != 4:
        raise ValueError(
            f"Expected 4 archaeal exact pairs; "
            f"found {len(archaea_exact)}"
        )

    # ==============================================================
    # 6. Build enriched output rows
    # ==============================================================

    enriched_rows = []

    for core_row in core_rows:

        template_id = core_row["template_id"]

        query_name = taxonomy_queries[
            core_row["organism_id"]
        ]["taxonomy_query_name"]

        tax = tax_by_name[query_name]

        locus = loci.get(template_id)

        pf = project_forward[template_id]
        pr = project_reverse[template_id]

        af = archaea_forward[template_id]
        ar = archaea_reverse[template_id]

        output_row = dict(core_row)

        # ----------------------------------------------------------
        # Taxonomy
        # ----------------------------------------------------------

        output_row.update(
            {
                "taxonomy_query_name":
                    query_name,

                "ncbi_taxid":
                    tax["Taxid"],

                "ncbi_tax_name":
                    tax["Tax name"],

                "ncbi_tax_rank":
                    tax["Rank"],

                "ncbi_domain":
                    tax["Domain/Realm name"],

                "ncbi_phylum":
                    tax["Phylum name"],

                "ncbi_class":
                    tax["Class name"],

                "ncbi_order":
                    tax["Order name"],

                "ncbi_family":
                    tax["Family name"],

                "ncbi_genus":
                    tax["Genus name"],

                "ncbi_species":
                    tax["Species name"],
            }
        )

        # ----------------------------------------------------------
        # Physical-locus provenance
        # ----------------------------------------------------------

        if locus is None:

            locus_values = {
                "physical_locus_provenance_available":
                    "no",
                "locus_parent_accession": ".",
                "locus_gene": ".",
                "locus_tag": ".",
                "locus_old_locus_tag": ".",
                "locus_start_1based": ".",
                "locus_end_1based": ".",
                "locus_strand": ".",
                "locus_provenance_source": ".",
            }

        else:

            locus_values = {
                "physical_locus_provenance_available":
                    "yes",

                "locus_parent_accession":
                    locus["parent_accession"],

                "locus_gene":
                    locus["gene"],

                "locus_tag":
                    locus["locus_tag"],

                "locus_old_locus_tag":
                    locus["old_locus_tag"],

                "locus_start_1based":
                    locus["start_1based"],

                "locus_end_1based":
                    locus["end_1based"],

                "locus_strand":
                    locus["strand"],

                "locus_provenance_source":
                    locus["provenance_source"],
            }

        output_row.update(locus_values)

        # ----------------------------------------------------------
        # Project-primer summary
        # ----------------------------------------------------------

        output_row.update(
            {
                "project_primer_terminal_window":
                    "100",

                "project_forward_best_mismatches":
                    pf["best_mismatches"],

                "project_forward_3prime_mismatches":
                    pf["three_prime_mismatches"],

                "project_forward_distance_from_expected_end":
                    pf["distance_from_expected_end"],

                "project_reverse_best_mismatches":
                    pr["best_mismatches"],

                "project_reverse_3prime_mismatches":
                    pr["three_prime_mismatches"],

                "project_reverse_distance_from_expected_end":
                    pr["distance_from_expected_end"],

                "project_exact_pair":
                    (
                        "yes"
                        if template_id in project_exact
                        else "no"
                    ),
            }
        )

        # ----------------------------------------------------------
        # Archaeal-specific-primer summary
        # ----------------------------------------------------------

        output_row.update(
            {
                "archaea_specific_forward_terminal_window":
                    "50",

                "archaea_specific_forward_best_mismatches":
                    af["best_mismatches"],

                "archaea_specific_forward_3prime_mismatches":
                    af["three_prime_mismatches"],

                "archaea_specific_forward_distance_from_expected_end":
                    af["distance_from_expected_end"],

                "archaea_specific_reverse_terminal_window":
                    "200",

                "archaea_specific_reverse_best_mismatches":
                    ar["best_mismatches"],

                "archaea_specific_reverse_3prime_mismatches":
                    ar["three_prime_mismatches"],

                "archaea_specific_reverse_distance_from_expected_end":
                    ar["distance_from_expected_end"],

                "archaea_specific_exact_pair":
                    (
                        "yes"
                        if template_id in archaea_exact
                        else "no"
                    ),
            }
        )

        enriched_rows.append(output_row)

    # ==============================================================
    # 7. Write without overwriting the validated core manifest
    # ==============================================================

    fieldnames = list(
        enriched_rows[0].keys()
    )

    with OUTPUT.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(enriched_rows)

    print(
        f"Wrote {len(enriched_rows)} enriched templates."
    )

    print(
        f"Organisms: {len(organism_ids)}"
    )

    print(
        f"97% OTUs: {len(otu_ids)}"
    )

    print(
        f"Genome-derived loci: {len(loci)}"
    )

    print(
        f"Project exact primer pairs: "
        f"{len(project_exact)}"
    )

    print(
        f"Archaeal exact primer pairs: "
        f"{len(archaea_exact)}"
    )

    print(
        f"Columns: {len(fieldnames)}"
    )

    print(
        f"Output: {OUTPUT}"
    )


if __name__ == "__main__":
    main()
