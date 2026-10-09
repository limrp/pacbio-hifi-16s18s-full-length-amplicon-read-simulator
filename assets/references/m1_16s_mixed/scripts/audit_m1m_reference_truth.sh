#!/usr/bin/env bash

set -euo pipefail

BASE="data/m1_work/mixed_16s"
TRUTH_BUILD="$BASE/final_truth_build"
QC="$BASE/final_panel_qc"
PRIMER_QC="$TRUTH_BUILD/primer_compatibility"

MANIFEST="$TRUTH_BUILD/m1m_truth_manifest_enriched.tsv"
CORE="$TRUTH_BUILD/m1m_truth_manifest.tsv"
FASTA="$TRUTH_BUILD/m1m_16s_templates.fasta"
HASHES="$TRUTH_BUILD/m1m_16s_templates_sha256.tsv"

FINAL_UC="$QC/clusters_97.uc"
STABILITY="$QC/order_stability_sha256.tsv"

SENS_FULL="$PRIMER_QC/sensitivity_full_length_97.uc"
SENS_AMP="$PRIMER_QC/sensitivity_amplicon_97.uc"


fail() {
    printf 'AUDIT FAIL: %s\n' "$1" >&2
    exit 1
}


printf '=== M1M reference-truth final audit ===\n'


# ------------------------------------------------------------
# Enriched-manifest structure
# ------------------------------------------------------------

rows=$(
    awk 'END {print NR - 1}' "$MANIFEST"
)

cols=$(
    awk -F '\t' 'NR == 1 {print NF; exit}' "$MANIFEST"
)

[[ "$rows" -eq 20 ]] \
    || fail "Expected 20 manifest rows; found $rows"

[[ "$cols" -eq 53 ]] \
    || fail "Expected 53 columns; found $cols"


# ------------------------------------------------------------
# Core biological truth
# ------------------------------------------------------------

organisms=$(
    awk -F '\t' '
    NR > 1 {
        x[$5] = 1
    }
    END {
        print length(x)
    }
    ' "$MANIFEST"
)

variants=$(
    awk -F '\t' '
    NR > 1 {
        x[$12] = 1
    }
    END {
        print length(x)
    }
    ' "$MANIFEST"
)

otus=$(
    awk -F '\t' '
    NR > 1 {
        x[$13] = 1
    }
    END {
        print length(x)
    }
    ' "$MANIFEST"
)

sequence_hashes=$(
    awk -F '\t' '
    NR > 1 {
        x[$11] = 1
    }
    END {
        print length(x)
    }
    ' "$MANIFEST"
)

[[ "$organisms" -eq 16 ]] \
    || fail "Expected 16 organisms; found $organisms"

[[ "$variants" -eq 20 ]] \
    || fail "Expected 20 exact variants; found $variants"

[[ "$otus" -eq 16 ]] \
    || fail "Expected 16 OTUs; found $otus"

[[ "$sequence_hashes" -eq 20 ]] \
    || fail "Expected 20 unique sequence hashes; found $sequence_hashes"


# ------------------------------------------------------------
# Domain structure
# ------------------------------------------------------------

archaea=$(
    awk -F '\t' '
    NR > 1 && $4 == "Archaea" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

bacteria=$(
    awk -F '\t' '
    NR > 1 && $4 == "Bacteria" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

[[ "$archaea" -eq 12 ]] \
    || fail "Expected 12 archaeal templates; found $archaea"

[[ "$bacteria" -eq 8 ]] \
    || fail "Expected 8 bacterial templates; found $bacteria"


# ------------------------------------------------------------
# Selected-final flag
# ------------------------------------------------------------

bad_selected=$(
    awk -F '\t' '
    NR > 1 && $16 != "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

[[ "$bad_selected" -eq 0 ]] \
    || fail "Some templates are not selected_final=yes"


# ------------------------------------------------------------
# Physical-locus provenance
# ------------------------------------------------------------

loci=$(
    awk -F '\t' '
    NR > 1 && $28 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

[[ "$loci" -eq 9 ]] \
    || fail "Expected 9 genome-derived loci; found $loci"


# ------------------------------------------------------------
# Primer-pair truth
# ------------------------------------------------------------

project_exact=$(
    awk -F '\t' '
    NR > 1 && $44 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

archaea_exact=$(
    awk -F '\t' '
    NR > 1 && $53 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

[[ "$project_exact" -eq 8 ]] \
    || fail "Expected 8 project exact pairs; found $project_exact"

[[ "$archaea_exact" -eq 4 ]] \
    || fail "Expected 4 archaeal exact pairs; found $archaea_exact"


project_archaea=$(
    awk -F '\t' '
    NR > 1 && $4 == "Archaea" && $44 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

project_bacteria=$(
    awk -F '\t' '
    NR > 1 && $4 == "Bacteria" && $44 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

archpair_archaea=$(
    awk -F '\t' '
    NR > 1 && $4 == "Archaea" && $53 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

archpair_bacteria=$(
    awk -F '\t' '
    NR > 1 && $4 == "Bacteria" && $53 == "yes" {
        n++
    }
    END {
        print n + 0
    }
    ' "$MANIFEST"
)

[[ "$project_archaea" -eq 0 ]] || fail "Unexpected exact project pair in Archaea"
[[ "$project_bacteria" -eq 8 ]] || fail "Expected 8 exact bacterial project pairs"
[[ "$archpair_archaea" -eq 4 ]] || fail "Expected 4 exact archaeal-specific pairs in Archaea"
[[ "$archpair_bacteria" -eq 0 ]] || fail "Unexpected exact archaeal-specific pair in Bacteria"


arch_ids=$(
    awk -F '\t' '
    NR > 1 && $53 == "yes" {
        print $1
    }
    ' "$MANIFEST" \
    | sort \
    | paste -sd, -
)

expected_arch_ids="M1M_16S_T004,M1M_16S_T005,M1M_16S_T007,M1M_16S_T008"

[[ "$arch_ids" == "$expected_arch_ids" ]] \
    || fail "Unexpected exact archaeal-primer template set: $arch_ids"


# ------------------------------------------------------------
# Stable FASTA
# ------------------------------------------------------------

read -r fasta_n fasta_sum fasta_min fasta_max fasta_ambiguous < <(
    seqkit fx2tab -i "$FASTA" \
    | awk -F '\t' '
    {
        sequence = toupper($2)
        len = length(sequence)

        n++
        total += len

        if (n == 1 || len < min)
            min = len

        if (len > max)
            max = len

        non_acgt = sequence
        gsub(/[ACGT]/, "", non_acgt)

        ambiguous += length(non_acgt)
    }

    END {
        print n, total, min, max, ambiguous
    }
    '
)

[[ "$fasta_n" -eq 20 ]] || fail "FASTA does not contain 20 templates"
[[ "$fasta_sum" -eq 29988 ]] || fail "FASTA total length changed"
[[ "$fasta_min" -eq 1436 ]] || fail "FASTA minimum length changed"
[[ "$fasta_max" -eq 1561 ]] || fail "FASTA maximum length changed"
[[ "$fasta_ambiguous" -eq 0 ]] || fail "FASTA contains ambiguous bases"


# ------------------------------------------------------------
# Enrichment must not alter original truth
# ------------------------------------------------------------

diff -q \
    "$CORE" \
    <(cut -f1-16 "$MANIFEST") \
    >/dev/null \
    || fail "Core 16-column truth changed during enrichment"


diff -q \
    "$HASHES" \
    <(cut -f1,9-11 "$MANIFEST") \
    >/dev/null \
    || fail "Stable-template hashes disagree with enriched manifest"


# ------------------------------------------------------------
# Final 97% clustering
# ------------------------------------------------------------

clusters=$(
    grep -c '^C' "$FINAL_UC"
)

[[ "$clusters" -eq 16 ]] \
    || fail "Expected 16 final 97% clusters; found $clusters"


# ------------------------------------------------------------
# Order-stability evidence
# ------------------------------------------------------------

stability_runs=$(
    wc -l < "$STABILITY"
)

distinct_stability_hashes=$(
    awk '{print $1}' "$STABILITY" \
    | sort -u \
    | wc -l
)

[[ "$stability_runs" -eq 8 ]] \
    || fail "Expected 8 order-stability runs; found $stability_runs"

[[ "$distinct_stability_hashes" -eq 1 ]] \
    || fail "Order-stability outputs are not byte-identical"


# ------------------------------------------------------------
# Primer-boundary sensitivity
# ------------------------------------------------------------

full_clusters=$(
    grep -c '^C' "$SENS_FULL"
)

full_hits=$(
    grep -c '^H' "$SENS_FULL" || true
)

amp_clusters=$(
    grep -c '^C' "$SENS_AMP"
)

amp_hits=$(
    grep -c '^H' "$SENS_AMP" || true
)

[[ "$full_clusters" -eq 4 ]] \
    || fail "Full-length sensitivity subset should have 4 OTUs"

[[ "$full_hits" -eq 0 ]] \
    || fail "Full-length sensitivity subset unexpectedly contains H records"

[[ "$amp_clusters" -eq 4 ]] \
    || fail "Amplicon sensitivity subset should have 4 OTUs"

[[ "$amp_hits" -eq 0 ]] \
    || fail "Amplicon sensitivity subset unexpectedly contains H records"


# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

printf '\n'
printf 'Templates:                     %s\n' "$rows"
printf 'Archaea:                       %s\n' "$archaea"
printf 'Bacteria:                      %s\n' "$bacteria"
printf 'Organisms:                     %s\n' "$organisms"
printf 'Exact variants:                %s\n' "$variants"
printf '97%% OTUs:                      %s\n' "$otus"
printf 'Genome-derived loci:           %s\n' "$loci"
printf 'Project exact primer pairs:    %s\n' "$project_exact"
printf 'Archaeal exact primer pairs:   %s\n' "$archaea_exact"
printf 'Order-stability runs:          %s\n' "$stability_runs"
printf 'Distinct stability hashes:     %s\n' "$distinct_stability_hashes"
printf 'Ambiguous FASTA bases:         %s\n' "$fasta_ambiguous"
printf '\nAUDIT PASS\n'
