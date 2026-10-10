# M1 R01 Final Freeze — Mixed 16S Controlled Simulation

Freeze date: 2026-10-08

## Scope

This freeze records the validated R01 outputs for the M1 mixed-16S
controlled PacBio HiFi/MHASS simulation.

Dataset metadata:

- marker: 16S
- dataset_id: toyctrl16
- simulation_replicate: R01
- master_seed: 16001
- truth templates: 20
- samples: S01, S02, S03
- intended reads per sample: 5120
- intended reads total: 15360
- realized reads total: 10228

Realized counts:

- S01: 3415
- S02: 3397
- S03: 3416

## MHASS-flank preparation policy

Fixed [35:-35] trimming was rejected because simulator/library-wrapper
indels shift the observed biological boundaries in many realized reads.

The final preparation method uses provenance-informed global alignment of
the observed read against the expected complete MHASS construct and
projects the known biological reference interval onto observed-query
coordinates using the validated CIGAR walker.

For boundary insertions:

- left boundary: use alignment_target_start_min
- right boundary: use alignment_target_end_max

Therefore boundary insertions are retained rather than silently deleted.

The observed read sequence is never replaced by the biological truth
sequence.

For reverse-oriented reads:

- sequence is reverse-complemented
- Phred quality scores are reversed

Sequence and quality strings are sliced using identical coordinates.

## Validation status

Read-level validation:

- raw reads: 10228
- diagnostic rows: 10228
- truth_reads.tsv rows: 10228
- prepared FASTQ reads: 10228
- total validation failures: 0
- validation status: PASS
- validator exit status: 0

Simulation-level audit:

- truth templates: 20
- intended molecules: 15360
- MHASS mapping rows: 15360
- realized reads: 10228
- zero-intended template/sample pairs: 7
- zero-intended pairs producing reads: 0
- total audit failures: 0
- simulation audit status: PASS
- validator exit status: 0

## Authoritative validation manifests

Frozen simulation design:

assets/simulation/m1_16s_mixed/SHA256SUMS

Raw R01 simulation:

data/m1_work/mixed_16s/r01_postrun_qc/R01_SHA256SUMS

Prepared R01 dataset:

data/m1_work/mixed_16s/r01_postrun_qc/prepared_full_R01_v01_SHA256SUMS

Complete authoritative method code:

data/m1_work/mixed_16s/r01_postrun_qc/M1_R01_METHOD_CODE_SHA256SUMS

Read-level validation evidence:

data/m1_work/mixed_16s/r01_postrun_qc/prepared_full_R01_v01_FINAL_EVIDENCE_SHA256SUMS

Simulation-level audit code:

data/m1_work/mixed_16s/r01_postrun_qc/test_simulation_R01_AUDIT_CODE_SHA256SUMS

Simulation-level audit evidence:

data/m1_work/mixed_16s/r01_postrun_qc/test_simulation_R01_AUDIT_EVIDENCE_SHA256SUMS

## Freeze rule

The artifacts referenced by these manifests are considered immutable for
M1 R01.

Do not silently modify or overwrite them after this freeze.

If a scientifically relevant change becomes necessary, create a new
versioned output/freeze rather than replacing this validated R01 package.

Smoke-test and exploratory diagnostic files may remain in the working
directory, but they are not the authoritative final R01 outputs unless
explicitly referenced by the manifests above.
