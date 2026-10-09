# M1 mixed 16S controlled simulation inputs

Frozen MHASS inputs for the mixed Bacteria + Archaea M1 controlled benchmark.

## Samples

- S01: balanced organism-level community
- S02: uneven abundance gradient with explicit absence
- S03: sparse E. coli / Salmonella 97%-OTU stress test

Each sample has 5,120 requested reads.

## Template truth

- 20 exact full-length 16S templates
- 12 Archaea
- 8 Bacteria
- 16 organisms
- 16 operational 97% OTUs

Counts are defined at the organism level first and then distributed among
the selected exact variants.

The fixed count matrix is hand-authored experimental truth.

## S03 stress test

E. coli:
- T016 -> M1M_OTU97_014
- T017 -> M1M_OTU97_015

Salmonella:
- T018 -> M1M_OTU97_014

For S03:
- M1M_OTU97_014 = 2,560 requested reads = 50%
- M1M_OTU97_015 = 1,280 requested reads = 25%
- combined E. coli + Salmonella = 3,840 requested reads = 75%

Thus organism abundance truth is deliberately not identical to
97%-OTU abundance truth.

## Simulation replicate

The same frozen truth_sequences.fasta and fixed_counts.tsv should be reused
unchanged across controlled 16S simulation replicates.

Planned seeds:
- R01 = 16001
- R02 = 16002
- R03 = 16003

Requested counts are intended template truth.
Realized MHASS reads must be recorded separately as read-level truth.
