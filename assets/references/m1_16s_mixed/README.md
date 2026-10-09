# M1 mixed full-length 16S reference truth

Frozen reference package for the controlled M1 benchmark.

## Reference design

- Marker: 16S rRNA
- Template definition: biological full-length 16S sequence
- Domains: Archaea + Bacteria
- Templates: 20
- Archaea templates: 12
- Bacteria templates: 8
- Biological organisms represented: 16
- Exact nucleotide variants: 20
- VSEARCH operational OTUs at 97% identity: 16
- Genome-derived annotated loci: 9
- Isolated 16S reference records: 11

## Operational 97% truth

VSEARCH:
- cluster_fast
- identity threshold: 0.97
- iddef: 2

The final 20-template panel produces:
- 16 OTUs
- 12 singleton OTUs
- 4 two-template OTUs

Eight tested FASTA input orders produced byte-identical UC outputs.

## Primer metadata

Primer compatibility is metadata describing the biological templates.
It does not define or alter the M1 simulated full-length sequences.

Project 16S pair:
- exact complete pairs: 8/20
- all 8 occur in the selected bacterial templates

Exploratory archaeal-specific pair:
- exact complete pairs: 4/20
- exact templates: M1M_16S_T004, M1M_16S_T005,
  M1M_16S_T007, M1M_16S_T008

For the archaeal-specific pair:
- forward-primer QC uses a 50-nt terminal search window
- reverse-primer QC uses a 200-nt terminal search window

## Primer-boundary sensitivity

The four exact archaeal-primer-compatible templates were compared as:

1. full-length biological 16S sequences
2. exact primer-bounded amplicons

Both definitions yielded four singleton 97% OTUs.

Thus primer bounding altered sequence length and some pairwise identities,
but did not alter the operational 97%-OTU interpretation in this tested subset.

## Important scope

This controlled benchmark evaluates sequence-inference recovery from known
biological 16S truth. It does not claim to model wet-lab PCR amplification
efficiency.

PCR primers are not synthetically appended to the biological reference
templates.
