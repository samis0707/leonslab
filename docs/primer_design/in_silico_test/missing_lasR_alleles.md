# Non-functional and absent lasR alleles

The curated `lasR.fasta` master file in R2 only contains records for isolates
whose lasR allele is functional at the protein level (canonical CDS: starts
ATG, ends with a single stop, mod-3 length, no internal stops). Isolates whose
lasR is truncated or absent were excluded.

`scripts/extract_truncated_lasR.py` recovers these alleles from the isolate
genomes by anchor-based homology search, classifies the truncation reason,
and writes per-isolate FASTAs with header field `functional=false` plus a
`truncation_reason=…` annotation. The primer designer reads this flag,
substitutes the reference CDS for scar logic, and surfaces a prominent
WARNING in the PDF / JSON output so the user knows the native protein product
is already absent or non-functional in this isolate.

## Round 1 results (16 isolates flagged as missing lasR)

### ✅ 13 successfully extracted as `functional=false`

| Isolate | Reason                              | Locus length |
|---------|-------------------------------------|-------------:|
| LB001   | premature_stop@codon_180            | 728 nt       |
| LB009   | insertion_in_locus_~1196nt          | 1916 nt      |
| LB020   | premature_stop@codon_232            | 718 nt       |
| LB038   | premature_stop@codon_114            | 719 nt       |
| LB039   | premature_stop@codon_114            | 719 nt       |
| LB040   | premature_stop@codon_114            | 719 nt       |
| LB042   | premature_stop@codon_114            | 719 nt       |
| LB043   | premature_stop@codon_114            | 719 nt       |
| LB052   | insertion_in_locus_~918nt           | 1638 nt      |
| LB069   | insertion_in_locus_~1116nt          | 1836 nt      |
| LB086   | premature_stop@codon_121            | 740 nt       |
| LB091   | insertion_in_locus_~1427nt          | 2147 nt      |
| LB117   | premature_stop@codon_114            | 719 nt       |

The 6 isolates that share a premature stop at codon 114 likely carry the same
ancestral mutation (a single SNP turning a Q codon into TAA/TAG/TGA) — worth
phylogenetic investigation when convenient.

### ❌ 3 isolates not designable: locus completely absent

LB045, LB046, LB076 have **no detectable lasR homology** anywhere in their
genomes — not just CDS but also the ~600 nt UP/DN flanks. Confirmed not a
species mismatch (lasB anchors hit cleanly in all three). Most likely a large
chromosomal deletion (>1.5 kb) that removed lasR plus surrounding context.

For these isolates ΔlasR primer design is not biologically meaningful — there
is no locus to delete. The tool will return `GeneRecordNotFound` if the user
selects this combination. The dropdown UI should ideally hide the "delete
lasR" option for these three isolates; for now the error message is the
guard.

## Re-running the extraction

```bash
set -a && source api/design/.env && set +a
python3 scripts/extract_truncated_lasR.py
```

The script is idempotent (overwrites existing files). Add new isolates by
editing the `MISSING_ISOLATES` list at the top of the script.
