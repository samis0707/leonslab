# Curated gene+flank reference records

Per-`(gene, isolate)` records produced by `scripts/build_gene_records.py`.
See `docs/primer_design/primer_design_skill_v2.md` §1.3 for the canonical
header format and sequence layout.

## Coverage convention — when an isolate is missing

For each gene there is one `<isolate>.fasta` per isolate where the gene
is present **and functional**. Missing files carry meaning:

| Gene | Convention for missing file |
|---|---|
| `lasB` | gene not present in this isolate's curated CDS list (assembly gap or true deletion) |
| `lasR` | **lasR is non-functional at the amino-acid level in this isolate** — explicit biological annotation by Leon. The tool refuses to design lasR primers for these isolates because the FASTA is absent, which matches the intent (no point tagging or expressing a broken gene from its native locus) |
| `lasI` | gene FASTA not yet uploaded to R2 (added later — see project_plan §6 step 1 backfill note) |

The tool surfaces missing-record errors via `gene_finder.GeneRecordNotFound`
with the isolate ID and gene name; downstream UIs can interpret these per
the table above.
