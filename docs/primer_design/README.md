# Primer Design Tool — Code Skeleton

Phase-2 starting point for the *P. aeruginosa* primer design tool. Read these documents in order before implementing:

1. `docs/primer_design/decisions_log.md` — fixed decisions, single source of truth
2. `docs/primer_design/project_plan.md` — architecture, modules, data model, roadmap
3. `docs/primer_design/primer_design_skill_v2.md` — algorithmic specification, pseudocode, worked examples

## Where things live

```
api/design/handler.py            Vercel serverless function entry point (POST /api/design)
lib/primer_design/               Core algorithm package
├── config.py                    constants, Tm parameters, polymerase offsets
├── types.py                     dataclasses (GeneRecord, PrimerSet, DesignRequest, ...)
├── exceptions.py                typed errors surfaced as JSON to the frontend
├── tags.py                      tag library (FLAG, His6, etc.) — IMPLEMENTED
├── vectors.py                   vector loading + cut detection + tail rules — PARTIAL
├── gene_finder.py               canonical record retrieval — IMPLEMENTED
├── storage_adapter.py           R2 access for full genomes — IMPLEMENTED
├── primer_designer.py           filters, Tm, scoring, search — PARTIAL
├── plasmid_builder.py           final assembly — STUB
├── verification.py              hard checks per D7.4 — PARTIAL
├── fasta_writer.py              combined FASTA output — STUB
├── pdf_writer.py                3 PDF templates — STUB
├── json_writer.py               manifest output — STUB
└── applications/
    ├── deletion.py              orchestration for pEXG2 deletion — STUB
    ├── tagging.py               orchestration for pEXG2 in-locus C-term — STUB
    └── expression.py            orchestration for pBBR1MCS2 expression — STUB

data/primer_design/
├── tail_conventions.json        per-(vector, enzyme, application) calibration — PROVIDED
├── codon_table_paeruginosa.json — PROVIDED
├── manifest.json                (placeholder — populated by upload_genomes.py)
├── genes/<gene>/<isolate>.fasta (populated by build_gene_records.py)
└── vectors/{pEXG2,pBBR1MCS2}.gb (populated in step 2 of roadmap)

scripts/
├── build_gene_records.py        local pipeline: fetches gene FASTAs + genomes from R2,
│                                  extracts flanks, writes canonical per-(gene, isolate) records
├── upload_genomes.py            one-time R2 upload (already done by user)
└── verify_manifest.py           SHA-256 verification of bundled vs R2 genomes

tests/
├── unit/                        per-module unit tests — STUB SIGNATURES
├── integration/                 end-to-end tests for 3 worked examples — STUB
├── golden/                      regression artefacts — placeholder
└── fixtures/                    mini test data — placeholder
```

## Module status legend

- **IMPLEMENTED** — Phase-2 should review and run; minimal further work expected
- **PARTIAL** — core algorithm written, edge cases / orchestration to fill in
- **STUB** — function signatures + docstrings + acceptance criteria; body to be written

## Getting started in Phase 2

Follow `project_plan.md` §6 roadmap. Step ordering matters because later modules depend on earlier ones.

```
Step 1 — build curated gene+flank records (run scripts/build_gene_records.py)
Step 2 — annotate vectors as GenBank, drop into data/primer_design/vectors/
Step 3 — calibrate tail conventions (verify against unit tests in tests/unit/test_vectors.py)
Step 4 — verify R2 manifest (run scripts/verify_manifest.py)
Step 5 — make tests/unit/ green
Step 6 — make tests/integration/test_deletion_LB001_lasB.py green
Step 7 — make tests/integration/test_expression_PA14_lasR.py green
Step 8 — make tests/integration/test_tagging_PA14_lasR_His6.py green
Step 9 — flesh out api/design/handler.py + frontend integration
Step 10 — deploy to Vercel preview
```

Acceptance gates per step are listed in `project_plan.md` §6.

## R2 bucket layout (verified 2026-05-04)

```
leonslab-pa-genomes/
├── Whole genome sequences/      ← all .fna files (case-sensitive, contains spaces)
│   ├── LB001.fna
│   ├── LB014.fna
│   └── ...
├── lasB.fasta                   ← gene CDS list, multi-record (one per isolate)
├── lasR.fasta
└── lasI.fasta                   (added later)
```

The `storage_adapter.py` knows about this layout via `R2_GENOMES_PREFIX = "Whole genome sequences/"`.

## Required environment variables

Set in Vercel for Production + Preview + Development:

```
R2_ACCOUNT_ID
R2_ACCESS_KEY_ID
R2_SECRET_ACCESS_KEY
R2_BUCKET_NAME      # = "leonslab-pa-genomes"
```

For local dev: `vercel env pull .env.local` then `python-dotenv` reads it.

## Conventions

- Python 3.11+ (matches Vercel current Python runtime)
- Type hints required on public APIs
- Docstrings: Google-style
- Lint: ruff (default ruleset)
- Format: black, line length 100
- Tests: pytest, no parallel runners
- All errors are typed exceptions; surface to handler as JSON error responses
