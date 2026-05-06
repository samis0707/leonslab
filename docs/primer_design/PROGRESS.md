# Phase 2 — Progress (handoff)

Living checkpoint for the *P. aeruginosa* primer design tool. Read this
*after* the three SSoT docs (`decisions_log.md`, `project_plan.md`,
`primer_design_skill_v2.md`) and the README in this folder.

> Tip for new chats: clone the repo, run `pytest tests/unit/`. Baseline as
> of 2026-05-06 is **56 passed, 1 skipped**; integration tests in
> `tests/integration/` add another 19 passed when `api/design/.env`
> carries the R2 credentials.

---

## Where we are

| # | Roadmap step (project_plan §6) | Status | Acceptance gate |
|---|---|---|---|
| 1 | Build curated gene+flank records | **DONE** | 110 lasB / 94 lasR / 2 lasI records; PA14 + PAO1 references added via `scripts/build_reference_records.py` |
| 2 | Annotate vectors as GenBank | **DONE** | `pEXG2.gb` and `pBBR1MCS2.gb` in `data/primer_design/vectors/`. All 4 tail-derivation tests green |
| 3 | Calibrate tail conventions | **DONE** | All 4 empirical tails reproduced |
| 4 | R2 manifest | **DONE** | 110 isolates incl. PA14/PAO1 with SHA-256 in `data/primer_design/manifest.json` |
| 5 | gene_finder + storage_adapter | **DONE** | |
| 6 | primer_designer.py core | **DONE** | `search_deletion_primers` byte-exact for LB001 lasB; off_target_scan via anchor-mismatch index |
| 7 | applications/deletion.py + integration | **DONE** | `tests/integration/test_deletion_LB001_lasB.py` 8/8 passed |
| 8 | applications/expression.py + integration | **DONE** | `tests/integration/test_expression_PA14_lasR.py` 6/6 passed; tiered relaxation for anchored primers (clamp last-resort) |
| 9 | applications/tagging.py + integration | **DONE** | `tests/integration/test_tagging_PA14_lasR_His6.py` 5/5 passed (His6 byte-exact + 3xFLAG/HiBiT rejection) |
| 10 | plasmid_builder + verification full | **DONE** | `count_recognition_sites` factored out; 3 previously-skipped tests now green |
| 11 | pdf_writer (3 templates) | **DONE** | A4 single-page deletion + tagging + expression PDFs in `examples/`. B7 polymerase only per user request |
| 12 | handler.py + frontend wiring | **DONE** | `api/design/handler.py` POST endpoint (already complete); Astro form at `/primer-design/` with build-time catalog, cascading dropdowns, fetch + result display |
| 13 | CI | **DONE** | `.github/workflows/ci.yml` runs Python unit tests + Astro build on push/PR |

## Test counts

```
pytest tests/unit/         →  56 passed, 1 skipped   (no R2 needed)
pytest tests/integration/  →  19 passed              (R2 + .env required)
```

The single remaining skip is an unrelated primer-dimer edge case in
`test_primer_filters.py::TestMax3PrimeSelfDimer::test_palindromic_3prime_high`.

## Algorithm validation snapshot (LB001 ΔlasB pEXG2 HindIII)

This is Step 6's empirical anchor — produced by `search_deletion_primers`:

```
N=16, C=2, scar = MKKVSTLDLLFVAIMGAL*
P1: CATAAATGTAAAGCA + GGTGTTCCAGCTGGTGCAG    (Tm 60.9 °C)
P2: CCGAGCTTACAACGC + ACCCATGATCGCAACGAACAAC (Tm 60.8 °C)
P3: GTTGCGATCATGGGT + GCGTTGTAAGCTCGGTGGTC   (Tm 61.0 °C)
P4: CGACCTGCAGAAGCT + GCCAGGTACTCGCCTTGC     (Tm 60.7 °C)
Tm spread: 0.21 °C, score 0.525, 117 ms
```

All four bodies match skill_v2 §4.4 expected output **exactly** (no ±2 nt
shift needed).

## Step 6 carry-over — what's still stubbed in primer_designer.py

- `search_tagging_primers`  (planned for Step 9)
- `search_expression_primers` (planned for Step 8)
- `off_target_scan` (needed by all three application orchestrators in Steps 7-9)

`off_target_scan` was originally listed as Step 6 but the integration tests
that actually need it live with the application orchestrators, so it's
deferred to whichever orchestrator hits it first (likely Step 7 / deletion).

## Local dev quick-start

```bash
git pull
pip install -r requirements.txt          # biopython, boto3, reportlab, pytest, …
cp api/design/.env.example .env          # then fill in R2_* values
                                          # (or `vercel env pull .env.local`)
pytest tests/unit/ -v                    # baseline: 53 passed, 4 skipped
```

To re-run Step 1 against R2 (after lasI lands, or when adding isolates):

```bash
set -a; source .env; set +a
python scripts/build_gene_records.py --genes lasI       # backfill missing gene
# or
python scripts/build_gene_records.py --isolates LB200   # backfill new isolate
```

## Conventions worth remembering

- `data/primer_design/genes/<gene>/<isolate>.fasta`: missing file ≢ "not yet built". For `lasR` specifically, missing means the isolate has a non-functional lasR at the amino-acid level (Leon's biological annotation). See `data/primer_design/genes/README.md`.
- `cds_seq` always includes the trailing stop codon. `len(cds) % 3 == 0`. Codon count `L = len(cds) // 3` is total codons *including* stop.
- `scar_orf` formula: `cds[:3*N] + cds[3*(L - C - 1):]`. The `-1` matters; the skill_v2 §4.1 pseudocode has an off-by-one but the worked example confirms this corrected form.
- Tail conventions in JSON name two rules per `(vector, enzyme, application)` triple. The rule name fixes the *physical span*, not the P1-vs-P4 role; the `p1_tail_rule` / `p4_tail_rule` keys decide which side gets which rule. See the `_geometry_note` field on the pBBR1MCS2 entry for the orientation reasoning.

## Next session pickup

- **Step 7** (deletion app + integration): wire `applications/deletion.py`. Needs:
  1. `off_target_scan` real implementation (skill_v2 §10)
  2. Plasmid assembly (`plasmid_builder.assemble`) at least to the point where the integration test can verify the final plasmid SHA-256 — see expected `final_plasmid_length_bp = 6155` in skill_v2 §15.1
  3. `verification.verify` for the deletion-specific gates (D7.4)
- The integration test stub already lists the expected primer values in `tests/integration/test_deletion_LB001_lasB.py`; the search core now reproduces them, so that file's `EXPECTED` dict is the spec.
