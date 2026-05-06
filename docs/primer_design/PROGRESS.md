# Phase 2 — Progress (handoff)

Living checkpoint for the *P. aeruginosa* primer design tool. Read this
*after* the three SSoT docs (`decisions_log.md`, `project_plan.md`,
`primer_design_skill_v2.md`) and the README in this folder.

> Tip for new chats: clone the repo, run `pytest tests/unit/`. If it's not
> 53 passed / 4 skipped (or higher), there's been drift since this note.

---

## Where we are

| # | Roadmap step (project_plan §6) | Status | Acceptance gate |
|---|---|---|---|
| 1 | Build curated gene+flank records | **PARTIAL** | 200 / ~360 records — `lasB` (108) + `lasR` (92) done; `lasI` deferred until that FASTA lands in R2 |
| 2 | Annotate vectors as GenBank | **DONE** | `pEXG2.fasta` from ENA `KM887143.1`; `pBBR1MCS2.gb` from Addgene #85168. `get_unique_cutters(pBBR1MCS2)` = 13 (D4.2). All 4 tail-derivation tests green |
| 3 | Calibrate tail conventions | **DONE** | All 4 empirical tails reproduced. Two bug fixes baked in: off-by-one in `_rule_left_arm_15nt_incl_5nt_recognition`, and P1/P4 swap in `tail_conventions.json` for `pBBR1MCS2|HindIII|expression` (Addgene's MCS reads `…SalI-ClaI-HindIII-EcoRV-EcoRI…`) |
| 4 | R2 manifest | **PARTIAL** | 108 isolates with SHA-256 in `data/primer_design/manifest.json`. Full panel will accumulate as more isolates touch `build_gene_records.py` |
| 5 | gene_finder + storage_adapter | **DONE** (Phase 1 implementation) | All gene_finder unit tests green |
| **6** | **primer_designer.py core** | **DONE** | `search_deletion_primers` reproduces the v1 empirical set for LB001 ΔlasB pEXG2 HindIII byte-exact (N=16, C=2, all 4 primer tails AND bodies). 117 ms. `off_target_scan` still stubbed — see "Step 6 carry-over" below |
| 7 | applications/deletion.py + integration | **DONE (modulo R2)** | `applications.deletion.run` orchestrates search → assemble → off-target → verify. `test_deletion_LB001_lasB.py`: 10/11 pass; `test_off_target_passes` skips locally (needs LB001.fna in cache or R2 creds). `final_plasmid.length = 6155 bp`, `AAGCTT count = 0`, all amplicon lengths byte-exact |
| 8 | applications/expression.py + integration | **NOT STARTED** | Includes `search_expression_primers` |
| 9 | applications/tagging.py + integration | **NOT STARTED** | Includes `search_tagging_primers` |
| 10 | plasmid_builder + verification full | **PARTIAL** | `assemble()` implemented for `site_destroyed` (deletion + tagging). 5 unit tests in `test_plasmid_builder.py`. Expression `site_partial_AAGCT` assembly still TODO (step 8). `RecognitionCounting` 3 unit tests still skipped (factor helper out for direct testing) |
| 11 | pdf_writer (3 templates) | **NOT STARTED** | |
| 12 | handler.py + frontend wiring | **NOT STARTED** | Astro placeholder page already routes; needs the POST endpoint + dropdowns |
| 13 | CI | **NOT STARTED** | |

## Test counts

```
pytest tests/             →  70 passed, 14 skipped
pytest tests/unit/        →  58 passed,  3 skipped
pytest tests/integration/ →  12 passed, 10 skipped (off-target / expression / tagging gated on R2 + steps 8-9)
```

Of the 14 skips:
- 3 `test_verification.py::TestRecognitionCounting` (Phase-2 helper-factoring task)
- 1 deletion `test_off_target_passes` (no LB001.fna in local cache and no R2 creds)
- 6 expression integration tests (step 8)
- 3 tagging integration tests (step 9)
- 1 lasI-deferred (Step 1)

To unskip the deletion off-target test locally, drop the genome at
`data/primer_design/genomes/LB001.fna` (or set `GENOME_LOCAL_CACHE_DIR`).

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
- ~~`off_target_scan`~~ — implemented in Step 7 (anchor walk + body extension +
  pair enumeration per skill_v2 §10).

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

- **Step 8** (expression app + integration): wire `applications/expression.py`.
  Needs `search_expression_primers` (2-primer search), and the
  `site_partial_AAGCT` branch of `plasmid_builder.assemble` for the asymmetric
  pBBR1MCS2 / HindIII junction. Test target:
  `tests/integration/test_expression_PA14_lasR.py`. PA14 isolate gene record
  must be present (currently only LB001-LB155 are bundled — PA14 needs a
  curated record built or substituted with a present isolate).
- **Step 9** (tagging app + integration): `search_tagging_primers`,
  `applications/tagging.py`, His6 cassette path. Test target:
  `tests/integration/test_tagging_PA14_lasR_His6.py`.
- **Step 10 finish**: factor `_check_recognition_count` into a public
  `count_recognition_sites(seq, enzyme)` helper and unskip the 3 unit tests in
  `test_verification.py::TestRecognitionCounting`.
- **Local off-target verification**: drop `LB001.fna` into
  `data/primer_design/genomes/` (or set `GENOME_LOCAL_CACHE_DIR`) to unskip
  `test_off_target_passes`.
