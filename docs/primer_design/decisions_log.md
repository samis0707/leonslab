# pEXG2 / pBBR1MCS2 Primer Design Tool — Decisions Log

**Document version:** 0.2
**Last updated:** 2026-05-04
**Project name (working title):** primer-design (sub-page of leonslab.vercel.app)
**Phase:** 1 (planning + skeleton). Phase 2 = implementation in Claude Code.

This document captures every architecture, algorithm, and scope decision made during the planning phase. Single source of truth for Phase 2. Append-only history; superseded decisions marked `SUPERSEDED` and kept for traceability. Cross-reference IDs in code comments.

Status legend: **DECIDED** • **PROPOSED** • **DEFERRED** • **OPEN** • **SUPERSEDED**

---

## 1 — Project scope

**D1.1.v2** Tool designs primers and assembles plasmid sequences for **three cloning applications** in v1:
- **Allelic exchange — in-frame deletion** via pEXG2 with In-Fusion 3-fragment assembly (vector + UP + DN).
- **Allelic exchange — in-locus C-terminal tag insertion** via pEXG2 with In-Fusion 3-fragment assembly (vector + UP + DN; tag cassette engineered into the P2/P3 junction overlap).
- **Constitutive expression** via pBBR1MCS2 with In-Fusion 2-fragment assembly (vector + insert), optionally with N- or C-terminal tag on the plasmid-encoded copy.

Status: DECIDED. Future applications kept extensible: complementation with native promoter, induced expression, point-mutation introduction, transcriptional reporter.

*D1.1.v1 SUPERSEDED 2026-05-04 — original two-application list (deletion + expression). Tagging promoted from "tag option within expression" to its own application type because the user prefers in-locus tagging at the native locus as default for tag experiments.*

**D1.2** Target genes for v1: `lasB`, `lasR`, `lasI`. Adding a new gene = drop a `(gene × isolate)` FASTA into `data/genes/`, no code changes.
Status: DECIDED.

**D1.3** Isolate panel: ~120 clinical *P. aeruginosa* isolates plus PA14 and PAO1 references.
Status: DECIDED.

---

## 2 — Workflow / user interaction

**D2.1** Operating mode: **lookup-only**. Curated database; no genome upload.
Status: DECIDED.

**D2.2.v2** User-facing inputs (UI dropdowns/radios):

1. **Isolate** (dropdown, ~120 entries)
2. **Gene** (`lasB` | `lasR` | `lasI`)
3. **Action** (radio):
   - "Delete this gene" → vector = pEXG2, application = `deletion`
   - "Tag this gene (C-terminal)" → vector = pEXG2, application = `tagging` *(default)*. Optional checkbox "Use plasmid expression instead" → switches to vector = pBBR1MCS2, application = `expression`, tag position = C
   - "Express this gene from plasmid" → vector = pBBR1MCS2, application = `expression`
4. **Tag** (when applicable): `none` | `FLAG` | `3xFLAG` | `His6` | `His8` | `HiBiT`. Default `none` for plasmid expression; default `His6` for in-locus tagging.
5. **Tag position** (plasmid expression only, gated by tag): `N` | `C`. In-locus tagging: always `C`.
6. **Restriction enzyme** (filtered to single-cutters of the chosen vector; default `HindIII`)
7. **Polymerase** (default `B7`)

No more RBS-mode dropdown — RBS is fixed (D6.1.v2).

Status: DECIDED.

*D2.2.v1 SUPERSEDED 2026-05-04 — RBS-mode dropdown removed; action-based UI replacing flat dropdown layout.*

**D2.3** Output artefacts per run:
- Combined annotated FASTA (primers with tail/body markup, amplicons, insert, full circular plasmid)
- PDF report (English, ReportLab)
- JSON parameter manifest
- Plasmid-map preview rendered live in UI (e.g. SeqViz)

Status: DECIDED.

---

## 3 — Hosting and architecture

**D3.1** Frontend: integrated as sub-page in existing `leonslab.vercel.app` navigation menu. Same Next.js project.
Status: DECIDED.

**D3.2** Backend: Python serverless function on Vercel at `/api/design`. Default 10 s timeout sufficient.
Status: DECIDED.

**D3.3** Reference data — gene + 600 bp UP + 600 bp DN per `(gene × isolate)`, ~5 MB total — bundled into function deployment.
Status: DECIDED.

**D3.4** Reference data — full genomes ~840 MB — stored in **Cloudflare R2** (free tier 10 GB, zero egress, S3-compatible). Private bucket, accessed via Vercel env vars.
- ENV vars: `R2_ACCOUNT_ID`, `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY`, `R2_BUCKET_NAME`. **All four set in Vercel as of 2026-05-04 (Production + Preview + Development, non-Sensitive).**
- Lazy fetch ~7 MB per request into `/tmp` (Vercel allows 500 MB), cached on warm instances.
- Fallback if R2 unavailable: Hugging Face Datasets. Storage adapter abstracted.

Status: DECIDED.

**D3.5** One-time setup: `scripts/upload_genomes.py` performs initial R2 upload. SHA-256 manifest committed to git.
Status: DECIDED.

---

## 4 — Vector database

**D4.1** Two vectors hardcoded in v1: pEXG2 (KM887143.1, 5084 bp), pBBR1MCS2 (5148 bp, sequence verified 2026-05-03).
Status: DECIDED.

**D4.2** Restriction-enzyme dropdown computed dynamically.

pBBR1MCS2 single-cutters (13): HindIII, EcoRI, BamHI, SalI, XhoI, KpnI, SacI, SpeI, XbaI, ApaI, ClaI, EcoRV, SmaI. Excluded: NotI, PstI (2 sites).

pEXG2: enumerate during DB construction in Phase 2.

Status: DECIDED.

**D4.3** **Tail convention** per `(vector × cut_enzyme × application)` triple, calibrated config:

| Triple | Convention | P1 tail | P4 tail | Recognition count in final plasmid |
|---|---|---|---|---|
| pEXG2 + HindIII + deletion | `site_destroyed` | last 15 nt of left arm before AAGCTT | RC of right arm 15 nt starting at second nt of AAGCTT | 0 |
| pEXG2 + HindIII + tagging | `site_destroyed` | same as deletion | same as deletion | 0 |
| pBBR1MCS2 + HindIII + expression | `site_partial_AAGCT` | RC of right arm 15 nt incl. 5 nt recognition (`AGCTT...`) | left arm 15 nt incl. 5 nt recognition (`...AAGCT`) | 1 (regenerated at one junction) |
| Other triples | `site_destroyed` (default) | by formula | by formula | 0 |

Insert orientation:
- Deletion / tagging: arbitrary, tool picks convention minimizing off-target risk.
- Expression: enforced so that vector promoter (Plac on bottom strand of pBBR1MCS2) reads through gene 5'→3'.

Status: DECIDED.

**D4.4** Validation per convention is hard (refuse output on deviation). See §7.4.
Status: DECIDED.

---

## 5 — Algorithm

**D5.1** Tm formula: **Nearest-Neighbor (Allawi & SantaLucia 1997, Na correction SantaLucia 1998)** via `Bio.SeqUtils.MeltingTemp.Tm_NN(Seq(body), Na=50, dnac1=500, dnac2=0)`.
Status: DECIDED.

**D5.2** Annealing temperature in PDF: Tm of weakest primer body, polymerase offsets:
- B7 (default) ±0 °C • Phusion +3 °C • Q5 +3 °C • Taq −5 °C

Status: DECIDED.

**D5.3** Primer body hard filters: 3' G/C clamp; no 4-homopolymer; GC body 40–70 %; length 18–28 nt; 3' self-dimer ≤ 4; Tm body 58–66 °C, prefer 60–64 °C. Selection score and Tm-spread cutoff (≤ 3 °C) per v1 skill, sections 5–7.
Status: DECIDED.

**D5.4** Scar design (deletion only): exhaustive `(N, C)` search, `N + C ≤ 24`, asymmetric splits allowed.
Status: DECIDED.

**D5.5** Off-target scan: target isolate genome only.
Status: DECIDED.

**D5.6** Primer3 secondary cross-validation: deferred to v1.1.
Status: DEFERRED.

---

## 6 — Application-specific configuration

**D6.1.v2** RBS strategy for **plasmid expression** (pBBR1MCS2): **fixed, no user choice in UI**.
- Tail = vector(15) + canonical Shine-Dalgarno `AGGAGG`(6) + synthetic 8-nt spacer (8). Body starts at ATG.
- **Default synthetic spacer = `acttgttc`** (matches working PA14 lasR/lasI primers #1611/#1612). Editable via config file, not via UI.
- Rationale: high translation prioritized over native regulatory context; canonical Shine-Dalgarno performs empirically more robustly across the isolate panel.

Status: DECIDED.

*D6.1.v1 SUPERSEDED 2026-05-04 — three-mode user choice (canonical+native_spacer / canonical+synthetic_spacer / native_RBS) collapsed to single fixed strategy.*

**D6.2.v2** Tag library:

| Tag | DNA (P. aeruginosa codon-optimized) | Protein | Allowed positions | In-locus capable? |
|---|---|---|---|---|
| FLAG | `GACTACAAGGACGACGATGACAAG` | DYKDDDDK | N or C | Yes (extended P2/P3 tails, body 24+9+3=36 nt) |
| 3×FLAG | `GACTACAAGGACCACGACGGCGACTACAAGGATCATGATATCGATTACAAGGATGACGATGACAAG` | DYKDHDGDYKDHDIDYKDDDDK | C only | **No** (78 nt cassette exceeds 4-primer design; plasmid expression only in v1) |
| His₆ | `CACCATCATCATCACCAC` | HHHHHH | N or C | Yes (18+9+3=30 nt, fits exactly) |
| His₈ | `CATCATCATCATCATCATCATCAT` | HHHHHHHH | N or C | Yes (24+9+3=36 nt, extended tails) |
| HiBiT | `GTGAGCGGCTGGCGGCTGTTCAAGAAGATCAGC` | VSGWRLFKKIS | C only | **No** (33+9+3=45 nt cassette; plasmid expression only in v1) |

Linker: `GGCGGCAGC` (= GGS).

**Position semantics:**
- **In-locus C-terminal tagging (pEXG2):** native stop codon replaced by `GGS-tag-newstop`. UP fragment carries gene up to last codon before stop; DN fragment starts at native stop position; tag cassette engineered into P2 + P3 junction overlap (extended to 36 nt for FLAG/His₈, 30 nt for His₆).
- **Plasmid N-terminal:** `M-tag-GGS-CDS_from_codon_2`. Native ATG dropped; tag's first ATG = new start.
- **Plasmid C-terminal:** `CDS_to_last_codon-GGS-tag-newstop`. Native stop dropped.

Reading-frame integrity verified post-assembly (no internal stops). Codon optimization uses *P. aeruginosa* high-frequency codon table, version-pinned in tool.

Status: DECIDED.

**D6.3** Expression-mode primer count: **2 primers** (P1 fwd at RBS+ATG, P2 rev at stop). 2-fragment In-Fusion.
Status: DECIDED.

**D6.4** **In-locus C-terminal tagging mechanics** (NEW in v0.2):

- 4-primer 3-fragment In-Fusion: linearized vector + UP + DN.
- UP fragment: 600 bp upstream + full gene CDS up to and including last codon BEFORE native stop.
- DN fragment: starts AT native stop position (the stop codon itself becomes the first 3 nt of DN if `site_destroyed` convention is used) + 600 bp downstream.
- Junction (P2 + P3 overlap): encodes `GGS linker + tag DNA + new stop codon` between UP and DN.
- Junction overlap length: 30 nt (His₆) or 36 nt (FLAG, His₈) — extends standard 30-nt In-Fusion overlap proportionally for longer tags.
- Tags exceeding 36 nt (3×FLAG, HiBiT): not supported in v1 in-locus mode. Tool surfaces a UI warning suggesting plasmid expression instead.
- 5-primer 4-fragment design with separate synthesized tag cassette (gBlock): deferred to v1.1.

Status: DECIDED.

---

## 7 — Output and validation

**D7.1** FASTA output records per run (single combined file): primers (with tail/body markup) → amplicons → insert (deletion/tagging only) → full circular plasmid → scar/cassette ORF + translation.
Status: DECIDED.

**D7.2** PDF output: English. **Three templates** (deletion / tagging / expression), structurally identical, application-specific sections vary. ReportLab + DejaVu fonts.
Status: DECIDED.

**D7.3** JSON output: complete reproducibility manifest (schema in skeleton).
Status: DECIDED.

**D7.4** Hard validation blocks (refuse output on any failure):
- Cut-enzyme recognition count in final plasmid matches convention.
- Off-target scan returns no unintended PCR products.
- Body reconstruction matches genome.
- For deletion: scar ORF starts M, ends single `*`, frame-correct.
- For tagging: assembled fusion starts M, ends single `*`, no internal `*`, GGS+tag at expected position before stop.
- For expression: cassette in-frame from start to stop, no internal stops, RBS at expected position, tag at expected position if requested.
- Junction overlap exact: last 30 (or 36) nt of UP amplicon == first 30 (or 36) nt of DN amplicon.

Status: DECIDED.

**D7.5** Soft validation warnings (mention in PDF, do not block): cut-site present in insert; Tm spread > 2 °C; GC body outside 45–55 %; native stop within tag-fusion region (expression).
Status: DECIDED.

---

## 8 — Phase-1 deliverables

**D8.1** Files exported from this chat session:

| File | Status |
|---|---|
| `decisions_log.md` (v0.2) | EXPORTED 2026-05-04 |
| `project_plan.md` | next turn |
| `primer_design_skill_v2.md` | turn after |
| `code_skeleton/` | final turn |

---

## 9 — Open items

| ID | Item | Note |
|---|---|---|
| O9.1 | pEXG2 single-cutter enumeration (D4.2) | Phase 2 vector-DB construction |
| O9.3 | Test-isolate regression suite | Phase 2; LB001-ΔlasB-pEXG2 = first canonical test |
| O9.10 | 5-primer 4-fragment in-locus tagging for 3×FLAG / HiBiT | v1.1 |

Closed:
- ~~O9.2~~ R2 setup credentials → DONE 2026-05-04 (env vars set in Vercel)
- ~~O9.5~~ Claude-Code scaffolding prompt → ISSUED in chat, awaiting execution

---

## 10 — Change log

| Version | Date | Change |
|---|---|---|
| 0.1 | 2026-05-04 | Initial document, all D1.x – D7.x DECIDED based on chat 2026-05-03/04 |
| 0.2 | 2026-05-04 | RBS simplified (D6.1.v1 → D6.1.v2): single fixed strategy, no UI choice. Application list grew from 2 to 3 (D1.1.v1 → D1.1.v2): in-locus C-terminal tagging via pEXG2 promoted to first-class application. UI redesigned around 3-action radio (D2.2.v1 → D2.2.v2). New D6.4 specifies in-locus tagging mechanics. Tag table updated with in-locus capability column. R2 env vars marked as set in Vercel. |
