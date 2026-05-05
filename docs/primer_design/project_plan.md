# pEXG2 / pBBR1MCS2 Primer Design Tool — Project Plan

**Version:** 0.1
**Last updated:** 2026-05-04
**Companion document:** `decisions_log.md` v0.2 (Single Source of Truth for scope decisions)
**Phase:** 1 (planning + skeleton, this document) → 2 (implementation in Claude Code)

This plan describes **how** the tool is built. The **what** (scope, algorithm choices, validation rules) is fixed in the decisions log; this document references those decisions via their `D-IDs`.

Reading order for Claude Code in Phase 2: (1) decisions_log.md, (2) this file, (3) primer_design_skill_v2.md, (4) the code skeleton.

---

## 1 — System architecture

### 1.1 Component diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                  leonslab.vercel.app (single Vercel project)     │
│                                                                  │
│   ┌──────────────────────┐         ┌────────────────────────┐    │
│   │  Frontend (Next.js)  │  POST   │ /api/design/handler.py │    │
│   │  /primer-design page │ ──────▶ │ (Python serverless fn) │    │
│   │  Dropdowns + preview │ ◀────── │                        │    │
│   └──────────────────────┘  JSON   └─────────┬──────────────┘    │
│                                              │                   │
│                                              ▼                   │
│                                ┌──────────────────────────────┐  │
│                                │  /lib/primer_design/         │  │
│                                │  ─────────────────────────   │  │
│                                │  applications/{deletion,     │  │
│                                │    tagging, expression}.py   │  │
│                                │  primer_designer.py          │  │
│                                │  plasmid_builder.py          │  │
│                                │  verification.py             │  │
│                                │  pdf_writer.py               │  │
│                                │  gene_finder.py              │  │
│                                │  storage_adapter.py          │  │
│                                │  vectors.py / tags.py        │  │
│                                │  config.py                   │  │
│                                └────┬───────────────────┬─────┘  │
│                                     │                   │        │
│                  bundled (~5 MB)    ▼                   ▼ lazy   │
│                          ┌─────────────────┐  ┌─────────────────┐│
│                          │ /data/          │  │  Cloudflare R2  ││
│                          │ primer_design/  │  │  (private,      ││
│                          │ ──────────────  │  │  ~840 MB total) ││
│                          │ genes/*.fasta   │  │  ─────────────  ││
│                          │ vectors/*.gb    │  │  {isolate}.fasta││
│                          │ codon_table.json│  │  manifest.json  ││
│                          │ tail_conventions│  │                 ││
│                          └─────────────────┘  └─────────────────┘│
└──────────────────────────────────────────────────────────────────┘
```

### 1.2 Per-request data flow

```
User clicks "Design"
    │
    ▼  POST /api/design  { isolate, gene, action, tag?, enzyme, polymerase }
handler.py
    │ 1. validate request schema
    │ 2. resolve action → application module
    ▼
applications/<deletion|tagging|expression>.py
    │
    ├──▶ gene_finder.get_gene_record(isolate, gene)
    │         reads bundled curated FASTA → GeneRecord(seq, up_flank, dn_flank, meta)
    │
    ├──▶ vectors.get_vector(name) + get_tail_convention(vec, enz, app)
    │         reads bundled GenBank → VectorRecord + TailConvention
    │
    ├──▶ primer_designer.search_<application>_primers(...)
    │         exhaustive search (D5.3, D5.4) → PrimerSet
    │
    ├──▶ storage_adapter.fetch_genome(isolate)
    │         lazy R2 fetch (cached on warm instance)
    │
    ├──▶ primer_designer.off_target_scan(primers, genome)  (D5.5)
    │
    ├──▶ plasmid_builder.assemble(vector, insert, convention) → CircularPlasmid
    │
    ├──▶ verification.verify(result, request)  (D7.4 hard blocks)
    │         pass → continue ; fail → raise + return error response
    │
    ├──▶ pdf_writer.render_<application>_pdf(result, /tmp/out.pdf)
    │
    └──▶ write FASTA + JSON to /tmp + return URLs in response
```

Total runtime budget: ≤ 3 s typical, ≤ 10 s worst-case (Vercel default timeout).

### 1.3 Deployment topology

| Component | Where | Format |
|---|---|---|
| Frontend page | `app/primer-design/page.tsx` | Next.js, deployed to Vercel |
| API endpoint | `api/design/handler.py` | Vercel Python runtime, auto-detected |
| Algorithm core | `lib/primer_design/` | Python module imported by handler |
| Reference data (small) | `data/primer_design/` | bundled in function deployment |
| Reference data (large) | Cloudflare R2 bucket `leonslab-pa-genomes` | accessed via boto3, env vars in Vercel |
| Output artefacts | `/tmp/<job_id>/` (ephemeral) returned as base64 in response | inline, no persistent storage |

Output strategy decision: artefacts are returned inline (base64-encoded FASTA/PDF/JSON) in the API response — no persistent file URLs needed. Frontend triggers browser download. This avoids needing additional storage for output files and keeps the Vercel function stateless. Trade-off: response size up to ~500 KB per design (acceptable; well under Vercel limits).

---

## 2 — Data model

### 2.1 Curated gene+flank reference (bundled)

One file per `(gene × isolate)` pair. Ships in function deployment.

**Path:** `data/primer_design/genes/{gene}/{isolate_id}.fasta`

**Format:** Single-record FASTA, header carries metadata as key=value pairs.

```
>{isolate_id}|{gene}|contig={contig_id}|start={1-based-start}|end={end}|strand={+|-}|cds_len={n_nt}|variant_group={group_id}|aa_changes={comma-list-or-none}
{1500 nt up-flank then CDS then 600 nt dn-flank, on top strand if strand=+, RC if strand=-}
```

Top strand with strand-normalization means the bundled record is always 5'→3' relative to the gene's reading direction. This simplifies all downstream primer logic — no strand handling in the algorithm.

**Sequence layout in record:**
```
[600 nt up-flank][CDS, multiple of 3][600 nt dn-flank]
                  └ starts with ATG ┘
```

**Construction rule (Phase-2 task):** for each `(gene, isolate)`, extract from the hybrid assembly using the variant-annotated coordinates. If neighbour gene starts within 600 bp, shrink flank to neighbour-boundary − 30 bp (minimum 300 bp).

### 2.2 Vectors (bundled)

**Paths:**
- `data/primer_design/vectors/pEXG2.gb` — GenBank format, KM887143.1
- `data/primer_design/vectors/pBBR1MCS2.gb` — GenBank format, hand-annotated from sequence verified 2026-05-03

GenBank chosen over FASTA so that MCS / sacB / KanR / GmR / ori features are queryable for sanity checks (e.g., "is the chosen cut site inside the MCS?").

### 2.3 Tail conventions (bundled JSON)

**Path:** `data/primer_design/tail_conventions.json`

```json
{
  "pEXG2|HindIII|deletion": {
    "name": "site_destroyed",
    "p1_tail_rule": "left_arm_15nt_excl_recognition",
    "p4_tail_rule": "rc_right_arm_15nt_starting_at_second_nt_of_recognition",
    "expected_recognition_count_in_final_plasmid": 0
  },
  "pEXG2|HindIII|tagging": {
    "name": "site_destroyed",
    "p1_tail_rule": "left_arm_15nt_excl_recognition",
    "p4_tail_rule": "rc_right_arm_15nt_starting_at_second_nt_of_recognition",
    "expected_recognition_count_in_final_plasmid": 0
  },
  "pBBR1MCS2|HindIII|expression": {
    "name": "site_partial_AAGCT",
    "p1_tail_rule": "rc_right_arm_15nt_incl_5nt_recognition",
    "p4_tail_rule": "left_arm_15nt_incl_5nt_recognition",
    "expected_recognition_count_in_final_plasmid": 1
  }
}
```

Tail rules are named computational primitives implemented in `vectors.py`. Adding a new `(vector, enzyme, application)` triple = one JSON entry + (if needed) one new rule function.

### 2.4 Tag library (bundled JSON or Python module)

**Path:** `lib/primer_design/tags.py` (Python module preferred over JSON because of codon-table dependency).

```python
LINKER_GGS = "GGCGGCAGC"

@dataclass(frozen=True)
class Tag:
    name: str
    dna: str
    protein: str
    n_term_ok: bool
    c_term_ok: bool
    in_locus_ok: bool          # False if cassette > 36 nt
    cassette_nt: int           # full cassette length: linker + dna + new stop

TAGS: dict[str, Tag] = {
    "FLAG":   Tag("FLAG",   "GACTACAAGGACGACGATGACAAG", "DYKDDDDK",
                  n_term_ok=True, c_term_ok=True, in_locus_ok=True, cassette_nt=36),
    "3xFLAG": Tag("3xFLAG", "GACTACAAG...GACAAG",         "DYKDHDGD...DDDDK",
                  n_term_ok=False, c_term_ok=True, in_locus_ok=False, cassette_nt=90),
    "His6":   Tag("His6",   "CACCATCATCATCACCAC", "HHHHHH",
                  n_term_ok=True, c_term_ok=True, in_locus_ok=True, cassette_nt=30),
    "His8":   Tag("His8",   "CATCATCATCATCATCATCATCAT", "HHHHHHHH",
                  n_term_ok=True, c_term_ok=True, in_locus_ok=True, cassette_nt=36),
    "HiBiT":  Tag("HiBiT",  "GTGAGCGGCTGGCGGCTGTTCAAGAAGATCAGC", "VSGWRLFKKIS",
                  n_term_ok=False, c_term_ok=True, in_locus_ok=False, cassette_nt=45),
}
```

### 2.5 R2 genomes

**Naming:** `{isolate_id}.fasta` (e.g. `LB001.fasta`, `PA14.fasta`)
**Format:** standard multi-record FASTA (one record per contig)
**Manifest:** `data/primer_design/manifest.json`, git-tracked, contains SHA-256 per genome:

```json
{
  "schema_version": 1,
  "generated": "2026-05-04T12:00:00Z",
  "genomes": [
    {"isolate_id": "LB001", "filename": "LB001.fasta", "sha256": "...", "size_bytes": 6660000, "n_contigs": 5}
  ]
}
```

`storage_adapter.fetch_genome()` verifies SHA-256 after download.

### 2.6 Codon table (bundled JSON)

**Path:** `data/primer_design/codon_table_paeruginosa.json` — *P. aeruginosa* high-frequency codons per amino acid, version-pinned. Used only for tag reverse-translation if non-default optimization desired (current tag DNA in §2.4 is pre-optimized).

### 2.7 Request / response schemas

**Request (POST /api/design):**

```typescript
type DesignRequest = {
  isolate_id: string;            // e.g. "LB001"
  gene: "lasB" | "lasR" | "lasI";
  action: "delete" | "tag" | "express";
  tag?: "FLAG" | "3xFLAG" | "His6" | "His8" | "HiBiT" | null;
  tag_position?: "N" | "C" | null;     // ignored for action="delete"
  use_plasmid_for_tag?: boolean;       // true → switch to expression path with C-term tag
  enzyme: string;                       // e.g. "HindIII", filtered single-cutter
  polymerase: "B7" | "Phusion" | "Q5" | "Taq";
}
```

**Response (success):**

```typescript
type DesignResponse = {
  status: "ok";
  tool_version: string;
  job_id: string;                       // UUID for traceability
  artefacts: {
    fasta: string;                      // base64 of combined FASTA file
    pdf: string;                        // base64 of PDF
    json: object;                       // structured manifest, inline (not base64)
  };
  summary: {
    primers: Array<{name, sequence, tail_len, body_len, tm_body_C, gc_body_pct}>;
    final_plasmid_length_bp: number;
    scar_or_cassette_translation: string;
    annealing_temp_C: number;
  };
  warnings: string[];
}
```

**Response (error):**

```typescript
type DesignError = {
  status: "error";
  error_code: string;       // e.g. "off_target_detected", "tag_too_long_for_in_locus"
  message: string;          // human-readable
  details?: object;
}
```

---

## 3 — Module structure (public APIs)

### 3.1 `api/design/handler.py`

```python
from http.server import BaseHTTPRequestHandler

def handler(event, context):
    """Vercel Python runtime entrypoint. Parses POST body, dispatches, returns JSON."""
    # 1. Parse + validate request against DesignRequest schema (use pydantic)
    # 2. action → applications.{deletion|tagging|expression}.run(request)
    # 3. Catch typed exceptions → DesignError responses
    # 4. Return JSON with base64'd FASTA + PDF and inline JSON manifest
```

### 3.2 `lib/primer_design/applications/`

Each submodule exposes:

```python
def run(request: DesignRequest) -> DesignResult:
    """Pipeline: gene_finder → primer_designer → off_target → assemble → verify → render."""
```

`DesignResult` is a dataclass with all fields needed by `pdf_writer` and the JSON manifest.

### 3.3 `lib/primer_design/gene_finder.py`

```python
def get_gene_record(isolate_id: str, gene: str) -> GeneRecord:
    """Read curated bundled FASTA, parse header metadata, return strand-normalized record."""

@dataclass
class GeneRecord:
    isolate_id: str
    gene: str
    cds_seq: str                # always 5'→3' on coding strand, starts with ATG
    up_flank: str               # 600 nt (or shorter), 5'→3' on coding strand
    dn_flank: str               # 600 nt (or shorter), 5'→3' on coding strand
    contig_id: str
    genome_start_1based: int
    genome_end_1based: int
    original_strand: Literal["+", "-"]
    variant_group: str | None
    aa_changes: list[str]
```

### 3.4 `lib/primer_design/primer_designer.py`

```python
def search_deletion_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
    config: PrimerConfig,
) -> DeletionPrimerSet:
    """Exhaustive (N, C) × P1 × P2 × P3 × P4 search with Tm-NN scoring (D5.1, D5.4).
       Returns top-1 with top-5 alternatives for inspection."""

def search_tagging_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
    tag_cassette: str,           # GGS + tag DNA + new stop
    config: PrimerConfig,
) -> TaggingPrimerSet:
    """Like search_deletion_primers but with tag cassette injected into P2/P3 junction.
       Junction overlap extended to len(tag_cassette) per D6.4."""

def search_expression_primers(
    gene: GeneRecord,
    vector: VectorRecord,
    convention: TailConvention,
    rbs_tail: str,               # AGGAGG + synthetic_spacer (D6.1.v2)
    tag: Tag | None,
    tag_position: Literal["N", "C"] | None,
    config: PrimerConfig,
) -> ExpressionPrimerSet:
    """2-primer search. P1 carries vector_tail + rbs_tail + body at ATG.
       P2 carries vector_tail + body at stop (RC). With tag, modifies as per §6.2 of decisions_log."""

def off_target_scan(
    primer_set: PrimerSet,
    genome_fasta: bytes,
    config: PrimerConfig,
) -> OffTargetReport:
    """Walk-anchor scan + full-body match per D5.5. Returns all expected and unexpected products."""
```

### 3.5 `lib/primer_design/plasmid_builder.py`

```python
def assemble(
    vector: VectorRecord,
    insert: str,                  # for deletion/tagging: UP+DN combined; for expression: PCR amplicon
    cut_position_in_vector: int,
    convention: TailConvention,
) -> CircularPlasmid:
    """Build the final circular plasmid by replacing the cut region with the insert.
       Verifies junction overlaps. Returns sequence (linear representation) + length + junction coords."""
```

### 3.6 `lib/primer_design/verification.py`

```python
def verify(result: DesignResult, request: DesignRequest) -> VerificationReport:
    """All hard checks per D7.4. Raises VerificationFailure on any failure."""

class VerificationFailure(Exception):
    failed_check: str
    details: dict
```

### 3.7 `lib/primer_design/pdf_writer.py`

```python
def render_deletion_pdf(result: DeletionResult, output_path: Path) -> None: ...
def render_tagging_pdf(result: TaggingResult, output_path: Path) -> None: ...
def render_expression_pdf(result: ExpressionResult, output_path: Path) -> None: ...
```

Common `_render_header()`, `_render_primer_table()`, `_render_off_target_table()`, `_render_pcr_conditions()` helpers shared between the three templates.

### 3.8 `lib/primer_design/storage_adapter.py`

```python
@functools.lru_cache(maxsize=4)
def fetch_genome(isolate_id: str) -> bytes:
    """Lazy R2 fetch with SHA-256 verification against manifest. Cached on warm instance."""

def _r2_client():
    """boto3 client with config from env vars R2_ACCOUNT_ID, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY."""
```

LRU cache size 4 trades memory (~28 MB peak in /tmp) for speed on repeated requests within a warm instance.

### 3.9 `lib/primer_design/vectors.py`

```python
def get_vector(name: Literal["pEXG2", "pBBR1MCS2"]) -> VectorRecord:
    """Read bundled GenBank, return parsed structure."""

def get_tail_convention(vector: str, enzyme: str, application: str) -> TailConvention:
    """Lookup from tail_conventions.json. Raises ConventionNotCalibrated on miss."""

def get_unique_cutters(vector: VectorRecord) -> list[str]:
    """Enumerate all unique restriction sites in the vector. Used by frontend to populate enzyme dropdown."""

def find_cut_position(vector: VectorRecord, enzyme: str) -> int:
    """Returns the 0-based nick position on top strand. Raises if not exactly one site."""
```

### 3.10 `lib/primer_design/config.py`

Constants. No imports beyond stdlib.

```python
TM_NN_PARAMS = {"Na": 50, "dnac1": 500, "dnac2": 0}      # D5.1
GC_RANGE = (0.40, 0.70)                                   # D5.3, P. aeruginosa
BODY_LEN_RANGE = (18, 28)
TM_PREFERRED_RANGE = (60, 64)
TM_HARD_RANGE = (58, 66)
TM_SPREAD_HARD_LIMIT = 3.0                                # °C
SCAR_MAX = 24                                             # D5.4
UP_FLANK = DN_FLANK = 600
JUNCTION_LEN_DEFAULT = 30
JUNCTION_LEN_EXTENDED = 36
SYNTHETIC_SPACER_DEFAULT = "ACTTGTTC"                     # D6.1.v2
CANONICAL_RBS = "AGGAGG"
LINKER_GGS = "GGCGGCAGC"
POLYMERASE_OFFSETS_C = {"B7": 0.0, "Phusion": 3.0, "Q5": 3.0, "Taq": -5.0}     # D5.2
```

---

## 4 — Workflow flow diagrams

### 4.1 Deletion (sequence)

```
Frontend → handler.py → applications/deletion.py
                                 │
                                 ├─▶ gene_finder.get_gene_record(isolate, gene)
                                 │     returns GeneRecord
                                 │
                                 ├─▶ vectors.get_vector("pEXG2")
                                 │   vectors.get_tail_convention("pEXG2", enzyme, "deletion")
                                 │
                                 ├─▶ primer_designer.search_deletion_primers(...)
                                 │     │ for N in 1..23:
                                 │     │   for C in 1..(24-N):
                                 │     │     for P1 in pool:
                                 │     │       for P2 in pool:
                                 │     │         ... score → keep top
                                 │     returns DeletionPrimerSet
                                 │
                                 ├─▶ storage_adapter.fetch_genome(isolate)
                                 │     ~7 MB R2 fetch + SHA-256 verify
                                 │
                                 ├─▶ primer_designer.off_target_scan(primers, genome)
                                 │     verifies P1+P2, P3+P4, P3+P2, P1+P4 outcomes (D5.5)
                                 │
                                 ├─▶ plasmid_builder.assemble(vector, UP+DN, cut, convention)
                                 │
                                 ├─▶ verification.verify(result, request)
                                 │     hard checks (D7.4) → raise on fail
                                 │
                                 └─▶ pdf_writer.render_deletion_pdf(result, out_path)
                                       returns artefacts → handler returns JSON
```

### 4.2 Tagging (delta from deletion)

Same as deletion **except**:
- `search_tagging_primers` builds `tag_cassette = LINKER_GGS + tag.dna + "TAA"` (~30 or 36 nt) and injects it into the P2/P3 junction overlap.
- Junction overlap length = `len(tag_cassette)` instead of fixed 30 nt.
- DN fragment starts AT native stop position (0-based offset 0 of `dn_flank` IS the native stop, gets replaced).
- Verification check: assembled fusion ORF reads `gene + GGS + tag` then stops, no internal stops.

### 4.3 Expression (sequence)

```
Frontend → handler.py → applications/expression.py
                                 │
                                 ├─▶ gene_finder.get_gene_record(isolate, gene)
                                 │
                                 ├─▶ vectors.get_vector("pBBR1MCS2")
                                 │   vectors.get_tail_convention("pBBR1MCS2", enzyme, "expression")
                                 │
                                 ├─▶ primer_designer.search_expression_primers(...)
                                 │     P1 = vector_tail + AGGAGG + ACTTGTTC + body_at_ATG
                                 │     P2 = vector_tail + body_at_stop (with tag if tag_position="C")
                                 │     If tag_position="N": P1 carries tag DNA between RBS and body
                                 │
                                 ├─▶ storage_adapter.fetch_genome(isolate)
                                 ├─▶ primer_designer.off_target_scan(...)
                                 ├─▶ plasmid_builder.assemble(...)
                                 ├─▶ verification.verify(...)
                                 └─▶ pdf_writer.render_expression_pdf(...)
```

---

## 5 — Test strategy

### 5.1 Unit tests (`tests/unit/`)

| Module | Test cases |
|---|---|
| `gene_finder.py` | round-trip: parse FASTA → GeneRecord → re-export → byte-equal. Strand normalization for - strand gene. Multiple-of-3 CDS check, single-stop check. |
| `primer_designer.py` (filters) | Known-good primer body passes; each hard filter (clamp, homopolymer, GC, length, dimer, Tm) tested with known violators. |
| `primer_designer.py` (Tm) | Tm_NN against published Allawi & SantaLucia values for ≥10 reference oligos, tolerance ±0.5 °C. |
| `vectors.py` (tail derivation) | pEXG2+HindIII produces `CATAAATGTAAAGCA` for P1 tail (matches v1 working primer). pBBR1MCS2+HindIII produces `CGGTATCGATAAGCT` for P1 tail (matches PA14 lasB primer #1573). |
| `tags.py` | Each tag's cassette length matches `cassette_nt`. In-locus capability flag matches `cassette_nt ≤ 36`. |
| `verification.py` | Each hard-block raises with correct error code; passing case completes silently. |

### 5.2 Integration tests (`tests/integration/`)

End-to-end through `applications/<X>.run()` with bundled mini-genome (single contig containing the gene + flanks, no R2 fetch). 

| Test | Application | Expected |
|---|---|---|
| `test_deletion_LB001_lasB_pEXG2_HindIII` | deletion | Final plasmid SHA-256 matches v1 reference (`pEXG2-lasB_delta_LB001.fasta`). 4 primers reproduce empirical sequences within ±2 nt body shift (algorithmic equivalence). |
| `test_expression_PA14_lasR_pBBR1MCS2_HindIII` | expression | P1 tail equals `cggtatcgataagctaggaggacttgttc`, P1 body starts with `ATGGCCTTGGTT...` (matches empirical #1612). |
| `test_tagging_PA14_lasR_pEXG2_HindIII_His6` | tagging | Junction overlap contains `GGCGGCAGC` + `CACCATCATCATCACCAC` + `TAA`. Final plasmid translates as LasR-GGS-HHHHHH-stop. |
| `test_tag_too_long_HiBiT_in_locus` | tagging | Returns error `tag_too_long_for_in_locus` with suggestion to switch to plasmid expression. |

### 5.3 Golden regression suite (`tests/golden/`)

Fully reproducible end-to-end runs with all artefacts (FASTA + PDF byte-equal up to timestamp; JSON byte-equal up to timestamp). One golden run per application + edge case. Updated explicitly via `pytest --update-golden` when intentional change to algorithm.

### 5.4 Empirical lab validation (post-deploy)

Phase 2 acceptance gate, **not automatable**: Leon (or designate) runs ≥3 designs through the deployed tool, orders the primers, runs PCR + assembly + sequencing, confirms ≥80 % first-pass success rate. If <80 %, parameters need tuning (likely: GC range, Tm range, off-target sensitivity).

---

## 6 — Phase-2 roadmap

Ordered list of implementation tasks. Each item ends with an acceptance gate.

| # | Task | Est. | Acceptance |
|---|---|---|---|
| 1 | Build curated gene+flank reference for 3 genes × ~120 isolates from existing assemblies | 1d | 360 FASTAs in `data/primer_design/genes/`, all CDS multiple-of-3, all start with ATG, all single internal stop |
| 2 | Annotate vectors as GenBank, load with Biopython, expose via `vectors.py` | 0.5d | `get_vector` returns parsed record; `get_unique_cutters` enumerates ≥10 sites for each vector |
| 3 | Calibrate tail conventions against empirical primers (pEXG2 v1, pBBR1MCS2 PA14 examples) | 0.5d | Unit tests in §5.1 row "vectors.py" pass |
| 4 | Run `scripts/upload_genomes.py` to push 120 genomes to R2; commit manifest.json | 1h | All 120 genomes in R2, manifest hashes verify |
| 5 | Implement `gene_finder.py` + `storage_adapter.py` | 1d | Unit tests pass |
| 6 | Implement `primer_designer.py` core (Tm, filters, scoring) | 2d | Unit tests pass |
| 7 | Implement `applications/deletion.py` + integration test | 1.5d | `test_deletion_LB001_lasB_pEXG2_HindIII` passes |
| 8 | Implement `applications/expression.py` + integration test | 1d | `test_expression_PA14_lasR_pBBR1MCS2_HindIII` passes |
| 9 | Implement `applications/tagging.py` + integration test | 2d | `test_tagging_PA14_lasR_pEXG2_HindIII_His6` passes |
| 10 | Implement `plasmid_builder.py` + `verification.py` | 1d | All hard checks raise correctly; passing cases produce SHA-equal plasmid to v1 reference |
| 11 | Implement `pdf_writer.py` (3 templates) | 2.5d | Visual inspection: PDF matches v1 PDF style; all sections render |
| 12 | Implement `handler.py` + frontend integration on /primer-design | 1.5d | End-to-end: dropdowns → POST → response with downloadable artefacts |
| 13 | Set up CI: run unit + integration + golden on every PR | 0.5d | GitHub Actions green |
| 14 | Empirical validation panel (ext.) | open | ≥3 designs ordered + tested |

**Total estimate:** ~16 days for Claude Code, parallelizable to ~10 calendar days with some idle on lab-validation feedback.

---

## 7 — User setup tasks

### Before Phase 2

- [x] Cloudflare R2 account, bucket `leonslab-pa-genomes` created
- [x] Vercel env vars `R2_*` set in Production + Preview + Development
- [x] Scaffolding prompt run in leonslab repo, folder structure in place
- [ ] Provide hybrid assemblies (`*.fasta` per isolate) for the 120 isolates as a single archive — Claude Code parses these in roadmap step 1 + 4

### During Phase 2

- [ ] Review each module PR before merging; especially the calibration step (#3) and the first plasmid-builder integration (#7)
- [ ] Run integration + golden tests locally before deploying
- [ ] Deploy to Vercel preview, click through `/primer-design` end-to-end at least once before Production

### After Phase 2 launch

- [ ] Order primers from a fresh design and verify in the lab (acceptance gate)
- [ ] When new isolates are sequenced: drop FASTA into `data/primer_design/genes/` (extracted), upload genome to R2, append manifest entry, commit
- [ ] When new genes added: same as above plus update `lib/primer_design/config.py` gene list

---

## 8 — Risks and mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Vercel cold start > 3 s, user perceives slowness | Medium | Low | Show frontend spinner with "warming up..." for first request; reuse warm instance after |
| R2 fetch latency on first per-instance request | High | Low | Acceptable (~200 ms); after warm cache, free |
| Tag cassette > 36 nt requested for in-locus | Certain | Low | Reject with clear suggestion to use plasmid (D6.4) |
| sacB-tolerant clinical isolate | Medium | Medium | PDF includes pre-test instructions; not the tool's responsibility to detect |
| Off-target primer in unusual genomic context | Medium | High | Per-isolate genome scan as hard block (D5.5, D7.4) |
| Algorithm finds zero valid primer sets | Medium | High | Return error with `relax_filters_suggestion` listing which filter is binding; documented in PDF as "if this happens..." |
| Vendor-lock to Vercel / R2 | Low | Low | Storage adapter abstracted; can switch to Hugging Face Datasets without algorithm changes |
| pEXG2 / pBBR1MCS2 sequence differs from canonical (cloning errors in lab stock) | Low | High | Vector loaded from bundled GenBank, hash-verified against manifest |

---

## 9 — Conventions

- Python 3.11+ (matches Vercel current Python runtime)
- Type hints required on all public APIs (`def f(x: int) -> str:`)
- Docstrings: Google-style
- Tests: `pytest`, no parallel runners
- Lint: `ruff` with default ruleset
- Format: `black`, line length 100
- Logging: stdlib `logging`, JSON formatter for production (Vercel parses)
- Errors: typed exceptions per module; surface to handler as JSON error responses
- Imports: stdlib → third-party → local, separated by blank lines (ruff enforces)
- No global mutable state; all config goes through `config.py` constants

---

## 10 — What this document is NOT

- Not the algorithm specification — that's `primer_design_skill_v2.md` (next).
- Not the API reference — that's auto-generated from docstrings during Phase 2.
- Not the user manual — that's a `/primer-design/help` page, written after launch.
- Not exhaustive on every edge case — those go in skill_v2 and individual module docstrings.

---

## 11 — Change log

| Version | Date | Change |
|---|---|---|
| 0.1 | 2026-05-04 | Initial document. Architecture, data model, module APIs, workflows, test strategy, Phase-2 roadmap, user setup tasks, risks. |
