# Primer Design Skill v2 — pEXG2 / pBBR1MCS2 In-Fusion

**Document version:** 2.0
**Last updated:** 2026-05-04
**Supersedes:** primer_design_skill_v1 (LB001-ΔlasB-pEXG2 demo skill)
**Companion:** `decisions_log.md` v0.2, `project_plan.md` v0.1

This is the **algorithmic specification** for the primer-design tool. It describes every computation the tool performs, with pseudocode at a level that maps 1:1 to the planned module structure (`project_plan.md` §3). `D-IDs` reference the decisions log.

Phase-2 implementation reads this document, follows the pseudocode, validates against the worked examples in §15.

---

## 0 — Overview

The tool produces ready-to-order primers, an assembled plasmid sequence, and a PDF report for one of three In-Fusion cloning applications:

| Application | Vector | Mode | Primers | Fragments | Junction overlap |
|---|---|---|---|---|---|
| Deletion | pEXG2 | in-frame in-locus deletion | 4 (P1, P2, P3, P4) | 3 (vector + UP + DN) | 30 nt |
| Tagging | pEXG2 | C-term tag at native locus | 4 (P1, P2, P3, P4) | 3 (vector + UP + DN) | 30 or 36 nt |
| Expression | pBBR1MCS2 | constitutive plasmid expression | 2 (P1, P2) | 2 (vector + insert) | 15 nt each side |

All three share: Tm-NN (Allawi & SantaLucia 1997, D5.1) for primer body Tm, identical hard filters on body composition (D5.3), per-`(vector, enzyme, application)` calibrated tail conventions (D4.3), the same off-target scan logic (D5.5), and the same hard validation gates (D7.4).

---

## 1 — Data ingestion

The tool ships with **canonical per-(gene, isolate) reference records** in `data/primer_design/genes/{gene}/{isolate}.fasta`. These are produced ONCE by a build pipeline that processes user-provided files (one-off or whenever new isolates arrive).

### 1.1 User-provided input

| Input | File pattern | Format | Source |
|---|---|---|---|
| Full genome assemblies | `LBxxx.fna`, `PA14.fna`, `PAO1.fna` | multi-record FASTA, one record per contig | hybrid assemblies |
| Per-gene CDS lists | `lasB.fasta`, `lasR.fasta`, `lasI.fasta` | multi-record FASTA, one record per isolate, **CDS only**, header carries isolate ID | curated |

Header format expected in gene FASTAs (flexible, parser handles common patterns):
```
>LB001|lasB
>LB014 lasB
>PA14_lasB
```
The parser extracts isolate ID by splitting on `|`, space, or first underscore-after-prefix; a regex fallback `(LB\d+|PA14|PAO1)` is the last line of defense.

### 1.2 Build pipeline (`scripts/build_gene_records.py`)

Pseudocode:

```python
def build_gene_records(
    gene_fasta_dir: Path,        # contains lasB.fasta, lasR.fasta, lasI.fasta
    genome_dir: Path,            # contains LB001.fna, ..., PA14.fna, PAO1.fna
    output_dir: Path,            # data/primer_design/genes/
    up_flank: int = 600,
    dn_flank: int = 600,
    min_flank: int = 300,
):
    for gene_file in gene_fasta_dir.glob("*.fasta"):
        gene = gene_file.stem                     # "lasB"
        for record in SeqIO.parse(gene_file, "fasta"):
            isolate_id = parse_isolate_id(record.id)   # "LB001"
            cds_seq = str(record.seq).upper()
            assert_cds_valid(cds_seq, isolate_id, gene)

            genome_path = genome_dir / f"{isolate_id}.fna"
            if not genome_path.exists():
                log.warning(f"genome missing for {isolate_id}, skipping {gene}")
                continue

            hit = locate_cds_in_genome(cds_seq, genome_path)
            if hit is None:
                log.warning(f"{gene} CDS not found in {isolate_id} genome — possible deletion or assembly gap")
                continue
            if hit.is_ambiguous:
                raise BuildError(f"{gene} matches multiple loci in {isolate_id} (paralog?); manual review needed")

            record_5to3 = strand_normalize(hit, genome_path, up_flank, dn_flank, min_flank)
            write_canonical_fasta(output_dir / gene / f"{isolate_id}.fasta", record_5to3, gene, isolate_id, hit)

    update_manifest(output_dir.parent / "manifest.json", genome_dir, gene_fasta_dir)
```

Key sub-functions:

```python
def assert_cds_valid(cds: str, isolate: str, gene: str):
    if len(cds) % 3 != 0:
        raise BuildError(f"{gene} CDS in {isolate} not multiple of 3")
    if not cds.startswith("ATG"):
        raise BuildError(f"{gene} CDS in {isolate} does not start with ATG")
    protein = Seq(cds).translate()
    if protein.count("*") != 1 or not protein.endswith("*"):
        raise BuildError(f"{gene} CDS in {isolate} has internal or no stop codon")

def locate_cds_in_genome(cds: str, genome_path: Path) -> GenomeHit | None:
    """Search both strands of all contigs for an exact CDS match.
       Exact match required — these are the curated CDS the user provided. No fuzzy match."""
    cds_rc = rc(cds)
    for contig in SeqIO.parse(genome_path, "fasta"):
        seq = str(contig.seq).upper()
        for start in find_all(seq, cds):
            yield GenomeHit(contig.id, start, start + len(cds), strand="+")
        for start in find_all(seq, cds_rc):
            yield GenomeHit(contig.id, start, start + len(cds), strand="-")
    # collapse: if multiple hits with same coordinates, fine; if multiple distinct loci, raise

def strand_normalize(hit, genome_path, up, dn, min_flank) -> str:
    """Extract gene + flanks, RC if on minus strand so result is always 5'→3' on coding strand.
       Shrink flank to min_flank if neighbour boundary closer."""
    contig = read_contig(genome_path, hit.contig_id)
    if hit.strand == "+":
        up_start = max(0, hit.start - up)
        dn_end = min(len(contig), hit.end + dn)
        sequence = contig[up_start:dn_end]
    else:
        up_start = max(0, hit.start - dn)              # in original genome, dn flank is "before" the gene
        dn_end = min(len(contig), hit.end + up)
        sequence = rc(contig[up_start:dn_end])
    enforce_min_flank(sequence, hit, min_flank)
    return sequence
```

### 1.3 Canonical record format on disk

`data/primer_design/genes/lasB/LB001.fasta`:
```
>LB001|lasB|contig=LB001_00001|start=1848921|end=1850417|strand=+|cds_len=1497|up=600|dn=600
ACGTACGTACGT...   (60 chars/line, total ~2700 nt for default flanks)
```

Sequence layout: `[up_flank ≤600 nt][CDS, multiple of 3, starts ATG][dn_flank ≤600 nt]`, all 5'→3' on the coding strand. `up_flank` and `dn_flank` headers tell the runtime how many flank nt are present (may be < 600 if neighbour-shortened).

### 1.4 Genome upload to R2

`scripts/upload_genomes.py`:

```python
def upload_genomes(genome_dir: Path):
    client = boto3.client("s3", endpoint_url=f"https://{os.environ['R2_ACCOUNT_ID']}.r2.cloudflarestorage.com",
                          aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"],
                          aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"])
    bucket = os.environ["R2_BUCKET_NAME"]
    manifest = []
    for fna in genome_dir.glob("*.fna"):
        with open(fna, "rb") as f: data = f.read()
        sha = hashlib.sha256(data).hexdigest()
        client.put_object(Bucket=bucket, Key=fna.name, Body=data)
        manifest.append({"isolate_id": fna.stem, "filename": fna.name,
                         "sha256": sha, "size_bytes": len(data),
                         "n_contigs": data.count(b">")})
    write_manifest(manifest)
```

Manifest is committed to git; bucket contents are not.

### 1.5 Adding more isolates later

```bash
# 1. Drop new genome files into the genome dir
cp NEW_LB200.fna NEW_LB201.fna data/genomes_input/

# 2. Run the upload script (only new ones uploaded)
python scripts/upload_genomes.py --only LB200 LB201

# 3. Re-extract gene records (only for these isolates; the gene .fasta files already contain entries for them)
python scripts/build_gene_records.py --isolates LB200 LB201

# 4. Commit data/primer_design/manifest.json + new files in data/primer_design/genes/{lasB,lasR,lasI}/LB20{0,1}.fasta
git add data/primer_design/ && git commit
```

For **adding a new gene** (e.g. `pqsR`):
1. Provide `pqsR.fasta` (multi-record, one per isolate)
2. Run `python scripts/build_gene_records.py --gene pqsR`
3. Add `"pqsR"` to `lib/primer_design/config.py:SUPPORTED_GENES`
4. Frontend gene-dropdown picks it up automatically (it reads from the same constant)

---

## 2 — Vector preparation

### 2.1 Vector loading

```python
def get_vector(name: Literal["pEXG2", "pBBR1MCS2"]) -> VectorRecord:
    gb_path = DATA_DIR / "vectors" / f"{name}.gb"
    record = next(SeqIO.parse(gb_path, "genbank"))
    return VectorRecord(
        name=name,
        sequence=str(record.seq).upper(),
        length=len(record.seq),
        features=record.features,             # parsed MCS, sacB, ori, KanR/GmR
    )
```

### 2.2 Cut site detection

```python
RESTRICTION_SITES = {
    "HindIII":  ("AAGCTT", 1),     # nick at A^AGCTT (between positions 0 and 1 of motif)
    "EcoRI":    ("GAATTC", 1),
    "BamHI":    ("GGATCC", 1),
    "SalI":     ("GTCGAC", 1),
    "XhoI":     ("CTCGAG", 1),
    "KpnI":     ("GGTACC", 5),     # nick at GGTAC^C
    "SacI":     ("GAGCTC", 5),
    "SpeI":     ("ACTAGT", 1),
    "XbaI":     ("TCTAGA", 1),
    "ApaI":     ("GGGCCC", 5),
    "ClaI":     ("ATCGAT", 2),
    "EcoRV":    ("GATATC", 3),     # blunt
    "SmaI":     ("CCCGGG", 3),     # blunt
}

def find_cut_position(vector: VectorRecord, enzyme: str) -> int:
    motif, nick_offset = RESTRICTION_SITES[enzyme]
    fwd_hits = list(find_all(vector.sequence, motif))
    rc_motif = rc(motif)
    rev_hits = [] if motif == rc_motif else list(find_all(vector.sequence, rc_motif))
    total = len(fwd_hits) + len(rev_hits)
    if total != 1:
        raise NotSingleCutter(f"{enzyme} cuts {vector.name} {total} times; need exactly 1")
    site_start = fwd_hits[0] if fwd_hits else rev_hits[0]
    return site_start + nick_offset    # 0-based position of the top-strand nick
```

`get_unique_cutters(vector)` enumerates the dropdown contents — used by frontend at page-load.

### 2.3 Tail derivation

Tail rules per `(vector, enzyme, application)` triple, from `tail_conventions.json`:

```python
TAIL_RULES = {
    "left_arm_15nt_excl_recognition": lambda v, motif, nick: v.sequence[nick - 15 : nick],
    "rc_right_arm_15nt_starting_at_second_nt_of_recognition":
        lambda v, motif, nick: rc(v.sequence[nick : nick + 15]),
    "rc_right_arm_15nt_incl_5nt_recognition":
        lambda v, motif, nick: rc(v.sequence[nick : nick + 15]),    # same span as above; convention name distinguishes intent
    "left_arm_15nt_incl_5nt_recognition":
        lambda v, motif, nick: v.sequence[nick - 10 : nick + 5],    # spans cut: 10 nt before + 5 nt of recognition
}

def derive_tails(vector, enzyme, application) -> tuple[str, str]:
    nick = find_cut_position(vector, enzyme)
    motif, _ = RESTRICTION_SITES[enzyme]
    convention = get_tail_convention(vector.name, enzyme, application)
    p1_tail = TAIL_RULES[convention.p1_tail_rule](vector, motif, nick)
    p4_tail = TAIL_RULES[convention.p4_tail_rule](vector, motif, nick)
    return p1_tail, p4_tail
```

**Calibration check (unit test):** for the three calibrated triples, the derivation must reproduce empirical primer tails:

| Triple | Expected P1 tail | Expected P4 tail |
|---|---|---|
| pEXG2+HindIII+deletion | `CATAAATGTAAAGCA` | `CGACCTGCAGAAGCT` (= RC of `AGCTTCTGCAGGTCG`) |
| pEXG2+HindIII+tagging | `CATAAATGTAAAGCA` | `CGACCTGCAGAAGCT` |
| pBBR1MCS2+HindIII+expression | `CGGTATCGATAAGCT` | `ATTCGATATCAAGCT` |

If any deviation: tail rule needs revision before any further work.

---

## 3 — Gene record retrieval

```python
def get_gene_record(isolate_id: str, gene: str) -> GeneRecord:
    path = DATA_DIR / "genes" / gene / f"{isolate_id}.fasta"
    if not path.exists():
        raise GeneRecordNotFound(f"{gene} not available for {isolate_id} — check curated dataset")
    record = next(SeqIO.parse(path, "fasta"))
    meta = parse_header_metadata(record.description)
    sequence = str(record.seq).upper()
    up = int(meta["up"])
    dn = int(meta["dn"])
    cds = sequence[up : len(sequence) - dn]
    return GeneRecord(
        isolate_id=isolate_id, gene=gene,
        cds_seq=cds,
        up_flank=sequence[:up],
        dn_flank=sequence[len(sequence)-dn:],
        contig_id=meta["contig"],
        genome_start_1based=int(meta["start"]),
        genome_end_1based=int(meta["end"]),
        original_strand=meta["strand"],
    )
```

Validation invariants (asserted on load):
- `len(cds) % 3 == 0`
- `cds.startswith("ATG")`
- `Seq(cds).translate().count("*") == 1`
- `Seq(cds).translate().endswith("*")`

---

## 4 — Application: deletion (pEXG2)

### 4.1 Scar design — exhaustive (N, C) search

For an in-frame deletion that retains N N-terminal codons and C C-terminal codons + stop:

```python
def search_scars(cds: str, scar_max: int = 24) -> list[tuple[int, int]]:
    """Yield all (N, C) pairs with N >= 1, C >= 1, N + C <= scar_max."""
    pairs = []
    for total in range(2, scar_max + 1):           # smallest scar = 2 aa total
        for N in range(1, total):
            C = total - N
            pairs.append((N, C))
    return pairs

def scar_orf(cds: str, N: int, C: int) -> str:
    """Return retained N codons + retained C codons + stop, as nucleotides."""
    last_codon_idx = len(cds) // 3 - 1
    return cds[:3*N] + cds[3*(last_codon_idx + 1 - C):]    # last C codons including stop

def assert_scar_valid(scar: str, N: int, C: int):
    prot = Seq(scar).translate()
    assert prot.startswith("M"), f"scar {N}+{C} doesn't start with M"
    assert prot.endswith("*"), f"scar {N}+{C} doesn't end with stop"
    assert prot.count("*") == 1, f"scar {N}+{C} has internal stop"
    assert len(scar) == (N + C + 1) * 3            # +1 for the stop codon (counted in C)
```

Wait — the stop codon is already part of the last codon of the CDS in our convention (CDS includes stop). So `cds[3*(L-C):]` where `L = len(cds)//3` includes the stop. The +1 is implicit. Re-verify:

```python
def scar_orf(cds: str, N: int, C: int) -> str:
    L = len(cds) // 3              # number of codons including stop
    return cds[:3*N] + cds[3*(L-C):]    # take N from N-term + C from C-term (including stop)

# For LB001 lasB (498 aa + stop = 499 codons, len(cds) = 1497):
# scar(N=16, C=2) = cds[:48] + cds[3*(499-2):] = cds[:48] + cds[1491:1497]
#                 = first 16 codons + last 2 codons (AL + stop)
#                 = "ATGAAGAAG..." + "GCGTTGTAA"   (= MKKVSTLDLLFVAIMG + AL*)
```

### 4.2 Primer construction

Primer geometry, all positions 0-based on top strand of the canonical gene record:

```
P1 (UP-Fwd, + strand): tail_p1 + body_p1
   body_p1 = up_flank[offset_p1 : offset_p1 + L_p1]          for offset_p1 in [0, 150], L_p1 in [18, 28]
   Anchored at 5' end of UP fragment.

P2 (UP-Rev, − strand): tail_p2 + body_p2
   tail_p2 = rc( first 15 nt of DN_segment )
            = rc( cds[3*(L-C) : 3*(L-C) + 15] )       — the C-terminal scar piece + maybe overrun into DN flank
   body_p2 = rc( last L_p2 nt of UP_segment )
            = rc( cds[:3*N][-L_p2:] )                  — anchored at the 3' end of UP fragment
   for L_p2 in [18, 28]

P3 (DN-Fwd, + strand): tail_p3 + body_p3
   tail_p3 = last 15 nt of UP_segment
            = (up_flank + cds[:3*N])[-15:]
   body_p3 = first L_p3 nt of DN_segment
            = (cds[3*(L-C):] + dn_flank)[:L_p3]        — anchored at the 5' end of DN fragment
   for L_p3 in [18, 28]

P4 (DN-Rev, − strand): tail_p4 + body_p4
   body_p4 = rc( dn_flank[end - L_p4 - offset_p4 : end - offset_p4] )
            for offset_p4 in [0, 150], L_p4 in [18, 28]
   Anchored at 3' end of DN fragment.
```

**Candidate-pool independence (key optimization):**
- `(P1 candidates)` depend only on `up_flank` → precompute pool once.
- `(P4 candidates)` depend only on `dn_flank` → precompute pool once.
- `(P2 candidates)` depend only on N (UP-end position).
- `(P3 candidates)` depend only on C (DN-start position).

Total search space: ~300 (N,C) pairs × ~10 (P2) × ~10 (P3) × O(log n) lookup of best (P1, P4) given target Tm. Milliseconds in Python.

### 4.3 Selection algorithm

```python
def search_deletion_primers(gene: GeneRecord, vector: VectorRecord,
                            convention: TailConvention, config: PrimerConfig) -> DeletionPrimerSet:
    p1_tail, p4_tail = derive_tails(vector, convention.enzyme, "deletion")
    cds = gene.cds_seq

    # Precompute P1 and P4 pools, sorted by Tm
    p1_pool = sorted(
        gen_p1_candidates(gene.up_flank, p1_tail, config),
        key=lambda c: c.tm
    )
    p4_pool = sorted(
        gen_p4_candidates(gene.dn_flank, p4_tail, config),
        key=lambda c: c.tm
    )
    if not p1_pool: raise NoCandidates("P1 — relax up-flank GC range or offset")
    if not p4_pool: raise NoCandidates("P4 — relax dn-flank GC range or offset")

    best = None
    top5 = []

    for N, C in search_scars(cds, scar_max=config.scar_max):
        scar = scar_orf(cds, N, C)
        try: assert_scar_valid(scar, N, C)
        except AssertionError: continue                # skip frame-broken scars

        up_segment = gene.up_flank + cds[:3*N]                  # what UP amplicon will contain
        dn_segment = cds[3*(len(cds)//3 - C):] + gene.dn_flank  # what DN amplicon will contain

        p2_pool = list(gen_p2_candidates(up_segment, dn_segment, config))
        p3_pool = list(gen_p3_candidates(up_segment, dn_segment, config))
        if not p2_pool or not p3_pool: continue

        for p2 in p2_pool:
            for p3 in p3_pool:
                target_tm = (p2.tm + p3.tm) / 2
                p1 = closest_by_tm(p1_pool, target_tm)
                p4 = closest_by_tm(p4_pool, target_tm)
                if p1 is None or p4 is None: continue
                tms = [p1.tm, p2.tm, p3.tm, p4.tm]
                spread = max(tms) - min(tms)
                if spread > config.tm_spread_hard_limit: continue        # D5.3
                score = compute_score(tms, [p1.gc, p2.gc, p3.gc, p4.gc])
                candidate = DeletionPrimerSet(
                    N=N, C=C, scar=scar,
                    p1=p1, p2=p2, p3=p3, p4=p4,
                    score=score, spread=spread,
                )
                top5 = update_top5(top5, candidate)
                if best is None or score < best.score:
                    best = candidate

    if best is None:
        raise NoCandidates("no (N,C,P1,P2,P3,P4) tuple passed all hard filters")
    return best, top5
```

### 4.4 Worked example: LB001 ΔlasB pEXG2 HindIII (must reproduce v1)

Input: `gene = lasB` for `isolate = LB001`; vector = pEXG2; enzyme = HindIII; application = deletion.

Expected output:
- N = 16, C = 2 (scar = `MKKVSTLDLLFVAIMGAL*`)
- P1 = `CATAAATGTAAAGCAGGTGTTCCAGCTGGTGCAG`
- P2 = `CCGAGCTTACAACGCACCCATGATCGCAACGAACAAC`
- P3 = `GTTGCGATCATGGGTGCGTTGTAAGCTCGGTGGTC`
- P4 = `CGACCTGCAGAAGCTGCCAGGTACTCGCCTTGC`
- UP amplicon length = 593 bp, DN amplicon length = 538 bp, insert after In-Fusion = 1101 bp, final plasmid = 6155 bp.

Algorithm output may differ in body shift by ±2 nt (as long as filters pass) but `(N, C)` and tail tabs must match exactly.

---

## 5 — Application: tagging (pEXG2 in-locus C-term)

### 5.1 Tag cassette construction

```python
def build_tag_cassette(tag: Tag, position: Literal["C"]) -> str:
    """For in-locus C-term: GGS-linker + tag DNA + new stop codon."""
    assert position == "C", "in-locus tagging is C-terminal only"
    if not tag.in_locus_ok:
        raise TagTooLongForInLocus(f"{tag.name} cassette {tag.cassette_nt} nt > 36 nt limit; use plasmid expression")
    return LINKER_GGS + tag.dna + "TAA"
```

For His6: `GGCGGCAGC` + `CACCATCATCATCACCAC` + `TAA` = 30 nt
For FLAG, His8: cassette = 36 nt
3xFLAG (90 nt), HiBiT (45 nt) → rejected with suggestion to use plasmid expression.

### 5.2 Junction overlap (extended)

Standard In-Fusion junction = 30 nt total (15 nt P2-tail + 15 nt P3-tail, overlapping). For tagging, the junction encodes the full tag cassette:

```
junction_total_len = len(tag_cassette)   # 30 or 36
overlap_len = junction_total_len // 2 + (junction_total_len % 2)  # 15 (His6) or 18 (FLAG, His8)
```

Specifically, the 30-nt His6 cassette splits 15+15 as in deletion. The 36-nt FLAG/His8 cassette splits 18+18 — both primer tails are 3 nt longer. This is the only deviation from deletion-mode.

### 5.3 Primer construction — tagging

UP fragment ends at: native CDS last codon BEFORE stop = `cds[: -3]` (the `-3` excludes native stop).
DN fragment starts at: native sequence right after native stop = `dn_flank[0:]`.

The cassette `tag_cassette = GGS + tag.dna + new_stop` is engineered between them:

```
final insert (after junction collapse):
[up_flank][cds without native stop][GGS-linker][tag DNA][new stop = TAA][dn_flank]
                                    └─────── tag cassette ────────────┘
```

```python
def search_tagging_primers(gene, vector, convention, tag, config):
    cassette = build_tag_cassette(tag, "C")
    p1_tail, p4_tail = derive_tails(vector, convention.enzyme, "tagging")
    cds_no_stop = gene.cds_seq[:-3]               # drop native stop codon

    overlap = len(cassette) // 2 + len(cassette) % 2     # 15 or 18

    # P2 tail = RC of last (overlap) nt of cassette + first 0 nt of dn_flank — but actually we want
    # P2 tail = RC of [first (overlap) nt of cassette] so that P2-tail RC + P3-tail = cassette
    # No — let's be careful. The junction in the assembled insert is:
    #   ...cds_no_stop[-X:][cassette[0:overlap_left]][cassette[overlap_left:]][dn_flank]
    # UP amplicon ends with: cds_no_stop[-some:] + cassette[:overlap_left]
    # DN amplicon starts with: cassette[overlap_left:] + dn_flank[:body_len]
    # The 15- (or 18-)nt overlap region between UP and DN must be identical → easiest: overlap = cassette[overlap_start : overlap_start + overlap_len].
    # We choose overlap to span the middle of the cassette.

    overlap_start = (len(cassette) - overlap) // 2
    overlap_end = overlap_start + overlap

    # P2 (UP-Rev): tail = RC of (cassette[0 : overlap_end])  ← brings the LEFT half of cassette + overlap
    #              body = RC of last L nt of cds_no_stop
    p2_tail = rc(cassette[0:overlap_end])

    # P3 (DN-Fwd): tail = cassette[overlap_start : ]  ← brings the RIGHT half of cassette + overlap
    #              body = first L nt of dn_flank
    p3_tail = cassette[overlap_start:]

    # P1 (UP-Fwd): same as deletion mode (tail_p1 + body anchored 5' of up_flank)
    # P4 (DN-Rev): same as deletion mode (tail_p4 + body anchored 3' of dn_flank)

    # Search/scoring identical to deletion mode, except:
    # — there is no (N, C) loop (junction is fixed by cassette)
    # — Tm constraints apply only to bodies, not tails (tails are fixed by cassette)
    return search_with_fixed_junction(p1_tail, p2_tail, p3_tail, p4_tail,
                                      gene, cassette, config)
```

### 5.4 Length validation

In-locus mode rejected if `tag.cassette_nt > 36`:

```python
if not tag.in_locus_ok:
    return DesignError(
        code="tag_too_long_for_in_locus",
        message=f"{tag.name} cassette ({tag.cassette_nt} nt) exceeds 36 nt limit for in-locus tagging.",
        suggestion="Use plasmid expression with C-terminal tag instead (action='express', tag_position='C')."
    )
```

UI surfaces this before submit (gates the Tag dropdown). API gives same error if it sneaks through.

---

## 6 — Application: expression (pBBR1MCS2)

### 6.1 RBS and spacer construction

Per D6.1.v2, fixed and not user-configurable:

```python
RBS_TAIL_PART = CANONICAL_RBS + SYNTHETIC_SPACER_DEFAULT     # "AGGAGG" + "ACTTGTTC" = 14 nt
```

For the forward primer (P1 in expression mode), full tail = `vector_p1_tail (15 nt) + RBS_TAIL_PART (14 nt)` = 29 nt total.

### 6.2 Primer construction — expression

**Hard invariant (D6.2.v3, 2026-05-11):** the PCR template is always **native
genomic DNA** of the requested isolate. It never contains a tag. Tag-encoding
nucleotides must therefore ride along inside the primer 5' overhangs; only the
body anneals to the genome. Earlier versions of this section anchored the body
inside the engineered tag cassette, which collapsed the genome-annealing
portion to ~0 nt and would not have primed at all on the bench.

```python
def search_expression_primers(gene, vector, convention, tag, tag_position, config):
    p1_vector_tail, p2_vector_tail = derive_tails(vector, convention.enzyme, "expression")

    # Final ORF (used for reporting, verification, plasmid assembly).
    if tag is None:
        coding_seq = gene.cds_seq
    elif tag_position == "N":
        if not tag.n_term_ok:
            raise InvalidTagPosition(f"{tag.name} cannot be N-terminal")
        coding_seq = "ATG" + tag.dna + LINKER_GGS + gene.cds_seq[3:]
    elif tag_position == "C":
        if not tag.c_term_ok:
            raise InvalidTagPosition(f"{tag.name} cannot be C-terminal")
        coding_seq = gene.cds_seq[:-3] + LINKER_GGS + tag.dna + "TAA"

    # What the primer bodies actually anneal to (native genomic DNA), and what
    # extra non-templated DNA we tack onto each tail to introduce the tag.
    if tag is None:
        native_template = gene.cds_seq
        p1_tag_overhang = ""
        p2_tag_overhang_rc = ""
    elif tag_position == "N":
        native_template = gene.cds_seq[3:]                                # drop native ATG
        p1_tag_overhang = "ATG" + tag.dna + LINKER_GGS                    # 5'→3' coding orientation
        p2_tag_overhang_rc = ""
    elif tag_position == "C":
        native_template = gene.cds_seq[:-3]                               # drop native stop
        p1_tag_overhang = ""
        p2_tag_overhang_rc = rc(LINKER_GGS + tag.dna + "TAA")

    p1_tail = p1_vector_tail + RBS_TAIL_PART + p1_tag_overhang
    p2_tail = p2_vector_tail + p2_tag_overhang_rc

    p1_pool = gen_body_candidates_at_5end(native_template, p1_tail, config)
    p2_pool = gen_body_candidates_at_3end(native_template, p2_tail, config)

    # Search: minimize Tm spread. Amplicon assembled downstream as
    #     amplicon = p1.tail + native_template + rc(p2.tail)
    # which collapses to p1_vector_tail + RBS + spacer + coding_seq + rc(p2_vector_tail).
    return select_best_pair(p1_pool, p2_pool, config)
```

For a C-terminal tag the P2 tail thus has three regions, 5'→3':
`vector_tail (15) | rc(new_stop) (3) | rc(tag_dna) (n) | rc(GGS_linker) (9)`,
followed by the native-annealing body. The tag overhang reverse-complements
the GGS linker into a GC-rich stretch immediately upstream of the body, so
the relaxed `gc5 ≤ 4` filter tier is commonly invoked on the reverse primer
of C-term-tagged designs — this is expected and acceptable.

### 6.3 Worked example: PA14 lasR pBBR1MCS2 HindIII (must reproduce empirical #1612 / #1608)

Input: gene = lasR for PA14; vector = pBBR1MCS2; enzyme = HindIII; application = expression; tag = None.

Expected:
- P1 starts: `cggtatcgataagctaggaggacttgttc` (29 nt tail) + `ATGGCCTTGGTT...` (body anchored at ATG, 22 nt)
- P2: `attcgatatcaagct` (15 nt tail) + `TCAGAGAGTAATAAGACCCAAATT` (body, RC of stop region, 24 nt)
- PCR amplicon ≈ 750 bp

Algorithm acceptance: P1 and P2 sequences within ±2 nt of empirical (body length tolerance), with full tail sequences exact.

---

## 7 — Primer body filters (shared)

Applied to every candidate body before scoring:

```python
def passes_hard_filters(body: str, config: PrimerConfig) -> bool:
    if not (config.body_len_min <= len(body) <= config.body_len_max):    return False
    if body[-1] not in "GC":                                              return False     # 3' GC clamp
    gc5 = sum(b in "GC" for b in body[-5:])
    if not (1 <= gc5 <= 3):                                               return False     # last-5 GC count
    if len(set(body[-4:])) == 1:                                          return False     # 4-identical at 3' end
    if any(body[i] == body[i+1] == body[i+2] == body[i+3] for i in range(len(body)-3)):
        return False                                                                       # any 4-homopolymer
    gc = sum(b in "GC" for b in body) / len(body)
    if not (config.gc_min <= gc <= config.gc_max):                        return False     # 0.40–0.70 default
    if max_3prime_self_dimer(body) >= 5:                                  return False     # ≥ 5-nt 3' self-dimer
    tm = Tm_NN(Seq(body), Na=50, dnac1=500, dnac2=0)                                       # D5.1
    if not (config.tm_min <= tm <= config.tm_max):                        return False     # 58–66 °C
    return True
```

`max_3prime_self_dimer(body)`: slide RC of body against itself in 8-nt windows, count contiguous 3'-end matches.

---

## 8 — Tm calculation

```python
from Bio.SeqUtils.MeltingTemp import Tm_NN

def primer_body_tm(body: str) -> float:
    return Tm_NN(Seq(body), Na=50, dnac1=500, dnac2=0)
```

Allawi & SantaLucia 1997 nearest-neighbor parameters with Na+ correction (SantaLucia 1998), as implemented by Biopython. This gives ±0.5 °C accuracy for primers in 18–28 nt range under In-Fusion buffer conditions (D5.1).

**Annealing temperature recommendation in PDF** (corrected 2026-07-27 — previously
omitted the base −5 °C, so B7 showed Ta = Tm_body with no offset):
```python
T_anneal_C = min(p.tm for p in primers) - 5.0 + POLYMERASE_OFFSETS_C[polymerase]
```

Base rule of thumb is Tm − 5 °C; POLYMERASE_OFFSETS_C then adjusts further per
enzyme. For B7 the offset is 0; Phusion / Q5 +3; Taq −5 (D5.2).

PDF additionally suggests gradient ±5 °C around the recommended Ta if first attempt fails.

---

## 9 — Scoring and selection

```python
def compute_score(tms: list[float], gcs: list[float]) -> float:
    tm_spread = max(tms) - min(tms)
    mean_tm = sum(tms) / len(tms)
    gc_spread = max(gcs) - min(gcs)
    gc_penalty = sum(0.1 * max(0, 0.45 - g) + 0.1 * max(0, g - 0.55) for g in gcs)
    return tm_spread + 0.25 * abs(mean_tm - 62.0) + 0.02 * gc_spread + gc_penalty

# Hard cutoff: tm_spread > 3 °C → reject (D5.3)
# Lower score is better
```

Selection: keep global minimum + top-5 for inspection (PDF "Alternativen"-Box).

---

## 10 — Off-target scan

Per D5.5, scan target isolate's full genome (lazy-fetched from R2):

```python
def off_target_scan(primers: list[Primer], genome_fasta: bytes, config) -> OffTargetReport:
    sites_per_primer = {}
    for primer in primers:
        body = primer.body
        anchor = body[-10:]                               # last 10 nt of body
        body_rc = rc(body)
        anchor_rc = body_rc[-10:]
        sites = []
        for contig_id, contig_seq in parse_fasta(genome_fasta):
            for start in find_with_mismatch(contig_seq, anchor, max_mm=1):    # ≤ 1 mm in last 10 nt
                aligned = contig_seq[start - len(body) + 10 : start + 10]
                if hamming(aligned, body) <= 4:
                    sites.append(Site(contig_id, start, "+", primer.name))
            for start in find_with_mismatch(contig_seq, anchor_rc, max_mm=1):
                aligned = rc(contig_seq[start - len(body_rc) + 10 : start + 10])
                if hamming(aligned, body) <= 4:
                    sites.append(Site(contig_id, start, "-", primer.name))
        sites_per_primer[primer.name] = sites

    products = enumerate_products(sites_per_primer, max_size_bp=6000, min_size_bp=50)
    return OffTargetReport(sites=sites_per_primer, products=products)
```

**Required outcomes** (PCR-pair-wise):

| Application | Pair | Required outcome |
|---|---|---|
| Deletion / Tagging | P1+P2 | exactly the intended UP amplicon |
| Deletion / Tagging | P3+P4 | exactly the intended DN amplicon |
| Deletion / Tagging | P1+P4 | exactly the WT locus (used in colony PCR for verification) |
| Deletion / Tagging | P3+P2 | zero products |
| Expression | P1+P2 | exactly the intended insert |

Any unexpected product = hard block, raise `OffTargetDetected` with details for the PDF report and developer log.

---

## 11 — Mandatory verification

All `D7.4` hard checks; raise `VerificationFailure` on any failure.

```python
def verify_deletion(result: DeletionResult, request: DesignRequest):
    # 1. Recognition site count in final plasmid matches convention
    enz_motif = RESTRICTION_SITES[request.enzyme][0]
    expected = result.convention.expected_recognition_count_in_final_plasmid
    actual = result.final_plasmid.count(enz_motif) + result.final_plasmid.count(rc(enz_motif))
    actual = actual if rc(enz_motif) != enz_motif else result.final_plasmid.count(enz_motif)
    assert actual == expected, f"{request.enzyme} recognition count {actual} != expected {expected}"

    # 2. Off-target scan returned no unintended products
    assert result.off_target.passed, f"off-target violations: {result.off_target.violations}"

    # 3. Body reconstruction matches genome
    for primer in [result.p1, result.p2, result.p3, result.p4]:
        reconstructed = derive_body_from_genome(primer, result.gene_record)
        assert primer.body == reconstructed, f"{primer.name} body drift: search vs reconstructed"

    # 4. Scar ORF valid
    scar_protein = Seq(result.scar).translate()
    assert scar_protein.startswith("M") and scar_protein.endswith("*") and scar_protein.count("*") == 1, \
        f"scar ORF malformed: {scar_protein}"
    assert len(result.scar) == (result.N + result.C) * 3, f"scar length wrong"

    # 5. Junction overlap exact (last 30 nt of UP == first 30 nt of DN)
    assert result.up_amplicon[-30:] == result.dn_amplicon[:30], "junction overlap mismatch"

    # 6. Final plasmid size = vector + insert - 30
    expected_size = len(result.vector.sequence) + len(result.insert) - 30
    assert len(result.final_plasmid) == expected_size, f"plasmid size wrong: {len(result.final_plasmid)} vs {expected_size}"

def verify_tagging(result, request):
    # Same as deletion plus:
    # 7. Assembled CDS-tag fusion translates correctly
    fusion_dna = result.up_segment[:-(result.cassette_overlap_left)] + result.cassette + result.dn_segment[result.cassette_overlap_right:]
    fusion_aa = Seq(extract_orf(fusion_dna)).translate()
    assert fusion_aa.startswith("M") and fusion_aa.endswith("*") and fusion_aa.count("*") == 1
    assert "GGS" in fusion_aa or "GGGS" in fusion_aa            # linker present
    expected_tag_aa = result.tag.protein
    assert expected_tag_aa in fusion_aa, "tag missing from fusion ORF"

def verify_expression(result, request):
    # Same as deletion §1–3, §5–6 (junction is single 15-nt overlap, not 30) plus:
    # 7. Insert in-frame, single stop, RBS present
    insert = result.amplicon
    rbs_pos = insert.find("AGGAGG")
    assert rbs_pos > 0 and rbs_pos < 50, "RBS not at expected position"
    atg_pos = insert.find("ATG", rbs_pos + 6)
    assert 8 <= (atg_pos - rbs_pos) <= 18, "RBS-ATG spacing out of range"
    cds = insert[atg_pos : extract_stop(insert, atg_pos)]
    fusion_aa = Seq(cds).translate()
    assert fusion_aa.startswith("M") and fusion_aa.endswith("*") and fusion_aa.count("*") == 1
    if request.tag:
        assert TAGS[request.tag].protein in fusion_aa, f"{request.tag} not in expressed protein"
```

**Soft warnings** (D7.5) accumulated and surfaced in PDF "Hinweise" but do not block:
- Restriction site present in insert (would prevent re-linearization)
- Tm spread > 2 °C (closer to the hard cutoff than ideal)
- Body GC outside 0.45–0.55 preferred range
- For expression: native stop codon detected within tag-fusion region (only relevant for unusual genes)

---

## 12 — Output assembly

### 12.1 Combined annotated FASTA

Single file, records in this order:

```
>P1_<gene>_UP_Fwd  application=<app>  tail_len=<n>  body_len=<m>  Tm_body=<x>C  GC_body=<y>%  tail_kind=<key>
<full primer sequence, tail in lowercase, body in UPPERCASE>

>P2_<gene>_UP_Rev  ... [if deletion or tagging]
<...>

>P3_<gene>_DN_Fwd  ... [if deletion or tagging]
<...>

>P4_<gene>_DN_Rev  ...
<...>

>UP_amplicon  len=<n>  ... [deletion / tagging]
<60 chars per line>

>DN_amplicon  len=<n>  ... [deletion / tagging]
<60 chars per line>

>Insert_after_InFusion  len=<n>  joining_collapsed=true  ... [deletion / tagging]
<60 chars per line>

>Final_plasmid_<vector>_<modification>  circular=true  len=<n>  validation=passed
<60 chars per line>

>Scar_or_Cassette_ORF  application=<app>  translation=<protein-string>
<DNA sequence of scar or fusion CDS region>
```

### 12.2 PDF templates

Three templates (`deletion`, `tagging`, `expression`), structurally:

| Section | Deletion | Tagging | Expression |
|---|---|---|---|
| 1. Konstrukt-Überblick | Y (with scar) | Y (with cassette) | Y (with RBS+tag) |
| 2. Primer-Tabelle | 4 primers | 4 primers | 2 primers |
| 3. Off-Target-Check | 4 pairs | 4 pairs | 1 pair |
| 4. PCR-Bedingungen | per polymerase | per polymerase | per polymerase |
| 5. Vektor-Vorbereitung + In-Fusion | 3-fragment | 3-fragment with extended junction | 2-fragment |
| 6. Cloning workflow | conjugation + sucrose CS + Gm sensitivity screen | same | transformation + Km selection (no sucrose) |
| 7. Validation suggestions | Sanger both junctions + colony PCR | Sanger junction + Western (anti-tag) | Sanger insert + Western (anti-tag) + activity assay |
| 8. Ordering checklist | yes | yes | yes |
| 9. Hinweise (soft warnings) | yes | yes | yes |

### 12.3 JSON manifest

Schema in `project_plan.md` §2.7. All fields are filled; missing optionals are explicit `null`.

---

## 13 — Failure modes and recovery hints

| Failure | Likely cause | Suggested fix |
|---|---|---|
| `NoCandidates: P1` | up-flank GC too high or homopolymer-rich | relax `gc_max` to 0.72; widen `offset_p1` to 0–200 |
| `NoCandidates: P4` | same on dn-flank | as above |
| Tm spread > 3 °C across all (N,C) | one anchor stuck in GC-rich patch | shift `(N, C)` ±1–3 aa; if persists, relax body length to 16–30 |
| Off-target product 50–6000 bp | repeat element near gene | increase C asymmetry to push P3 out of repeat; or pick different scar (N, C) |
| `TagTooLongForInLocus` | 3xFLAG, HiBiT in-locus mode | switch to plasmid expression with C-term tag |
| Verification: AAGCTT count off | tail derivation mismatched convention | review `tail_conventions.json` calibration; check enzyme is single-cutter |
| Verification: junction mismatch | algorithm bug, P2 tail and P3 tail not derived from same scar | bug; report with full primer set + (N,C) |
| sacB-tolerant clinical isolate (lab observation) | some biofilm-adapted PA strains | use I-SceI counterselection or Gm-loss screen; not a tool issue |

---

## 14 — Host-specific defaults

The tool defaults to *P. aeruginosa* settings (GC 66 %). For other hosts, adjust `config.PrimerConfig`:

| Host | Genome GC | gc_min | gc_max | offset_max | Notes |
|---|---|---|---|---|---|
| *E. coli* | 51 % | 0.40 | 0.60 | 100 | tightest defaults |
| **P. aeruginosa** (default) | 66 % | 0.40 | 0.70 | 150 | relaxed upper |
| *Burkholderia* | 67 % | 0.40 | 0.72 | 150 | very similar to PA |
| *Streptomyces* | 72 % | 0.45 | 0.75 | 200 | smaller candidate pools |
| *S. aureus* | 33 % | 0.30 | 0.55 | 100 | lower clamp; many candidates |

In v1, only *P. aeruginosa* host is supported (matches the curated dataset). Other hosts deferred.

---

## 15 — Worked acceptance examples (regression suite)

These three runs are the canonical acceptance tests for Phase 2. The tool must reproduce the exact tail tabs and within ±2 nt body shift the empirical primers; final plasmid SHA-256 must match a hand-validated reference.

### 15.1 Deletion: LB001 ΔlasB pEXG2 HindIII

```python
request = {"isolate_id": "LB001", "gene": "lasB", "action": "delete",
           "enzyme": "HindIII", "polymerase": "B7"}
expected = {
    "N": 16, "C": 2,
    "scar_protein": "MKKVSTLDLLFVAIMGAL*",
    "p1_tail": "CATAAATGTAAAGCA",
    "p4_tail": "CGACCTGCAGAAGCT",
    "final_plasmid_length_bp": 6155,
    "up_amplicon_length_bp": 593,
    "dn_amplicon_length_bp": 538,
    "insert_length_bp": 1101,
    "AAGCTT_count_in_final_plasmid": 0,
}
```

### 15.2 Expression: PA14 lasR pBBR1MCS2 HindIII (untagged)

```python
request = {"isolate_id": "PA14", "gene": "lasR", "action": "express",
           "enzyme": "HindIII", "polymerase": "B7"}
expected = {
    "p1_tail": "CGGTATCGATAAGCT" + "AGGAGG" + "ACTTGTTC",   # 29 nt
    "p2_tail": "ATTCGATATCAAGCT",                            # 15 nt
    "p1_body_starts": "ATGGCCTTGGTT",                        # ATG of lasR
    "p2_body_starts": "TCAGAGAGTAATAAGACCCAAATT",            # RC of lasR stop region
    "amplicon_length_bp": 750,
    "AAGCTT_count_in_final_plasmid": 1,                      # site_partial_AAGCT regenerates one site
}
```

### 15.3 Tagging: PA14 lasR pEXG2 HindIII His6 (in-locus C-term)

```python
request = {"isolate_id": "PA14", "gene": "lasR", "action": "tag",
           "tag": "His6", "tag_position": "C", "use_plasmid_for_tag": False,
           "enzyme": "HindIII", "polymerase": "B7"}
expected = {
    "tag_cassette": "GGCGGCAGC" + "CACCATCATCATCACCAC" + "TAA",       # 30 nt
    "junction_overlap": 15,                                            # split 15+15 (His6 fits exactly)
    "fusion_protein_endswith": "GGSHHHHHH*",
    "AAGCTT_count_in_final_plasmid": 0,
}
```

### 15.4 Tag-too-long rejection

```python
request = {"isolate_id": "PA14", "gene": "lasR", "action": "tag",
           "tag": "HiBiT", "tag_position": "C",
           "enzyme": "HindIII", "polymerase": "B7"}
# Expected: DesignError(code="tag_too_long_for_in_locus",
#                       suggestion="Use plasmid expression instead")
```

---

## 16 — Glossary

- **In-Fusion** — Takara's homology-based DNA assembly using 5' exonuclease + ligase. Requires 15 nt of identical sequence at fragment junctions.
- **Junction overlap** — the 15-nt (or longer) region of identical sequence that two adjacent fragments share. After assembly, this region appears once in the final construct.
- **Scar** — for in-frame deletion: the protein sequence remaining after the deletion. N N-terminal aa + C C-terminal aa + stop.
- **Tail** — the non-templated 5' portion of a primer that mediates In-Fusion homology. Vector-tails join to vector arms; junction-tails join two PCR fragments.
- **Body** — the 3' portion of a primer that anneals to template during PCR. Subject to all hard filters (§7).
- **Cut nick position** — 0-based index on top strand where a restriction enzyme cuts. For HindIII A^AGCTT cut at recognition position p (0-based), nick = p + 1.
- **CDS** — coding sequence; from ATG through stop codon, inclusive.

---

## 17 — Change log

| Version | Date | Change |
|---|---|---|
| 1.0 | 2026-04 | Original LB001-ΔlasB-pEXG2 demo skill, deletion-only |
| 2.0 | 2026-05-04 | Generalized to 3 applications. Added data-ingestion pipeline (.fna + .fasta inputs). Added tagging (in-locus C-term) and expression (pBBR1MCS2) workflows. Tail derivation generalized over (vector × enzyme × application) triples. Three PDF templates. RBS strategy fixed per D6.1.v2. Three worked acceptance examples for regression. |
