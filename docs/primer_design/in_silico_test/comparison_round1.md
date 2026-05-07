# In-silico-Vergleich Runde 1 — Tool vs. unabhängiges Design

Datum: 2026-05-07. Spec-Referenz: `docs/primer_design/decisions_log.md` v0.2.

## Test-Set

| Case | Tool-Output | Claude-Output |
|---|---|---|
| LB014 ΔlasR | ✅ PASS | ✅ designed |
| LB014 ΔlasB | ✅ PASS | ✅ designed |
| LB020 ΔlasB | ✅ PASS | ✅ designed |
| LB056 ΔlasB | ✅ PASS | ✅ designed |
| LB056 ΔlasR | ❌ verification_failure (HindIII=1) | ✅ designed (siehe Diskussion) |
| LB020 ΔlasR | — FASTA fehlt im Repo | — |

## Globale Beobachtungen

### 1. Mein P1-Tail ist falsch (Bug bei mir, nicht beim Tool)

| Quelle | P1-Tail | Mechanik |
|---|---|---|
| Mein Designer | `GCATAAATGTAAAGC` | letzte 15 nt **vor** dem AAGCTT (Spec wörtlich) |
| Tool | `CATAAATGTAAAGCA` | letzte 15 nt des linearisierten Vektor-3'-Endes (inkl. 1. A von AAGCTT, der nach HindIII-Schnitt am linken Arm bleibt) |

HindIII schneidet `A^AGCTT`. Nach dem Schnitt endet der lineare Vektor-Top-Strand auf `…AAATGTAAAGCA` — die letzten 15 nt sind also `CATAAATGTAAAGCA`, nicht `GCATAAATGTAAAGC`. Die Tool-Version ist mechanisch korrekt; die Spec-Formulierung „last 15 nt of left arm before AAGCTT" ist missverständlich. **Tool-Bug? Nein. Spec-Wording-Bug — sollte präzisiert werden.** Ich fixe meinen Designer.

P4-Tail (`CGACCTGCAGAAGCT`) übereinstimmend, kein Diff.

### 2. Scar-Strategie: Tool ist deutlich klüger

| Case | Tool-Scar | Mein Scar |
|---|---|---|
| LB014 ΔlasR | **15 N + 2 C** = MALVDGFLELERSSG\|TL\* | 3 N + 3 C = MALTL\* |
| LB014 ΔlasB | **16 N + 2 C** = MKKVSTLDLLFVAIMG\|AL\* | 3 N + 3 C = MKKAL\* |
| LB020 ΔlasB | 16 N + 2 C (identisch) | 3 N + 3 C |
| LB056 ΔlasB | **15 N + 5 C** = MKKVSTLDLLFVAIM\|CPSAL\* (anders!) | 3 N + 3 C |

→ Das Tool macht eine asymmetrische, oft große (N,C)-Suche und schiebt die Junction in GC-reiche, gut amplifizierbare Regionen hinein. Mein Designer wählt naiv das kleinste (N,C) ≥ (3,3) und musste dann den P2-Body 39 nt nach 5' retrahieren (76 nt Primer!). Der Tool-Ansatz ist korrekt und mein nächster Designer-Iteration anzupassen.

→ **Bei LB056 ΔlasB schiebt das Tool C von 2 auf 5** — eine SNP im (16,2)-Bereich macht offensichtlich die Default-(16,2)-Junction unbrauchbar (Tm/GC oder Sekundärstruktur), und der Optimizer findet (15,5). Gut.

### 3. Tm-Optimierung sehr eng

Tool spread: 0.16 / 0.21 / 0.21 / 0.28 °C. Spec erlaubt ≤ 3 °C, Tool zielt anscheinend auf < 0.5 °C — exzellent.

### 4. LB014 vs. LB020 ΔlasB: Primer **identisch**, finale Plasmidgröße 6155 vs. 6156 bp

→ 1-bp-Indel in der DN-Region zwischen den beiden Isolaten, sonst sequenzidentisch. Tool reproduziert dies konsistent. ✓

---

## Root-Cause-Analyse: LB056 ΔlasR `verification_failure: HindIII=1`

**Diagnostik:**
- Im LB056-lasR-Locus selbst kommt **kein** AAGCTT vor (weder Top noch Bottom Strand) ✓
- Keine (N,C)-Scar-Kombination erzeugt AAGCTT in der Insert-Junction ✓
- Mein Designer-Scan der P4-Body-Kandidaten in der LB056-lasR DN-Region zeigt: **die zwei Top-Tm-Kandidaten beginnen beide mit `T`** (`TATCGAGAATTCGCCAGCAACCG`, `TGGTATCGAGAATTCGCCAGCAAC`).

**Mechanismus der Site-Regeneration:**

Im finalen Plasmid endet das Insert (Top-Strand) auf `…RC(P4_body) + AGCTTCTGCAGGTCG` (RC der P4-Tail). Wenn `P4_body` mit `T` **beginnt**, endet `RC(P4_body)` mit `A`. Konkateniert:

```
…RC(P4_body) | AGCTTCTGCAGGTCG…
…X-X-X-X-A   | A-G-C-T-T-…
                ↑↑↑↑↑↑
              AAGCTT regeneriert!
```

→ HindIII-Site re-erstellt am rechten Vektor-Junction.

**Warum LB014 ΔlasR durchläuft, LB056 ΔlasR nicht:** LB014 und LB056 lasR unterscheiden sich in 10 SNPs über UP/CDS/DN. Das Tool wählt Bodies isolatweise neu; bei LB014 fiel die Wahl auf `CCAGCCTTTGCGCTCCTTG` (beginnt mit `C`, kein Trap), bei LB056 vermutlich auf einen `T…`-Body wegen marginaler Tm-Unterschiede.

**Bug-Klassifikation:** Algorithmus-Bug, **nicht** Verifikations-Bug. Verifikation funktioniert korrekt (fängt's ab). Der Bug ist in `primer_designer.py` (oder gleichwertig): Body-Filter prüft nicht, ob die Kombination `(body_5'_terminal_base, vector_tail_3'_end)` die Cut-Site regeneriert.

## Vorgeschlagene Fixes

### Fix #1 — Hard-Filter für Junction-Site-Regeneration (kritisch)

In Body-Selektion zusätzlich filtern:
- **P1-Body:** darf nicht mit `AGCTT` beginnen, wenn Tail auf `…AAAGCA` endet (HindIII-Convention `site_destroyed`)
- **P4-Body:** darf nicht mit `T` beginnen, wenn Tail auf `…AGCTT…` (im finalen Plasmid) konstruiert wird

Generisch: für jedes Cut-Enzym + Convention das verbotene 5'-Terminal-Pattern für P1/P4 vorberechnen und in `body_passes()` einbauen.

### Fix #2 — Spec-Wording präzisieren (kosmetisch)

`decisions_log.md` D4.3, Zeile zu P1-Tail: „last 15 nt of left arm before AAGCTT" → „last 15 nt of the linearized vector top strand at the 3' end (= 15 nt up to and including the first A of AAGCTT, which remains on the left arm after HindIII cleaves A^AGCTT)". Symmetrisch P4-Tail-Wording sanity-checken.

### Fix #3 — Mein Designer (out of scope, intern)

- P1-Tail auf `CATAAATGTAAAGCA` korrigieren
- Scar-Suche: kleinstes (N,C) durch optimale (N,C) ersetzen, mit Score = Tm-Spread + Body-Quality

---

## Was ich für den nächsten Test brauche

1. **LB020 lasR FASTA** ins Repo (fehlt noch unter `data/primer_design/genes/lasR/`).
2. **Tool-Output für LB056 ΔlasR**: Auch wenn die Verifikation blockt — die generierten Primer + Trace, damit ich das `T…`-Body-Hypothese definitiv bestätigen kann.
3. **Liste der „weiteren Bugs"**, die dir aufgefallen sind — am besten kurz mit Case + erwartetem vs. beobachtetem Verhalten.
