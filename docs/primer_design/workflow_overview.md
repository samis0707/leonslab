# Primer Design Tool — Workflow hinter den Kulissen

Kurzfassung: Was passiert, wenn jemand auf der `/primer-design`-Seite auf
**„Design primers"** klickt. Gedacht zum Erklären in 5 Minuten.

## Big Picture in einem Satz

Browser schickt die Eingaben als JSON an eine Vercel-Serverless-Funktion
(`/api/design`); ein Python-Modul sucht passende Primer, baut das fertige
Plasmid im Speicher zusammen, prüft es gegen das Genom, und schickt
PDF + FASTA + JSON zurück an die Seite zum Download.

## Die Schritte im Detail

### 1. Frontend — Formular & Request

- Datei: `src/views/PrimerDesignView.astro`
- Der User wählt: Anwendung (delete / tag / express), Gen, Isolat, Enzym,
  optional Tag.
- Beim Submit wird ein JSON-Body via `fetch('/api/design', { method: 'POST' })`
  geschickt. Mehr macht der Browser nicht — keine Logik, kein Genom im
  Frontend.

### 2. API Handler — Eingang & Validierung

- Datei: `api/design/handler.py`
- Vercel ruft die `handler`-Klasse als Serverless-Funktion auf.
- `_validate_request` checkt Pflichtfelder und erlaubte Werte
  (action, tag, polymerase …). Bei Verstoß: 400 mit `error_code`.
- Gültige Requests gehen an `dispatch(request)`.

### 3. Dispatch — Welche Anwendung?

- Datei: `lib/primer_design/applications/__init__.py`
- Routet je nach `request.application` zu einem von drei Orchestrierern:
  - `deletion.py` — In-frame Deletion auf pEXG2 (4 Primer)
  - `tagging.py`  — In-locus C-terminal Tag auf pEXG2 (4 Primer + Cassette)
  - `expression.py` — Plasmid-Expression auf pBBR1MCS2 (2 Primer)

### 4. Pipeline pro Anwendung

Für alle drei Anwendungen läuft im Prinzip dieselbe Sequenz:

1. **Gen laden** (`gene_finder.get_gene_record`)
   FASTA aus `data/primer_design/genes/<gene>/<isolate>.fasta` lesen → liefert
   CDS + Up/Down-Flanken.
2. **Vektor + Tail-Convention laden** (`vectors`)
   Welche Schwänze (homologe Überhänge) müssen die Primer für In-Fusion
   haben? Hängt ab von Vektor (pEXG2 / pBBR1MCS2), Enzym (Schnittstelle) und
   Anwendung (deletion / tagging / expression).
3. **Primer-Suche** (`primer_designer.search_*_primers`)
   - Generiert Kandidaten aus den Flanken / CDS-Enden.
   - Berechnet Tm (Allawi & SantaLucia, via Biopython `Tm_NN`).
   - Hard Filter: Länge, GC-Gehalt, 3'-G/C-Clamp, Homopolymer-Run, Primer-Dimer.
   - Score-Funktion wählt das beste Set; Tm-Spread aller Primer ≤ Schwellwert.
4. **Amplicon-Aufbau** (`_build_amplicons`)
   - Setzt UP- und DN-PCR-Produkte rechnerisch zusammen: `Tail + Body + Junction`.
   - Bei Tagging: dazwischen sitzt die Tag-Cassette.
   - Bei Expression: nur ein Amplicon (`P1.tail + CDS + RC(P2.tail)`).
5. **Genom holen** (`storage_adapter.fetch_genome`)
   - Lazy: das Genom liegt **nicht im Repo**, sondern als FASTA in einem
     Cloudflare-R2-Bucket. Erst jetzt wird es geladen (SHA256 gegen das
     Manifest geprüft).
6. **Off-Target-Scan** (`primer_designer.off_target_scan`)
   - Anchor-Walk + Body-Extension gegen das echte Isolat-Genom.
   - Vergleicht gefundene PCR-Produkte mit `expected_products`. Alles
     Unerwartete → `OffTargetDetected` (400).
7. **Plasmid-Assembly** (`plasmid_builder.assemble`)
   - Vektor wird am Enzym-Schnittpunkt aufgemacht, Insert eingefügt.
8. **Verifikation** (`verification.verify`)
   - Harte Checks: keine zusätzlichen Schnittstellen, Plasmid-Größe im
     erlaubten Fenster, Anwendungs-spezifische Frame-Checks (z. B. Scar
     in-frame, Tag-Cassette in-frame, Expression-CDS in-frame & sauber
     gestoppt).

### 5. Artefakte erzeugen

Zurück in `handler.run_design`:

- **FASTA** (`fasta_writer.write_fasta`) — alle Primer + finales Plasmid.
- **JSON-Manifest** (`json_writer.write_manifest`) — vollständiger
  reproduzierbarer Bericht (`tool_version`, `job_id`, alle Primer-Metriken,
  Off-Target-Report, Sequenzen).
- **PDF** (`pdf_writer.render_pdf`) — Bench-Sheet zum Ausdrucken,
  geht erst nach `/tmp`, wird dann gelesen.

PDF und FASTA werden **base64-codiert** in die JSON-Antwort eingebettet,
damit sie über eine einzige HTTP-Antwort zurückkommen.

### 6. Antwort & Anzeige

- Server schickt `{ status, job_id, artefacts: { fasta, pdf, json }, summary, warnings }`.
- Frontend rendert die Primer-Tabelle aus `summary.primers` und baut drei
  Download-Buttons aus den base64-Strings (`PDF`, `FASTA`, `JSON`).

## Mentales Bild zum Erklären

```
Browser                Vercel Function           lib/primer_design/
──────────             ───────────────           ──────────────────
Formular  ──POST──►   handler.py                  applications/{deletion,
   ▲                   │  validate                              tagging,
   │                   │  dispatch ───────────►   expression}.run()
   │                   │                            │
   │                   │                            ├─ gene_finder    (FASTA aus data/)
   │                   │                            ├─ vectors        (Tail-Convention)
   │                   │                            ├─ primer_designer (Tm, Filter, Score)
   │                   │                            ├─ storage_adapter (Genom aus R2)
   │                   │                            ├─ off_target_scan
   │                   │                            ├─ plasmid_builder
   │                   │                            └─ verification
   │                   │  write_fasta + write_manifest + render_pdf
   └──JSON+base64──────┘
```

## Was wo liegt — Quick Reference

| Zweck                       | Datei                                         |
| --------------------------- | --------------------------------------------- |
| UI-Formular                 | `src/views/PrimerDesignView.astro`            |
| HTTP-Eingang + Validierung  | `api/design/handler.py`                       |
| Anwendungs-Dispatch         | `lib/primer_design/applications/__init__.py`  |
| Pipeline pro Anwendung      | `lib/primer_design/applications/*.py`         |
| Primer-Algorithmus          | `lib/primer_design/primer_designer.py`        |
| Vektoren + Tail-Conventions | `lib/primer_design/vectors.py`                |
| Tags / Cassettes            | `lib/primer_design/tags.py`                   |
| Gen-Records                 | `lib/primer_design/gene_finder.py` + `data/`  |
| Genome aus R2               | `lib/primer_design/storage_adapter.py`        |
| Plasmid-Assembly            | `lib/primer_design/plasmid_builder.py`        |
| Endkontrolle                | `lib/primer_design/verification.py`           |
| Artefakte (FASTA/JSON/PDF)  | `fasta_writer.py`, `json_writer.py`, `pdf_writer.py` |
