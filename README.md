# LeonsLab

A creative-space portfolio for Leon Rösch — built with Astro, deployed on Vercel.

## Local development

```bash
npm install
npm run dev      # → http://localhost:4321
npm run build    # → static output to ./dist
```

## Project layout

```
src/
  components/   Reusable UI (Header, Footer, Hero, Cartoon/Sticker blocks…)
  views/        Page-level views composed from components
  pages/        Astro file-based routes (DE at root, EN under /en/)
  layouts/      BaseLayout
  content/      Astro content collections (cartoons, stickers)
  i18n/         Bilingual UI dictionary + helpers
  styles/       Tailwind v4 + global tokens
public/         Static assets (optimised images live here)
```

## Primer Design Tool

A sub-page under `/primer-design` ships an In-Fusion primer-design tool for
*P. aeruginosa* cloning experiments. The tool covers three applications:

1. **In-frame deletion** on the suicide vector pEXG2 — 4 primers, 3-fragment
   In-Fusion assembly.
2. **In-locus C-terminal tagging** on pEXG2 — 4 primers, 3-fragment assembly
   with extended tag-junction homology.
3. **Constitutive plasmid expression** on pBBR1MCS2 — 2 primers, 2-fragment
   assembly.

The page itself lives in `src/views/PrimerDesignView.astro`. The Python
backend (Vercel serverless function) sits at `api/design/handler.py` and
delegates into `lib/primer_design/`. Genome FASTAs are streamed from a
Cloudflare R2 bucket — they are **not** committed.

The architectural contract (scope, hard-validation rules, application
shapes) lives in [`docs/primer_design/decisions_log.md`](docs/primer_design/decisions_log.md).

### Phase status

- **Phase 1 (current):** scaffolding only — folder tree, placeholder
  modules, navigation entry, deployment config. No algorithm yet.
- **Phase 2:** algorithm implementation, R2 genome upload, PDF writer,
  end-to-end wiring.

## Elastase Assay Calculator

A sub-page under `/elastase-assay` calculates volumes for a DQ Elastin
elastase assay from a single input (number of samples) and returns a PDF
report.

Since the tool has only one integer input, every possible report
(`num_samples` 1–1000) is pre-rendered offline and stored in the same R2
bucket used for the primer-design genomes, under the `elastase-assay/reports/`
prefix — regenerate/upload with `scripts/upload_elastase_reports.py`. The
runtime path (`api/elastase-assay.js`) is a small Node.js function that just
streams the matching object back; there's no Python involved at request
time, which avoids running two Python serverless functions in the same
deployment (`api/design` is the only one). The generation logic itself lives
in `lib/elastase_assay/` for local regeneration.
