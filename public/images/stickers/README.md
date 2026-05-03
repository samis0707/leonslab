# Stickers — image upload guide

This folder is where the sticker images belong. Until they are uploaded
the page renders clean placeholder tiles with the sticker's name; as soon
as a file lands here with the correct name, it appears on the next build.

## Expected files

| File                                            | What it is                                                       |
|-------------------------------------------------|------------------------------------------------------------------|
| `hero-plants.jpg`                               | Page hero — the photo of both stickers together on plants        |
| `01-pseudomonas-park-cutout.png`                | Pseudomonas Park sticker, transparent background (cutout)        |
| `01-pseudomonas-park-context.jpg`               | Pseudomonas Park sticker on the raw exposed concrete wall        |
| `02-pipetboy-cutout.png`                        | Pipetboy sticker, transparent background (cutout)                |
| `02-pipetboy-context.jpg`                       | Pipetboy sticker held in hand, building facade in the background |

## Recommendations

- **Cutouts**: PNG with transparent background. Long edge ~2000 px is plenty.
- **Context shots**: JPG, long edge 2400–3000 px. The page hero is shown at
  a 16:9 aspect — pick a crop that works horizontally.
- **Color**: keep the originals — the page background is concrete-tinted and
  picks up the sticker colors well.

## Uploading via GitHub web UI

1. Open the repo on GitHub on the branch `claude/create-portfolio-site-1nPvv`.
2. Navigate to `public/images/stickers/`.
3. Click **Add file → Upload files** and drop the images in (named exactly
   as above).
4. Commit on the branch.

## Adding a third (or more) sticker

1. Create `src/content/stickers/03-<slug>.md` with the same frontmatter
   shape as the existing ones (`order`, `title`, `year`, `size`, `material`,
   `image`, `imageContext`, `alt`, `altContext`, `aspect`, `storyEn`,
   `storyDe`).
2. Drop `03-<slug>-cutout.png` (and optional `-context.jpg`) into this
   folder.
3. The page picks the new sticker up automatically and alternates the
   layout direction.
