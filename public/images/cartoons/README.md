# Cartoons — image upload guide

This folder is where the **18 cartoon images** belong. Until they are
uploaded here, the gallery on the site shows clean placeholder tiles with
just the title — once a file is dropped in with the correct name, it appears
automatically on the next build.

## Expected filenames

Drop each cartoon as a `.png` (preferred for line art with transparent
background) or `.jpg` with this exact name:

| #  | Filename                              | Title                          |
|----|---------------------------------------|--------------------------------|
| 01 | `01-under-observation.png`            | Under Observation              |
| 02 | `02-bacillus-family.png`              | The Bacillus Family            |
| 03 | `03-duplication-rate.png`             | Duplication Rate               |
| 04 | `04-beach-day.png`                    | Beach Day                      |
| 05 | `05-cultural-differences.png`         | Cultural Differences           |
| 06 | `06-miniprep-fun.png`                 | Who said Miniprep is fun?      |
| 07 | `07-stop-bacteria-testing.png`        | Stop Bacteria Testing          |
| 08 | `08-sterilizing-experience.png`       | A Sterilizing Experience       |
| 09 | `09-resisdance.png`                   | ResisDANCE                     |
| 10 | `10-amazing-arachnobacterium.png`     | The Amazing Arachnobacterium   |
| 11 | `11-thoughts-of-a-scientist.png`      | Thoughts of a Scientist        |
| 12 | `12-pure-evil.png`                    | The Pure Evil                  |
| 13 | `13-parental-advice.png`              | Parental Advice                |
| 14 | `14-adaptive-response.png`            | Adaptive Response              |
| 15 | `15-tool-for-every-situation.png`     | A Tool for Every Situation     |
| 16 | `16-group-photo.png`                  | Group Photo                    |
| 17 | `17-big-culture-life.png`             | Big-Culture-Life               |
| 18 | `18-pore-opening-experience.png`      | A Pore-Opening Experience      |

## Recommendations

- **Format**: PNG with transparent background, otherwise JPG.
- **Resolution**: long edge between 2000–3000 px is ideal.
- **Aspect ratio**: keep the original — the gallery is masonry-style and
  every cartoon's declared `aspect` in `src/content/cartoons/*.md` reserves
  the correct space.

## Uploading via GitHub web UI

1. Open the repo on GitHub on the branch `claude/create-portfolio-site-1nPvv`.
2. Navigate to `public/images/cartoons/`.
3. Click **Add file → Upload files**.
4. Drag the 18 images in (named as above), commit on the branch.

## If you need to use `.jpg` instead of `.png`

Edit `src/content/cartoons/<slug>.md` and change the `image:` field's
extension. The metadata is the only source of truth.
