import { S3Client, GetObjectCommand } from '@aws-sdk/client-s3';

// The elastase assay calculator has a single integer input, so every possible
// report (num_samples 1..MAX_SAMPLES) is pre-rendered and stored in R2 by
// scripts/upload_elastase_reports.py. This function just streams the right
// object back — no Python involved at request time.
const MAX_SAMPLES = 100;
const R2_PREFIX = 'elastase-assay/reports/';

let client;
function r2Client() {
  if (!client) {
    client = new S3Client({
      region: 'auto',
      endpoint: `https://${process.env.R2_ACCOUNT_ID}.r2.cloudflarestorage.com`,
      credentials: {
        accessKeyId: process.env.R2_ACCESS_KEY_ID,
        secretAccessKey: process.env.R2_SECRET_ACCESS_KEY,
      },
    });
  }
  return client;
}

export default async function handler(req, res) {
  res.setHeader('Access-Control-Allow-Origin', '*');

  if (req.method === 'OPTIONS') {
    res.setHeader('Access-Control-Allow-Methods', 'GET, OPTIONS');
    res.setHeader('Access-Control-Allow-Headers', 'Content-Type');
    res.status(204).end();
    return;
  }

  const raw = req.query?.num_samples ?? new URL(req.url, 'http://localhost').searchParams.get('num_samples');
  const numSamples = Number(raw);

  if (!Number.isInteger(numSamples) || numSamples <= 0 || numSamples > MAX_SAMPLES) {
    res.status(400).json({
      status: 'error',
      message: `num_samples must be an integer between 1 and ${MAX_SAMPLES}`,
    });
    return;
  }

  const key = `${R2_PREFIX}${numSamples}.pdf`;

  try {
    const obj = await r2Client().send(
      new GetObjectCommand({ Bucket: process.env.R2_BUCKET_NAME, Key: key }),
    );
    const bytes = await obj.Body.transformToByteArray();

    res.setHeader('Content-Type', 'application/pdf');
    res.setHeader(
      'Content-Disposition',
      `attachment; filename="Elastase_Assay_Report_${numSamples}_samples.pdf"`,
    );
    res.setHeader('Cache-Control', 'public, max-age=31536000, immutable');
    res.status(200).send(Buffer.from(bytes));
  } catch (err) {
    if (err?.name === 'NoSuchKey') {
      res.status(404).json({ status: 'error', message: 'Report not found' });
      return;
    }
    console.error('elastase-assay: R2 fetch failed', err);
    res.status(500).json({ status: 'error', message: 'Internal error' });
  }
}
