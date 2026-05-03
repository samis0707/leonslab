import { defineCollection, z } from 'astro:content';
import { glob } from 'astro/loaders';

const cartoons = defineCollection({
  loader: glob({ pattern: '*.md', base: './src/content/cartoons' }),
  schema: z.object({
    order: z.number(),
    title: z.string(),
    tag: z.enum(['lab-life', 'slice-of-life', 'pop-parody', 'pseudomonas-park']),
    image: z.string(),
    alt: z.string(),
    aspect: z.string().default('1/1'),
    year: z.number().optional(),
    captionEn: z.string().optional(),
    captionDe: z.string().optional(),
  }),
});

const stickers = defineCollection({
  loader: glob({ pattern: '*.md', base: './src/content/stickers' }),
  schema: z.object({
    order: z.number(),
    title: z.string(),
    year: z.number(),
    size: z.string(),
    material: z.string(),
    image: z.string(),
    imageContext: z.string().optional(),
    alt: z.string(),
    altContext: z.string().optional(),
    aspect: z.string().default('1/1'),
    storyEn: z.string(),
    storyDe: z.string(),
  }),
});

export const collections = { cartoons, stickers };
