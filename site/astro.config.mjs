import { defineConfig } from 'astro/config';
import { unified } from '@astrojs/markdown-remark';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

export default defineConfig({
  site: 'https://lescailab.github.io',
  base: '/chorale',
  outDir: './dist',
  trailingSlash: 'ignore',
  markdown: {
    processor: unified({
      remarkPlugins: [remarkMath],
      rehypePlugins: [rehypeKatex],
    }),
    shikiConfig: { // The high-contrast pair. Every token colour either of them assigns
    // clears WCAG AA against this site's two page grounds; the plain
    // github-light and github-dark do not.
    themes: { light: 'github-light-high-contrast', dark: 'github-dark-high-contrast' }, wrap: true },
  },
});
