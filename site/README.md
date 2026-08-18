# CHORALE documentation site

The documentation published at <https://lescailab.github.io/chorale>, built with
[Astro](https://astro.build).

## Everything on the site comes from the package

`site/sync.R` generates the site's content from the package, so the two cannot diverge:

| Source | Becomes |
|---|---|
| `vignettes/*.Rmd` | the narrative pages, and the navigation order |
| `man/*.Rd` | the function reference and its index |

The vignettes are the source of record. Files under `site/src/pages/` are generated and are not
committed, so editing one has no effect beyond the next build. Edit the vignette instead. The same
files ship inside the package, so a reader offline sees what a reader on the site sees.

## Building

```bash
Rscript site/sync.R   # from the package root, needs R
cd site
npm install           # once
npm run dev           # local server with live reload
npm run build         # static output in site/dist
```

`.github/workflows/site.yaml` runs the same two steps on every push to `main` and publishes the
result.

## Notes

Mathematics is rendered at build time by KaTeX on the site, and to inline MathML by pandoc in the
vignettes, so neither runs a script in the reader's browser for it.

Adding a page means adding a vignette. `site_order` in its YAML sets where it appears in the
navigation, and `site_wide: true` widens the measure for a page carrying tables or equations.
