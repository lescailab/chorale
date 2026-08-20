# CHORALE documentation site

The documentation published at <https://lescailab.github.io/chorale>, built with
[Astro](https://astro.build).

## Where the content comes from

Nothing on the generated narrative and reference pages is written on the site.

| Step | Source | Becomes | Run by |
|---|---|---|---|
| `site/sync.R` | `vignettes/*.Rmd`, `man/*.Rd` | the site pages and the function reference | before every build |

The repository vignettes are the canonical narrative sources. No build step reads files outside the
repository.

Files under `site/src/pages/` are generated and are not committed, so editing one has no effect
beyond the next build. The navigation is read from those pages rather than from a generated list,
so `npm run build` needs `site/sync.R` to have run but nothing else.

## Building

```bash
Rscript site/sync.R                   # from the package root, needs R
cd site
npm install                           # once
npm run dev                           # local server with live reload
npm run build                         # static output in site/dist
```

`.github/workflows/site.yaml` runs the same two steps on every push to `main` and publishes the
result.

## Notes

Mathematics is rendered at build time by KaTeX on the site, and to inline MathML by pandoc in the
vignettes, so neither runs a script in the reader's browser for it.

Adding a narrative page means adding a vignette. `site_order` in its YAML sets where it appears in
the navigation, and `site_wide: true` widens the measure for a page carrying tables or equations.
