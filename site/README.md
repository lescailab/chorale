# CHORALE documentation site

The documentation published at <https://lescailab.github.io/chorale>, built with
[Astro](https://astro.build).

## Where the content comes from

Nothing on the site is written on the site. Two steps carry it from the sources:

| Step | Source | Becomes | Run by |
|---|---|---|---|
| `tools/manual_to_vignettes.R` | `knowledge/MANUAL.md` | `vignettes/{why,how-it-works,outputs,methods}.Rmd` | a person, when the manual changes |
| `site/sync.R` | `vignettes/*.Rmd`, `man/*.Rd` | the site pages, the navigation, the function reference | the workflow, on every push |

The manual is the source of record for the four narrative vignettes, and it holds the full account
that a methods section is written from. It sits outside the repository, so the vignettes it produces
are committed and continuous integration builds from those.

`vignettes/tutorial.Rmd` is written by hand and is never generated. It teaches the workflow rather
than describing the method, and has no counterpart in the manual.

Files under `site/src/pages/` are generated and are not committed, so editing one has no effect
beyond the next build.

## Building

```bash
Rscript tools/manual_to_vignettes.R   # only after editing the manual
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

Adding a narrative page means adding a part to the manual and a mapping in
`tools/manual_to_vignettes.R`. Adding a page with no counterpart in the manual, as the tutorial has
none, means adding a vignette by hand. `site_order` in a vignette's YAML sets where it appears in
the navigation, and `site_wide: true` widens the measure for a page carrying tables or equations.
