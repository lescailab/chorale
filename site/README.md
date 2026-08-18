# CHORALE documentation site

The narrative documentation published at <https://lescailab.github.io/chorale>. Built with
[Astro](https://astro.build); the generated function reference is produced separately by pkgdown
and published beneath this site at `/reference/`.

```bash
npm install     # once
npm run dev     # local server with live reload
npm run build   # static output in dist/
npm run preview # serve the built output
```

`.github/workflows/site.yaml` builds and publishes on every push to `main`. Both this workflow and
the pkgdown one deploy to the `gh-pages` branch with `clean: false` and write to disjoint paths, so
neither removes the other's files.

Prose in `src/pages/why.md`, `how-it-works.md`, `outputs.md` and `methods.md` is derived from
`knowledge/MANUAL.md` in the project root, which is the source of record. Edit the manual and carry
the change across rather than letting the two drift. `tutorial.md` is written for the site and has
no counterpart in the manual; the code in it runs as shown.

Mathematics is rendered at build time by KaTeX, so no script runs in the reader's browser for it.
