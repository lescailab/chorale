# Production-shaped benchmark

Profiling identified repeated multivariable factor-profile fitting as the dominant reusable
calculation in phenotype, bootstrap, and conditional permutation loops. `production.R` compares the
previous factor-by-factor regressions with the vectorised multi-response implementation on 2,000
samples, 20 factors, and mixed continuous and categorical covariates.

Run with:

```bash
Rscript benchmarks/production.R
```

The script fails when speedup is below 2x. End-to-end small-input regression remains covered by the
package test suite; the benchmark isolates the confirmed bottleneck so ICA variability and unrelated
I/O do not obscure the optimisation being tested.

Reference run with the portable R 4.5 environment: 0.202 seconds for the legacy loop and 0.068
seconds for the vectorised implementation, a 2.97x speedup over 30 repetitions. Timings vary by
machine; the executable 2x assertion is the acceptance criterion.
