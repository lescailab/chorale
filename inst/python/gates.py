"""The two distributional conditions of applicability.

Both are evaluated on independent components, and both are calibrated rather
than tested in absolute terms.

The non-Gaussianity condition asks whether a modality's recovered components
are further from normal than components recovered the same way from data that
is Gaussian by construction. Independent component analysis maximises
non-Gaussianity, so an unconditional normality test on its output rejects
almost always and carries no evidence. Each modality is therefore paired with
surrogates of identical dimension and covariance, pushed through the identical
pipeline, and the condition is whether the real components stand outside that
distribution.

The modality-difference condition asks whether the component distributions of
two modalities differ from one another. Identification rests on the modalities
not being copies of one another in distribution, so agreement here is a
failure, not a success.

Components are recovered by a callable supplied by the caller, so the estimator
that decides these conditions is the estimator that will consume them. Nothing
in this file fits a component itself.
"""

import numpy as np
from scipy import stats


def anderson_critical_value(result, alpha):
  """The Anderson-Darling critical value at a level, found by its level.

  `scipy` returns critical values alongside the significance levels they
  belong to. Reading one out by position assumes an ordering the library does
  not promise, so the level is matched instead and a level it does not carry
  is an error rather than a neighbouring value silently used.
  """
  levels = np.asarray(result.significance_level, dtype=float)
  wanted = float(alpha) * 100.0
  hit = np.flatnonzero(np.isclose(levels, wanted))
  if hit.size == 0:
    raise ValueError(
      f"the Anderson-Darling table carries levels {sorted(levels)} per cent, "
      f"not {wanted}; choose an alpha the table carries")
  return float(np.asarray(result.critical_values, dtype=float)[hit[0]])


def component_stats(sources, alpha=0.05):
  """Per-component distributional summaries of a samples-by-components array.

  Returns plain Python values rather than a data frame: this is the boundary
  the caller reads across, and a scalar that arrives as a foreign array is a
  scalar the caller cannot compare.
  """
  rows = []
  for j in range(sources.shape[1]):
    v = np.asarray(sources[:, j], dtype=float)
    v = v[np.isfinite(v)]
    ad = stats.anderson(v, dist="norm")
    rows.append({
      # One-based, because the caller counts factors from one.
      "component": j + 1,
      "anderson_A2": float(ad.statistic),
      "anderson_critical_value": anderson_critical_value(ad, alpha),
      "excess_kurtosis": float(stats.kurtosis(v, fisher=True)),
      "skew": float(stats.skew(v)),
      "skewtest_p": float(stats.skewtest(v).pvalue) if len(v) >= 8 else float("nan"),
    })
  return rows


def _as_sources(ica_fn, x, k, seed):
  """Call the caller's estimator and return standardised sources.

  The matrix crosses a language boundary in both directions, so its
  orientation is checked rather than assumed. A transposed return would be
  silently accepted as a different set of components, and every verdict below
  would be computed on it.
  """
  x = np.asarray(x, dtype=float)
  if x.ndim != 2:
    raise ValueError(f"expected a samples-by-features matrix, got {x.shape}")
  s = np.asarray(ica_fn(x, int(k), int(seed)), dtype=float)
  if s.ndim != 2 or s.shape != (x.shape[0], int(k)):
    raise ValueError(
      f"the component estimator returned {s.shape}, expected "
      f"({x.shape[0]}, {k}) samples by components")
  sd = s.std(axis=0)
  sd[sd == 0] = 1.0
  return (s - s.mean(axis=0)) / sd


def gate_nongaussianity(x, k, label, ica_fn, seed=1, n_surrogate=100,
                        alpha=0.05):
  """Is the modality further from normal than a Gaussian of its covariance?

  The summary statistic is the median Anderson-Darling statistic over a
  modality's components, compared with the distribution of that statistic
  across surrogates. Ranking a modality's components against one surrogate's
  instead would tie the resolution of the test to how many components the
  modality admits, and a modality admitting two could not be assessed at any
  conventional level whatever its distributions.
  """
  rng = np.random.default_rng(seed)
  real = component_stats(_as_sources(ica_fn, x, k, seed), alpha)
  for row in real:
    row["source"] = "observed"
  observed = float(np.median([row["anderson_A2"] for row in real]))

  n, p = x.shape
  cov = np.cov(x, rowvar=False)
  evals, evecs = np.linalg.eigh(cov)
  evals = np.clip(evals, 0, None)
  root = (evecs * np.sqrt(evals)) @ evecs.T

  null = np.empty(n_surrogate)
  detail = list(real)
  for b in range(n_surrogate):
    surrogate = rng.standard_normal((n, p)) @ root
    stats_b = component_stats(_as_sources(ica_fn, surrogate, k, seed), alpha)
    null[b] = float(np.median([row["anderson_A2"] for row in stats_b]))
    if b == 0:
      for row in stats_b:
        row["source"] = "gaussian_surrogate"
      detail.extend(stats_b)

  p_value = float((1 + np.sum(null >= observed)) / (1 + n_surrogate))
  for row in detail:
    row["modality"] = label

  summary = {
    "modality": label,
    "k_components": int(k),
    "n_surrogate": int(n_surrogate),
    "median_A2_observed": round(observed, 4),
    "median_A2_surrogate": round(float(np.median(null)), 4),
    "p_value": p_value,
    "pct_components_rejecting_normality": round(100.0 * float(np.mean(
      [row["anderson_A2"] > row["anderson_critical_value"] for row in real])), 1),
    "pct_components_asymmetric": round(100.0 * float(np.mean(
      [row["skewtest_p"] < alpha for row in real])), 1),
    "verdict": "pass" if p_value < alpha else "fail",
  }
  return summary, detail


def gate_modality_difference(sources, alpha=0.05):
  """Do the component distributions of two modalities differ?

  Each pair is compared by pooling its components' values and by counting the
  component pairs a two-sample test cannot separate. A pair passes only when
  the pooled comparison separates them in both directions of the pairing, so a
  single extreme component cannot carry the verdict.
  """
  names = list(sources)
  rows = []
  for i in range(len(names)):
    for j in range(i + 1, len(names)):
      a, b = np.asarray(sources[names[i]]), np.asarray(sources[names[j]])
      pooled = stats.ks_2samp(a.ravel(), b.ravel())
      pooled_neg = stats.ks_2samp(a.ravel(), -b.ravel())
      indistinguishable = 0
      total = 0
      for ca in range(a.shape[1]):
        for cb in range(b.shape[1]):
          total += 1
          if stats.ks_2samp(a[:, ca], b[:, cb]).pvalue >= alpha:
            indistinguishable += 1
      rows.append({
        "pair": f"{names[i]} vs {names[j]}",
        "pooled_KS_D": round(float(pooled.statistic), 4),
        "pooled_KS_p": float(pooled.pvalue),
        "pooled_KS_p_sign_flipped": float(pooled_neg.pvalue),
        "n_component_pairs": int(total),
        "pct_pairs_indistinguishable": round(
          100.0 * indistinguishable / total, 1) if total else float("nan"),
        "verdict": ("pass" if (pooled.pvalue < alpha and
                               pooled_neg.pvalue < alpha) else "fail"),
      })
  return rows
