# RMTIF 0.3.1

* Add published foundations for the restricted-time estimand and pairwise
  regression to DESCRIPTION with DOI links in CRAN format.
* Remove unsolicited progress printing and worker startup output from
  oncology studies. Results, diagnostics, and run times remain available in
  the returned objects; checkpoint writing remains available on request.
* Preserve the statistical calculations, random-seed streams, public
  function signatures, and equivalence of standalone and package sources.

# RMTIF 0.3.0

* Synchronize the package with the standalone version 1.2 functions.
* Provide estimators for complete ordered-state paths, terminal-aware IPCW,
  pairwise regression, pointwise curves, and Monte Carlo summaries.
* Add parameterized simulations and comparison methods for oncology endpoints.
* Document every public function and provide reproducible short examples.
