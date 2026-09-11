## Submission summary

This is a new submission of PhysioAnalysis, a package providing analysis and
visualization functions for physiological signal data stored in
'PhysioExperiment' objects (built on Bioconductor's 'SummarizedExperiment').

## Test environments

* Local: Ubuntu 24.04, R 4.5.2 (release)
* win-builder: R-devel and R-release (planned)
* macbuilder: R-release (planned)
* R-hub v2 (GitHub Actions): linux, windows, macos (planned)

## R CMD check results

`R CMD check --as-cran` was run locally (R 4.5.2). There are 0 ERRORs and no
package-level WARNINGs. The check reports the standard "New submission" NOTE
under CRAN incoming feasibility.

The compiled-code check ("checking compiled code") passes cleanly.

### Notes / expected messages

* checking CRAN incoming feasibility ... NOTE
  New submission.

  This is a legitimate first submission to CRAN. On CRAN's own incoming
  pipeline this may additionally surface "unable to verify current time" and
  "possibly misspelled words in DESCRIPTION" (domain terms such as FFT,
  topoplots, multichannel, Bioconductor); these are expected for a new
  submission and are false positives / standard domain vocabulary.

* "Strong dependencies not in mainstream repositories: PhysioCore"

  PhysioCore is the ecosystem's base package and is submitted to CRAN before
  PhysioAnalysis (see "Reverse dependencies" below). This message resolves once
  PhysioCore is available on CRAN.

* The local check also emits "'qpdf' is needed for checks on size reduction of
  PDFs". This is purely a limitation of the local machine (qpdf is not
  installed); the vignettes are HTML and the package ships no PDFs. qpdf is
  available on CRAN, win-builder and R-hub, where this message does not appear.

## Reverse dependencies

This package has no reverse dependencies on CRAN yet.

PhysioAnalysis is part of the PhysioExperiment ecosystem. The intra-ecosystem
sibling packages are submitted in dependency order:

    PhysioCore -> PhysioIO / PhysioPreprocess -> PhysioAnalysis

PhysioAnalysis depends on PhysioCore and is submitted after PhysioCore is
available on CRAN.
