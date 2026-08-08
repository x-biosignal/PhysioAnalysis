.inference_repo <- function() {
  candidates <- c(
    Sys.getenv("PHYSIO_REPO_ROOT"),
    file.path(testthat::test_path(), "..", "..", "..", ".."),
    getwd()
  )
  for (candidate in candidates[nzchar(candidates)]) {
    candidate <- normalizePath(candidate, mustWork = FALSE)
    if (file.exists(file.path(
      candidate, "physio-ecosystem", "validation", "inference", "surface.csv"
    ))) return(candidate)
  }
  NULL
}

.require_inference_repo <- function() {
  repo <- .inference_repo()
  if (is.null(repo)) {
    skip("central inference fixtures are available only in the monorepo")
  }
  repo
}

test_that("WS8 PhysioAnalysis exports pass the offline parity gate", {
  repo <- .require_inference_repo()
  script <- file.path(
    repo, "physio-ecosystem", "validation", "inference", "run_parity.R"
  )
  output <- tempfile("physioanalysis-parity-")
  log <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      shQuote(script), "--packages", "PhysioAnalysis",
      "--output", shQuote(output), "--fail-fast"
    ),
    stdout = TRUE, stderr = TRUE
  )
  expect_null(attr(log, "status"), info = paste(log, collapse = "\n"))
  results <- utils::read.csv(
    file.path(output, "inference_parity.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(nrow(results) >= 25L)
  expect_true(all(results$status == "PASS"))
})

test_that("asymmetric SPM fixtures detect a reversed contrast", {
  repo <- .require_inference_repo()
  root <- file.path(
    repo, "physio-ecosystem", "validation", "inference"
  )
  input <- readRDS(file.path(
    root, "fixtures", "spm", "spm-v1", "input.rds"
  ))
  forward <- spmTTest(input$two_sample, input$group1, input$group2)
  reverse <- spmTTest(input$two_sample, input$group2, input$group1)
  expect_equal(forward$t, -reverse$t, tolerance = 1e-12)
  expect_gt(abs(forward$t[5]), 0.1)
})
