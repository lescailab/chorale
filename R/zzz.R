.onLoad <- function(libname, pkgname) {
  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  if (nzchar(conda_prefix)) {
    reticulate::use_condaenv(conda_prefix, required = FALSE)
  }
}

#' Set up the Python environment for chorale
#'
#' Where chorale is loaded inside an already-active conda environment, that
#' environment's interpreter is used directly (see `R/zzz.R` `.onLoad`) and
#' this function is not needed. Otherwise, it creates a managed conda
#' environment from the package's `inst/python/environment.yml` the first
#' time it is called, and binds reticulate to it. Delegating environment
#' creation to basilisk is deliberately excluded, since its isolated
#' environment would duplicate and shadow the interpreter a conda
#' installation already provides.
#'
#' @param envname Name of the managed conda environment to create or reuse.
#' @returns Invisibly, the path to the conda environment's Python
#'   interpreter.
#' @export
#' @examplesIf interactive()
#' chorale_python_setup()
chorale_python_setup <- function(envname = "r-chorale") {
  rlang::check_installed("reticulate")

  yml <- system.file("python", "environment.yml", package = "chorale")
  if (!nzchar(yml)) {
    rlang::abort("Could not locate inst/python/environment.yml in the installed package.")
  }

  existing <- tryCatch(reticulate::conda_list()$name, error = function(e) character())
  if (!envname %in% existing) {
    reticulate::conda_create(envname, environment = yml)
  }
  reticulate::use_condaenv(envname, required = TRUE)
  invisible(reticulate::conda_python(envname))
}
