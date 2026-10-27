# Compatibility helpers for the pre-Calibre EpiLink API.

# Results created by older code store their payload under @misc$EpiLink. The
# canonical Calibre API stores the same payload under @misc$Calibre, while
# mirroring the old key so old callers can continue to inspect returned
# Seurat objects during the migration period.
.calibre_sync_object <- function(object) {
  if (!inherits(object, "Seurat")) return(object)
  if (is.null(object@misc$Calibre) && !is.null(object@misc$EpiLink)) {
    object@misc$Calibre <- object@misc$EpiLink
  }
  if (!is.null(object@misc$Calibre)) {
    object@misc$EpiLink <- object@misc$Calibre
  }
  object
}

# Historical method IDs remain epilink_* in old result tables. New code uses
# calibre_* IDs and can normalize old IDs at public entry points.
calibre_output_method_id <- function(method_id) {
  sub("^epilink_", "calibre_", as.character(method_id))
}

calibre_legacy_method_id <- function(method_id) {
  sub("^calibre_", "epilink_", as.character(method_id))
}

calibre_normalize_method_ids <- function(method_ids) {
  calibre_output_method_id(method_ids)
}

calibre_result_directory <- function(root, method_id) {
  ids <- unique(c(
    as.character(method_id),
    calibre_output_method_id(method_id),
    calibre_legacy_method_id(method_id)
  ))
  directories <- file.path(root, ids)
  existing <- directories[dir.exists(directories)]
  if (length(existing)) existing[[1L]] else directories[[1L]]
}

calibre_object_path <- function(root, method_id) {
  directory <- calibre_result_directory(root, method_id)
  candidates <- file.path(directory, c("calibre_object.rds", "epilink_object.rds"))
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1L]] else candidates[[1L]]
}

# Public API aliases retained for existing scripts and serialized workflows.
EpiLink_single <- function(...) Calibre_single(...)
EpiLink <- function(...) .calibre_sync_object(Calibre(...))
EpiLinkDynamic <- function(...) .calibre_sync_object(CalibreDynamic(...))
EpiLinkDiff <- function(...) .calibre_sync_object(CalibreDiff(...))
EpiLinkFilter <- function(object, ...) {
  CalibreFilter(.calibre_sync_object(object), ...)
}

epilink_pair_seeds <- function(...) calibre_pair_seeds(...)
with_epilink_seed <- function(seed, code) with_calibre_seed(seed, code)
simulate_epilink_data <- function(...) simulate_calibre_data(...)
epilink_simulation_methods <- function(...) calibre_simulation_methods(...)
run_epilink_simulation <- function(...) run_calibre_simulation(...)
simulate_epilink_dynamic_data <- function(...) simulate_calibre_dynamic_data(...)
fit_empirical_epilink_family <- function(...) fit_empirical_calibre_family(...)
fit_da_epilink <- function(...) fit_da_calibre(...)
