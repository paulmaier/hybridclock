#' hybridclock: How long have hybrids been forming?
#'
#' Estimates the history of hybrid populations from diagnostic markers that
#' need no reference genome: the average age of the hybrids, how long ago the
#' rarer species' ancestry entered, how many times the hybrids have
#' backcrossed, and whether new hybridization is still happening. Population
#' samples are compared with simulated histories that include drift, small
#' breeding groups, sibling mating and separate breeding pools, so that
#' inbreeding is not mistaken for age. Recent hybrids can also be
#' classified individually.
#'
#' The usual workflow is [diagnostic_markers()] and [hybrid_scores()], then
#' [classify_hybrids()] for individuals, and [build_reference()],
#' [train_clock()] and [estimate_history()] for populations.
#'
#' @keywords internal
#' @importFrom Cubist cubist
#' @importFrom ranger ranger
"_PACKAGE"
