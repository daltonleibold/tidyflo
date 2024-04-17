# func.R

#-------------------
# statistics
#-------------------
#----
#' @title geo.mean
#' @param x a vector
#' @description calculates the geometric mean a robust metric of the central tendency of a population. applicable to values on a log-normal scale is robust to outliers will only ever be at most equal to the arithmetic mean
#' @export 
geo.mean <- function(x, na.rm = T){
    10^base::mean(log.sc(x), na.rm = na.rm)
}
#----
#' @title se
#' @param x a numeric vector
#' @description calculates the standard error
#' @export 
se <- function(x, na.rm = T){
  stats::sd(x, na.rm = na.rm) / base::sqrt(base::length(x))
}
#----
#' @title r.sd
#' @param x a numeric vector
#' @description calculates the robust standard deviation; based on median absolute deviation. 
#' @export
r.sd <- function(x, na.rm = T) {
  1.4826 * stats::mad(x, constant = 1.4826, na.rm = na.rm)
}
#----
#' @title r.se
#' @param x a numeric vector
#' @description calculates the robust standard error
#' @export
r.se <- function(x, na.rm = T) {
  r.se = r.sd(.val, na.rm = na.rm) / base::sqrt(.n)
}
#----
#' @title log.sc
#' @param x a numeric vector
#' @description scales a vector to log10 allows for negative or zero values
#' @export
log.sc <- function(x) {
  base::suppressWarnings(base::ifelse(x <= 0, -base::log10(base::abs(x) + 1), base::log10(x + 1)))
}

#----
#' @title summary.fcs
#' @param data a dataframe of fcs data
#' @param na.rm a logical of whether or not to remove missing values from calculations. defaults to TRUE.
#' @param batch a logical of whether to batch process or not. if TRUE, groups by every column except value (.val) and cell number (.n)
#' @description a tidyverse wrapper for calculating common summary statistics for flow cytometry data. includes number of cells (n), median (median), geometric mean (geo.mean), arithmetic mean (mean), standard deviation (sd), standard error (se), robust standard deviation (r.sd), and robust standard error (r.se).
#' 
summary.fcs <- function(data, na.rm = TRUE){
  data %>%
    dplyr::reframe(.n = dplyr::n(),
            median = stats::median(.val, na.rm = na.rm),
            geo.mean = geo.mean(.val, na.rm = na.rm),
            mean = base::mean(.val, na.rm = na.rm),
            sd = stats::sd(.val, na.rm = na.rm),
            se = se(.val, na.rm = na.rm),
            r.sd = r.sd(.val, na.rm = na.rm),
            r.se = r.sd(.val, na.rm = na.rm) / base::sqrt(.n))
}

#----
#' @title re.code
#' @param x a character or factor vector
#' @param mapping a vector of matched channel names to label names
#' @description function for re-labelling channels with new parameter names
#' @export
re.code <- function(x, mapping) {
  labels <- base::unname(mapping)
  base::names(labels) <- base::names(mapping)
  out <- labels[base::match(x, base::names(mapping))]
  out[base::is.na(out)] <- x[base::is.na(out)]
  return(out)
}

