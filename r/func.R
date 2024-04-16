# func.R

#-------------------
# statistics
#-------------------
#----
#' @name geo.mean
#' @param x a vector
#' @description
#' calculates the geometric mean
#' a robust metric of the central tendency of a population. 
#' applicable to values on a log-normal scale
#' is robust to outliers
#' will only ever be at most equal to the arithmetic mean
geo.mean <- function(x, na.rm = T){
    10^base::mean(log.sc(x), na.rm = na.rm)
}
#----
#' @name se
#' @param x a numeric vector
#' @description
#' calculates the standard error
se <- function(x, na.rm = T){
  stats::sd(x, na.rm = na.rm) / base::sqrt(base::length(x))
}
#----
#' @name r.sd
#' @param x a numeric vector
#' @description 
#' calculates the robust standard deviation; based on median absolute deviation. 
r.sd <- function(x, na.rm = T) {
  1.4826 * stats::mad(x, constant = 1.4826, na.rm = na.rm)
}
#----
#' @name r.se
#' @param x a numeric vector
#' @description 
#' calculates the robust standard error
r.se <- function(x, na.rm = T) {
  r.se = r.sd(.val, na.rm = na.rm) / base::sqrt(.n)
}
#----
#' @name log.sc
#' @param x a numeric vector
#' @description 
#' scales a vector to log10
#' allows for negative or zero values
log.sc <- function(x) {
  base::suppressWarnings(base::ifelse(x <= 0, -base::log10(base::abs(x) + 1), base::log10(x + 1)))
}
#----
#' @name summary.fcs
#' @param data a dataframe of fcs data
#' @param na.rm a logical of whether or not to remove missing values from calculations. defaults to TRUE.
#' @param batch a logical of whether to batch process or not. if TRUE, groups by every column except value (.val) and cell number (.n)
#' @description 
#' a tidyverse wrapper for calculating common summary statistics for flow cytometry data. includes number of cells (n), median (median), geometric mean (geo.mean), arithmetic mean (mean), standard deviation (sd), standard error (se), robust standard deviation (r.sd), and robust standard error (r.se).
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

#-------------------
# data processing
#-------------------
#----
#' @name read.fcs
#' @param file a path to a *.fcs file
#' @description 
#' simple wrapper that reads a *.fcs to a dataframe
read.fcs <- function(path){
  data <- base::list()
  if(stringr::str_detect(path, ".fcs$") == T){
    file <- path
    name <- file
  } else {
    file <- base::list.files(path, ".fcs$", full.names = T)
    name <- base::list.files(path, ".fcs$", full.names = F)
  }
  for(i in 1:base::length(file)){
    df <- base::data.frame(flowCore::read.FCS(file[i], transformation = F)@exprs)
    data[[i]] <- df
    base::names(data)[[i]] <- stringr::str_remove(name[i], ".fcs$")
  }
  data <- plyr::ldply(data)
  return(data)
}
#----
#' @name tidy.fcs
#' @param data a dataframe read in from .fcs
#' @description 
#' a tidyr pipeline for standardizing .fcs data
tidy.fcs <- function(data){
  data %>%
    dplyr::rename_with(.cols = tidyselect::everything(),
                .fn = tidyselect::str_to_lower) %>%
    dplyr::select(-time) %>%
    dplyr::group_by(.id) %>%
    dplyr::mutate(.n = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(-c(.n, .id),
                 names_to = ".ch",
                 values_to = ".val") %>%
    dplyr::mutate(.param = stringr::str_extract(.ch, ".[a|h|w]$"),
           .ch = stringr::str_remove(.ch, ".[a|h|w]$")) %>%
    dplyr::select(.id, .n, .ch, .param, .val) %>%
    dplyr::arrange(.id, .n, .ch, .param, .val)
}
#----
#' @name process.fcs
#' @param data a dataframe
#' @description 
#' a pipeline for setting a standard processing protocol on flow cytometry data. log transforms all values, filters offscale events, then gates to most "average" cells in the population - within 1 s.d. of the median in all channels.
process.fcs <- function(data){
  data %>% 
    dplyr::mutate(.val = log.sc(.val)) %>%
    tidyr::pivot_wider(names_from = ".ch",
                values_from = ".val") %>%
    dplyr::filter(dplyr::if_all(.cols = -c(.id, .n, .param),
                  .fns = ~.x >= 0 & .x <= 5)) %>%
    dplyr::group_by(.id) %>%
    dplyr::filter(dplyr::if_all(.cols = -c(.n, .param),
                  .fns = ~.x <= stats::median(.x) + (stats::median(.x) * r.sd(.x)) &
                    .x >= stats::median(.x) - (stats::median(.x) * r.sd(.x)))) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = -c(.id, .n, .param),
                 names_to = ".ch",
                 values_to = ".val")
}
#----
#' @name re.code
#' @param x a character or factor vector
#' @param mapping a vector of matched channel names to label names
#' @description 
#' function for re-labelling channels with new parameter names
re.code <- function(x, mapping) {
  labels <- base::unname(mapping)
  base::names(labels) <- base::names(mapping)
  out <- labels[base::match(x, base::names(mapping))]
  out[base::is.na(out)] <- x[base::is.na(out)]
  return(out)
}
#----
#' @name extract.scc
#' @param model a brms object
#' @param neg a character
#' @description 
#' function for extracting posterior distributions generated by a simple intercept model
#' tailored for single color controls
extract.scc <- function(model, neg = "af"){
  post <- model %>% 
    brms::as_draws_df() %>%
    tibble::as_tibble() %>%
    dplyr::select(base::grep("^b_", base::colnames(.))) %>%
    dplyr::rename_with(.cols = tidyselect::everything(),
                .fn = ~stringr::str_remove(.x, "^b_"))
  # true negative
  # autofluorescence channel in the unstained sample
  a. <- post %>%
    reframe(.ch = neg,
            .samp = neg,
            .val = Intercept)
  # autofluorescence spillover / channel background
  # fluorescent signal in each channel of the unstained control
  b. <- post %>%
    select(Intercept, grep("^.ch", colnames(.))) %>%
    select(-grep(":", colnames(.))) %>%
    rename_with(.cols = everything(),
                .fn = ~str_remove(.x, "^.ch")) %>%
    pivot_longer(-Intercept,
                 names_to = ".ch",
                 values_to = ".val") %>%
    mutate(.val = Intercept + .val,
           .samp = "af") %>%
    select(.ch, .samp, .val)
  # dye-induce autofluorescence
  # autofluorescent channel in each single color control
  c. <- post %>%
    select(Intercept, grep("^.samp", colnames(.))) %>%
    select(-grep(":", colnames(.))) %>%
    rename_with(.cols = everything(),
                .fn = ~str_remove(.x, "^.samp")) %>%
    pivot_longer(-Intercept,
                 names_to = ".samp",
                 values_to = ".val") %>%
    mutate(.val = Intercept + .val,
           .ch = neg) %>%
    select(.ch, .samp, .val)
  # fluorescent spillover
  # each of the fluorescent channels for each of the single color controls
  d. <- post %>%
    select(Intercept, grep(":", colnames(.))) %>%
    rename_with(.cols = everything(),
                .fn = ~str_remove(.x, "^.samp")) %>%
    rename_with(.cols = everything(),
                .fn = ~str_remove(.x, ".ch")) %>%
    pivot_longer(-Intercept,
                 names_to = ".samp.ch",
                 values_to = ".val") %>%
    separate(col = .samp.ch,
             into = c(".samp", ".ch"),
             sep = ":") %>%
    select(.ch, .samp, .val)
  # bind the extracted distributions back together
  data <- rbind(a., b., c., d.)
  return(data)
}
#-------------------
# graphing
#-------------------
#----
#' @title theme_custom
#' @description a custom theme for ggplots that makes publication-ready figures
theme_custom <- theme_classic(base_size = 24) +
  # adjust axis title position
  theme(axis.title.y=element_text(vjust=1.5), 
        axis.title.x=element_text(vjust=0.2)) + 
  # adjust plot margins and line element size
  theme(plot.margin = unit(c(.3,.3,.6,.6), "cm"), 
        line = element_line(linewidth = 1.25)) + 
  # draw x and y axes
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")) + 
  #put margins around axis labels so that nothing overlaps
  theme(axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) + 
  # move tickmarks inside the axes and paint black
  theme(axis.ticks.length =unit(-0.3, "cm")) + 
  #spread out facets
  theme(panel.spacing = unit(2, units = "lines")) + 
  #make tick marks black
  theme(axis.ticks = element_line(color = "black")) + 
  #remove border from facet labels
  theme(strip.background = element_blank()) 
#----
#' @name plot.wave
#' @param data a dataframe of flow cytometry data
#' #' @param batch a logical of whether to batch process or not. if TRUE, adds additional facets by sample id. defaults to FALSE.
#' @description
#' wrapper for making a ggplot object that plots the fluorescent intensities in each channel of a dataframe as density plots (waves). 
plot.wave <- function(data){
  data %>%
    ggplot2::ggplot(aes(x = .val, fill = .ch, alpha = 0.1)) +
    ggplot2::geom_density() +
    ggplot2::facet_grid(rows = vars(.ch)) +
    ggplot2::guides(alpha = "none", fill = "none") +
    ggplot2::theme_classic()
}
#----
#' @name plot.density
#' @param data a dataframe of flow cytometry data
#' @param x an object; the channel to plot on the x axis; defaults to fsc (forward scatter)
#' @param y an object; the channel to plot on the y axis; defaults to ssc (side scatter)
#' @param batch a logical of whether to batch process or not. if TRUE, adds additional facets by sample id. defaults to FALSE.
#' @description 
#' a wrapper for making a ggplot object of the density of cells on two axes. creates a density column that estimates distance of each observation from the center of the channel. the product of density on both channels defines the density of cells. 
plot.density <- function(data, x = fsc, y = ssc){
  data %>%
    tidyr::pivot_wider(names_from = ".ch",
                values_from = ".val") %>%
    dplyr::mutate(.d.x = base::scale({{x}}, scale = F),
           .d.y = base::scale({{y}}, scale = F),
           .d.z = base::sqrt(base::abs(.d.x * .d.y))) %>%
    dplyr::mutate(dplyr::across(.cols = .d.x:.d.z,
                  .fns = ~base::abs(.x))) %>%
    ggplot2::ggplot(aes(x = {{x}},
               y = {{y}},
               color = .d.z,
               alpha = 0.1)) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::guides(color = "none", alpha = "none") +
    ggplot2::scale_color_gradient2(low = "darkred",
                          mid = "gold",
                          high = "darkblue",
                          midpoint = 0.5) +
    ggplot2::theme_classic()
}
