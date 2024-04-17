
#' @title plot.wave
#' @description wrapper for making a ggplot object that plots the fluorescent intensities in each channel of a dataframe as density plots (waves). 
#' @param data a dataframe of flow cytometry data
#' @examples \dontrun{
#' }
#' @export

plot.wave <- function(data){
  data %>%
    ggplot2::ggplot(aes(x = .val, fill = .ch, alpha = 0.1)) +
    ggplot2::geom_density() +
    ggplot2::facet_grid(rows = vars(.ch)) +
    ggplot2::guides(alpha = "none", fill = "none") +
    ggplot2::theme_classic()
}