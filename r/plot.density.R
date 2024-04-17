

#' @title plot.density
#' @description a wrapper for making a ggplot object of the density of cells on two axes. creates a density column that estimates distance of each observation from the center of the channel. the product of density on both channels defines the density of cells. 
#' @param data a dataframe of flow cytometry data
#' @param x an object; the channel to plot on the x axis; defaults to fsc (forward scatter)
#' @param y an object; the channel to plot on the y axis; defaults to ssc (side scatter)
#' @param batch a logical of whether to batch process or not. if TRUE, adds additional facets by sample id. defaults to FALSE.
#' @return a ggplot object
#' @examples \donrun{plot.density(data = data)}
#' @export 
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
