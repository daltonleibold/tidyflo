
#----
#' @title process.fcs
#' @description a pipeline for setting a standard processing protocol on flow cytometry data. log transforms all values, filters offscale events, then gates to most "average" cells in the population - within 1 s.d. of the median in all channels.
#' @param data a dataframe
#' @return a dataframe
#' @examples \dontrun{
#' }
#' @export 
process.fcs <- function(data){
  data %>% 
    dplyr::mutate(.val = log.sc(.val)) %>%
    tidyr::pivot_wider(names_from = ".ch",
                values_from = ".val") %>%
    dplyr::filter(dplyr::if_all(.cols = -c(.id, .n, .param),
                  .fns = ~.x >= 0 & .x <= 5)) %>%
    dplyr::group_by(.id) %>%
    dplyr::filter(dplyr::if_all(.cols = -base::c(.n, .param),
                  .fns = ~.x <= stats::median(.x) + (stats::median(.x) * r.sd(.x)) &
                    .x >= stats::median(.x) - (stats::median(.x) * r.sd(.x)))) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = -base::c(.id, .n, .param),
                 names_to = ".ch",
                 values_to = ".val")
}