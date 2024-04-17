
#' @title tidy.fcs
#' @description a tidyr pipeline for standardizing .fcs data
#' @param data a dataframe read in from .fcs
#' @return a dataframe with standardized column names and a row for each observation
#' @examples \dontrun{
#' }
#' @export

tidy.fcs <- function(data){
  data %>%
    dplyr::rename_with(.cols = tidyselect::everything(),
                .fn = stringr::str_to_lower) %>%
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