
#-------------------
# data processing
#-------------------
#----
#' @title read.fcs
#' @description simple wrapper that reads a *.fcs to a dataframe
#' @param file a character string that indicates a path to a *.fcs file. If a directory is given, all *.fcs files in the directory will be read.
#' @return a dataframe that contains the filename, ID, channel name and corresponding values
#' @examples \dontrun{
#' }
#' @export

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