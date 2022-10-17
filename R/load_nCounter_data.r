#' Load nCounter files 
#'
#' This loads count and target file info (internal function only)
#' @param pathtoRCC  matrix raw or normalized (default is current working directory)
#' @param target_file file with experimental information (meta data csv file)
#' @keywords Ncounter RCC files
#' @export
#' @examples
#' load_nCounter_files(pathtoRCC="./", target_file="targetfile.csv")
#' 
load_nCounter_files <- function(pathtoRCC="./", target_file) {
    target = read.csv(target_file)

}