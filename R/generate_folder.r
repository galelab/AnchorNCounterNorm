#' generates folder to put results
#'
#' @param foldername  name of folder (directory) to generate
#' @examples
#' generate_folder("directoryname")
#' 
#' 


generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}