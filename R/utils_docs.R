#' A dummy function to store shared documentation items.
#' 
#' @import data.table
#' 
#' @param use_log_file logical. If TRUE, information about data processing
#' will be saved to a file.
#' @param append logical. If TRUE, information about data processing will be added
#' to an existing log file.
#' @param verbose logical. If TRUE, information about data processing wil be printed
#' to the console.
#' @param log_file_path character. Path to a file to which information about 
#' data processing will be saved. 
#' If not provided, such a file will be created automatically.
#' If `append = TRUE`, has to be a valid path to a file.
#' 
#' @return NULL.
#' @keywords internal
.documentFunction = function(fewMeasurements, 
                             useUniquePeptide,
                             summaryforMultipleRows, 
                             removeProtein_with1Feature,
                             removeProtein_with1Protein,
                             removeOxidationMpeptides,
                             removeMpeptides) {
  NULL
}