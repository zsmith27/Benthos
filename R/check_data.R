#==============================================================================
# Check Data
#==============================================================================
#'Check Data
#'
#'@param el.
#'@return 
#'@export
#'


check_data <- function(long.df, master.df){
  
  benthos.cols <- c("EVENT_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                "SAMPLE_NUMBER", "FINAL_ID", "REPORTING_VALUE")
  
  if(any((benthos.cols %in% names(long.df)) == FALSE)){
    missing.cols <- benthos.cols[!benthos.cols %in% names(long.df)]
    error.1 <- paste("You are missing the following column(s):",
                      paste(missing.cols, collapse = ", "))
    stop(error.1)
  }
  
  #if((!class(long.df[, c("EVENT_ID", "STATION_ID", "FINAL_ID")]) %in% c("character", "factor")) == FALSE){}
  
  #If(length(long.df[!long.df$FINAL_ID %in% master.df$FINAL_ID, ]) > 0) {
  #  missing.taxa <- sort(unique(long.df[!long.df$FINAL_ID %in% master.df$FINAL_ID, "FINAL_ID"]))
  #  error.2 <- paste0("The following taxon or taxa were not found in the Master Taxa List: ",
  #                   paste(missing.taxa, collapse = ", "), ". Check the spelling of each name,
  #                   an exact match is required.")
  #  print(error.2)
  #}
  
}
