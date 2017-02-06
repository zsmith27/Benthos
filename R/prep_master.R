#setwd("C:/Users/Owner/Desktop/Benthos/Benthos")
#data(master)
#master <- read.csv("MTL.csv", stringsAsFactors = FALSE)
#master$AGENCY_ID <- trimws(master$AGENCY_ID)
#master$AGENCY_ID <- gsub(" ","_", master$AGENCY_ID)
#check.master <- master[, c("AGENCY_ID", "AGENCY_ID2")]
#devtools::use_data(master, overwrite = TRUE)


#devtools::use_package("data.table")
#devtools::use_package("vegan")
#devtools::use_package("tidyr")
