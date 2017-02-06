#==============================================================================
# METRICS FOR ASSESSING TAXA BY PREDEFINED GROUPS
# Functional Feeding Groups, Habitat, Composition, etc.
#==============================================================================
#'Data Frame of a specific list of taxa
#'
#'@param NameList = uninque list of taxa.
#'@param Taxa.df = Wide data frame format of taxonomic counts.
#'@return Data frame of a specific list of taxa. Used in habit and functional
#'feeding group (FFG) fuctions to extract only the taxa specified in the object
#'NameList from the wide data frame of taxa (i.e. Family taxa.rank or Genus taxa.rank)
#'@export

group_taxa <- function(NameList, Taxa.df){
  ID <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  taxa.list <- as.character(unlist(NameList))
  #taxa.list[which(c(1, diff(taxa.list)) != 0)]
  #idx <- match(taxa.list, names(Taxa.df))
  #idx <- idx[! is.na(idx)]
  taxa_list.df <- data.frame(Taxa.df[, names(Taxa.df) %in% c(ID, taxa.list)])
  taxa_list.df <- taxa_list.df[, !(names(taxa_list.df) %in% "UNIDENTIFIED")]
  taxa_list.df[is.na(taxa_list.df)] <- 0 #NA = zero
  if(ncol(taxa_list.df) < 6) {
    final.vec <- 0
  } else {
    if(ncol(taxa_list.df) == 6) {
      final.vec <- taxa_list.df[, 6]
    } else {
      if(ncol(taxa_list.df) > 6)
       final.vec <- rowSums(taxa_list.df[, 6:ncol(taxa_list.df)])
    }
  }
  return(final.vec)
}

#==============================================================================
#'Vector of taxa richness for a specific list of taxa
#'
#'@param NameList = uninque list of taxa.
#'@param Taxa.df = Wide data frame format of taxonomic counts.
#'@return A vector of taxa richness for a specific list of taxa representing
#'each sampling event. NameList from the wide data frame of taxa
#' (i.e. Family taxa.rank or Genus taxa.rank)
#'@export

group_rich <- function(NameList, Taxa.df){
  ID <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  taxa.list <- as.character(unlist(NameList))
  taxa_list.df <- data.frame(Taxa.df[, names(Taxa.df) %in% c(ID, taxa.list)])
  taxa_list.df <- taxa_list.df[, !(names(taxa_list.df) %in% "UNIDENTIFIED")]
  taxa_list.df[is.na(taxa_list.df)] <- 0 #NA = zero
  if(ncol(taxa_list.df) < 6) {
    final.vec <- 0
  } else {
    final.vec <- vegan::specnumber(taxa_list.df[, 6:ncol(taxa_list.df)])
  }
  return(final.vec)
}
#==============================================================================
#Only taxa present in the data are used to create each data frame
#Therefore, absent taxa do not have a column of counts for each station
#The metrics require a value for each variable called on
#This function asks if a taxon is present in the data frame and thus present in at least one of your stations
#If the taxon is not present a temporary column of zeros is formed to represent the taxon for a specific metric
#'Temporary vector with zero values
#'

#'@param taxon = The name of the taxon necessary for calculating a metric but
#'is not found in the data set.
#'@param taxa.rank = Taxanomic taxa.rank (e.g. Class, Order, Family, Genus, etc.).
#'@return Creates a temporary vector containing zero values the length of the
#'data frame specified by the object taxa.rank.  Prevents errors created by missing taxon names
#'necessary for metric calculations. For example, if no Trichoptera taxa were observed
#'in the data set then there would be no column named "Trichoptera" in the wide order
#'taxa.rank data frame.  When metrics, such as pct_ept or pct_trichoptera, search for a column
#'named "Trichoptera" a null value is returned and the function fails. Therefore,
#'this function temporarly fills the missing columns with zeros and the metric
#'can be calculated.
#'@export

blank_col <- function(taxon, taxa.rank){
  if(taxon %in% colnames(taxa.rank)) {
    final.vec <- taxa.rank[, taxon] 
  } else {
    final.vec <- rep(0, nrow(taxa.rank))
  }
  return(final.vec)
}

#==============================================================================
#'The percent of a group
#'
#'@param Taxa.df = Wide data frame format of taxonomic counts.
#'@param master.df.df = taxonomic attributes table.
#'@param Group = The taxonomic group to be assessed
#'@param Group_taxa.rank = The specific taxa.rank or taxa.ranks of the group to be assessed.
#'@param taxa.rank = The taxonomic taxa.rank used during the assessment.
#'@return The percentage of taxa representing a predefined group. Typically,
#'this function is used to assess functional feeding group and habits.
#'@export

pct_attribute <- function(Taxa.df, master.df, Group, Group_taxa.rank, taxa.rank = "FAMILY"){
  #split.taxa <- split(master.df[, taxa.rank], master.df[, Group])
  #name.list <- split.taxa[Group_taxa.rank]
  new.group <- c(Group_taxa.rank)
  grep.taxa <- master.df[grepl(paste(new.group,collapse="|"), master.df[, Group]), ]
  name.list <- as.list(unique(grep.taxa$FINAL_ID))
  group.taxa <- group_taxa(name.list, Taxa.df)
  if(sum(group.taxa) == 0){
    final.vec <- 0
  }else{
    final.vec <- (group.taxa / rowSums(Taxa.df[, 6:ncol(Taxa.df)])) * 100
  }
  return(final.vec)
}

#==============================================================================
#'The richness of a group
#'
#'@param taxa.wide = Wide data frame format of taxonomic counts.
#'@param master = taxonomic attributes table.
#'@param attribute.column = The name of the column that contains the
#'attribute of interest.
#'@param attribute.interest = The specific attribute of interest
#'@param rank = The taxonomic taxa.rank used during the assessment.
#'@return The richness of taxa representing a predefined group. Typically,
#'this function is used to assess functional feeding group and habits.
#'@export

rich_attribute <- function(taxa.wide, master.df = BIBI::master, attribute.column,
                           attribute.interest, rank = "FAMILY"){
  new.group <- c(attribute.interest)
  grep.taxa <- master.df[grepl(paste(new.group, collapse="|"), master.df[, attribute.column]), ]
  name.list <- as.list(grep.taxa$FINAL_ID)
  ID <- c("EVENT_ID", "STATION_ID", "DATE", "AGENCY_CODE", "SAMPLE_NUMBER")
  taxa.list <- as.character(unlist(name.list))
  group.rich <- group_rich(name.list, taxa.wide)
  final.vec <- group.rich
  #taxa_list.df <- data.frame(taxa.wide[, names(taxa.wide) %in% c(ID, taxa.list)])
  #taxa_list.df[is.na(taxa_list.df)] <- 0 #NA = zero
  #if(ncol(taxa_list.df) < 6) {
  #  final_taxa.wide <- 0
  #} else {
  #  final_taxa.wide <- vegan::specnumber(taxa_list.df[, 6:ncol(taxa_list.df)])
  #}
  return(final.vec)
}

#==============================================================================
#'The percent of the most dominant group
#'
#'@param long.df = long.df data frame format of taxonomic counts.
#'@param master.df = taxonomic attributes table.
#'@param Group = The taxonomic group to be assessed
#'@param taxa.rank = The taxonomic taxa.rank used during the assessment.
#'@return The percentage of taxa representing by the most dominant (abundant)
#'group. Typically, this function is used to assess the functional feeding
#'groups and habits.
#'@export

pct_dom1_group <- function(long.df, master.df, Group, taxa.rank){
  taxa.info <- master.df[, c("FINAL_ID", Group)]
  merged <- merge(long.df, taxa.info, by.x = taxa.rank, by.y = "FINAL_ID",
                  all.x = TRUE)
  wide.df <- wide(merged, Group)
  final.vec <- pct_dom(wide.df, 1)
  return(final.vec)
}
