#==============================================================================
#'Percentage of Amphipod Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as amphipods
#' (Order: Amphipoda).  This metric will typically increase with degradation.
#'@export

pct_amphipoda <- function(order.wide) {
  final.vec <- blank_col("AMPHIPODA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Chironomid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as chironomids
#'(Family: Chironomidae).
#'This metric will typically increase with degradation.
#'@export

pct_chironomidae <- function(family.wide) {
  final.vec <- blank_col("CHIRONOMIDAE", family.wide) /
    rowSums(family.wide[, 8:ncol(family.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Corbiculid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as corbiculids
#'(Family: Corbiculidae).
#'@export

pct_corbiculidae <- function(family.wide) {
  final.vec <- blank_col("CORBICULIDAE", family.wide) /
    rowSums(family.wide[, 8:ncol(family.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Bivalve Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as bivalves
#'(Class: Bivalvia).
#'@export

pct_bivalvia <- function(class.wide) {
  final.vec <- blank_col("BIVALVIA", class.wide) /
    rowSums(class.wide[, 8:ncol(class.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Unionoid Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as unionoids
#' (Order: Unionoida).
#'@export

pct_unionoida <- function(order.wide) {
  final.vec <- blank_col("UNIONOIDA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Dipteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as dipterans
#' (Order: Diptera). This metric will typically increase with degradation.
#'@export

pct_diptera <- function(order.wide) {
  final.vec <- blank_col("DIPTERA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Coleopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as coleopterans
#'(Order: Coleoptera).
#'@export

pct_coleoptera <- function(order.wide) {
  final.vec <- blank_col("COLEOPTERA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Odonate Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as odonates (Order: Odonata).
#'@export

pct_odonata <- function(order.wide) {
  final.vec <- blank_col("ODONATA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}


#==============================================================================
#'Percentage of Ephemeropteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as Ephemeropterans
#' (Order: Ephemeroptera).  This metric typically decreases with degradation.
#'@export
pct_ephemeroptera <- function(order.wide) {
  final.vec <- blank_col("EPHEMEROPTERA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Retreat-Making Trichopteran Individuals
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@return The percentage of individuals identified as retreat-making
#' trichopterans (Suborder: Annulipalpia).
#'@export
pct_retreat_trichoptera <- function(long) {
  sub.ord <- wide(long, "SUBORDER")
  final.vec <- blank_col("ANNULIPALPIA", sub.ord) /
    rowSums(sub.ord[, 8:ncol(sub.ord)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Hydropsychid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as hydropsychids
#' (Family: Hydropsychidae).
#'@export

pct_hydropsychidae <- function(family.wide) {
  final.vec <- blank_col("HYDROPSYCHIDAE", family.wide) /
    rowSums(family.wide[, 8:ncol(family.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Non-Insect Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals that were not identified as insects
#'(i.e., not identified as class Insecta).
#'@export

pct_non_insect <- function(class.wide){
  Non_Insect <- rowSums(class.wide[, 8:ncol(class.wide)]) -
    blank_col("INSECTA", class.wide)
  final.vec <- Non_Insect / rowSums(class.wide[, 8:ncol(class.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Oligochaet Individuals
#'
#'@param class.wide = Taxonomic counts aggregated at the class level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as oligochaets
#'(Class: Oligochaeta).
#'@export

pct_oligochaeta <- function(class.wide){
  final.vec <- blank_col("OLIGOCHAETA", class.wide) /
    rowSums(class.wide[, 8:ncol(class.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Plecopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as plecopterans
#'(Order: Plecoptera).
#'@export

pct_plecoptera <- function(order.wide) {
  final.vec <- blank_col("PLECOPTERA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Trichopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as trichopterans
#'(Order: Trichoptera).
#'@export

pct_trichoptera <- function(order.wide) {
  final.vec <- blank_col("TRICHOPTERA", order.wide) /
    rowSums(order.wide[, 8:ncol(order.wide)]) * 100
  return(final.vec)
}

#==============================================================================
#'Percentage of Non-Hydropsychid Trichopteran Individuals
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as trichopterans
#'(Order:Trichoptera) excluding hydropsychids (Family: Hydropsychidae).
#'@export

pct_non_hydrop_trichoptera <- function(order.wide, family.wide) {
  trichop.all <- blank_col("TRICHOPTERA", order.wide)
  trichop.no.hydro <-  trichop.all - blank_col("HYDROPSYCHIDAE", family.wide)
  final.vec <- ifelse(trichop.all == 0, 0, (trichop.no.hydro / trichop.all) * 100)
  return(final.vec)
}

#==============================================================================
#'Percentage of Simuliid Individuals
#'
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as simuliids
#'(Family: Simuliidae).
#'@export

pct_simuliidae <- function(family.wide) {
  final.vec <- blank_col("SIMULIIDAE", family.wide) /
    rowSums(family.wide[, 8:ncol(family.wide)]) * 100
  return(final.vec)
}


#============================================================================
#'Ephemeropteran Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. This
#'function requires a rank below the Order level taxonomic classification.
#'@return The number of taxa identified as ephemeropterans (Order: Ephemeroptera).
#'@export

rich_ephemeroptera <- function(long, rank = "FAMILY"){
  Order <- split(long[, rank], long$ORDER)
  taxa.list <- unique(Order$EPHEMEROPTERA)
  taxa.wide <- wide(long, rank)
  final.vec <- group_rich(taxa.list, taxa.wide)
  return(final.vec)
}

#============================================================================
#'Plecopteran Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. This
#'function requires a rank below the Order level taxonomic classification.
#'@return The number of taxa identified as plecopterans (Order: Plecoptera).
#'@export

rich_plecoptera <- function(long, rank = "FAMILY"){
  Order <- split(long[, rank], long$ORDER)
  taxa.list <- unique(Order$PLECOPTERA)
  taxa.wide <- wide(long, rank)
  final.vec <- group_rich(taxa.list, taxa.wide)
  return(final.vec)
}

#============================================================================
#'Trichopteran Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. This
#'function requires a rank below the Order level taxonomic classification.
#'@return The number of taxa identified as tichopterans (Order: Trichoptera).
#'@export

rich_trichoptera <- function(long, rank = "FAMILY"){
  Order <- split(long[, rank], long$ORDER)
  taxa.list <- unique(Order$TRICHOPTERA)
  taxa.wide <- wide(long, rank)
  final.vec <- group_rich(taxa.list, taxa.wide)
  return(final.vec)
}














