#==============================================================================
# Composition Metrics
#==============================================================================
#'Percentage of a Specified Taxon
#'@description Calculate the percentage of each sample represented by the
#'specified taxon or taxa.
#'@param long.df Taxonomic counts arranged in a long data format.
#'@param unique.id.col The name of the column that contains a unique sampling
#'event ID.
#'@param count.col The name of the column that contains taxanomic counts.
#'@param taxon.col The name of the column that contains the taxon or taxa 
#'of interest.
#'@param taxon The taxon or taxa of interest. To specify more than one taxa 
#'use: c("TAXA1", "TAXA2", "TAXA3").
#'@return A numeric vector of percentages.
#'@export

pct_taxon <- function(long.df, unique.id.col, count.col, taxon.col, taxon) {
  # Prep.
  unique.id.col <- enquo(unique.id.col)
  taxon.col <- enquo(taxon.col)
  count.col <- enquo(count.col)
  #----------------------------------------------------------------------------
  # Calculate the percentage of the specified taxon.
  final.vec <- long.df %>% 
    group_by(!!unique.id.col) %>% 
    summarise(TOTAL = sum(!!count.col),
              INDV = sum(UQ(count.col)[UQ(taxon.col) %in% taxon]),
              PCT = INDV / TOTAL * 100) %>% 
    pull(PCT)
  #----------------------------------------------------------------------------
  return(final.vec)
}

#==============================================================================
#'Proportion of Gastropoda, Oligochaeta, and Dipteran Individuals
#'@description The percentage of individuals not represented by 
#' gastropods (Class: Gastropoda), oligochaetes (Class: Oligochaeta),
#'  and dipteran (Order: Diptera). This metric typically decreases with 
#'  degradation.
#'@param long.df Taxonomic counts arranged in a long data format.
#'@param unique.id.col The name of the column that contains a unique sampling
#'event ID.
#'@param count.col The name of the column that contains taxanomic counts.
#'@param class.col The name of the column with class-level taxonomic information.
#'@param order.col The name of the column with order-level taxonomic information.
#'@return A numeric vector containing the percentage of each sample 
#'represented by the specified taxon or taxa.
#'@return Numeric vector of percentages.
#'@export

gold <- function(long.df, unique.id.col, count.col, class.col, order.col) {
  # Prep.
  unique.id.col <- enquo(unique.id.col)
  count.col <- enquo(count.col)
  class.col <- enquo(class.col)
  order.col <- enquo(order.col)
  #----------------------------------------------------------------------------
  # Calculate the percentage of Gastropoda, Oligochaeta, and Diptera.
  god <- pct_taxon(long.df,
                   unique.id = UQ(unique.id.col),
                   count.col = UQ(count.col),
                   taxon.col = UQ(class.col),
                   taxon = c("GASTROPODA", "OLIGOCHAETA")) +
    pct_taxon(long.df, 
              unique.id = UQ(unique.id.col),
              count.col = UQ(count.col),
              taxon.col = UQ(order.col),
              taxon = "DIPTERA")
  #----------------------------------------------------------------------------
  # GOLD is the percentage of taxa not represented by Gastropoda, Oligochaeta,
  # and Diptera.
  final.vec <- (100 - god)
  #----------------------------------------------------------------------------
  return(final.vec)
}

#==============================================================================

