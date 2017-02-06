#==============================================================================
#The percent of the sample represented by each taxa per taxon level
#==============================================================================
#'The percentage each taxon makes up of a sample
#'
#'@param wide.df = Data in a wide data format.
#'@param taxa.rank = Taxonomic rank (all caps).
#'@param master.df = Master taxa list
#'@return The percent of the sample represented by each taxa per taxon level.
#'@export
#==============================================================================
calc_pct_taxa <- function(wide.df, taxa.rank, master.df = BIBI::master){
  wide.df[, 6:ncol(wide.df)] <- (wide.df[, 6:ncol(wide.df)] /
                                   rowSums(wide.df[, 6:ncol(wide.df)])) * 100
  t_list <- unique(toupper(master.df[, taxa.rank]))
  shorten <- wide.df[, colnames(wide.df) %in% t_list]
  cn <- colnames(shorten)
  list_taxa.df <- if(length(cn) == 0){
    list_taxa.df <- data.frame(wide.df[, c("EVENT_ID", "STATION_ID",
                                           "DATE", "SAMPLE_NUMBER",
                                           "AGENCY_CODE")])
    list_taxa.df$NO_MATCH <- 0
    colnames(list_taxa.df) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                "AGENCY_CODE", "NO_MATCH")
    list_taxa.df
  }else{
    list_taxa.df <- data.frame(cbind(wide.df[, 1:5], shorten))
    colnames(list_taxa.df) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                "AGENCY_CODE", paste("PCT", cn, sep = "_"))
    list_taxa.df
  }
  final.df <- data.frame(list_taxa.df)
  return(data.frame(list_taxa.df[order(list_taxa.df$EVENT_ID),]))
}
#==============================================================================
#'The percentage each taxon makes up of a sample
#'
#'@param long.df = Data in a long data format
#'@param master.df = Master taxa list
#'@return The percent of the sample represented by each taxon per taxonomic level.
#'@export
#==============================================================================

seq_pct_taxa <- function(long.df, master.df = BIBI::master){
  
  phylum <- if("PHYLUM" %in% names(long.df)) wide(long.df, "PHYLUM")
  subphylum <- if("SUBPHYLUM" %in% names(long.df)) wide(long.df, "SUBPHYLUM")
  class <- if("CLASS" %in% names(long.df)) wide(long.df, "CLASS")
  subclass <- if("SUBCLASS" %in% names(long.df)) wide(long.df, "SUBCLASS")
  order <- if("ORDER" %in% names(long.df)) wide(long.df, "ORDER")
  suborder <- if("SUBORDER" %in% names(long.df)) wide(long.df, "SUBORDER")
  family <- if("FAMILY" %in% names(long.df)) wide(long.df, "FAMILY")
  subfamily <- if("SUBFAMILY" %in% names(long.df)) wide(long.df, "SUBFAMILY")
  tribe <- if("TRIBE" %in% names(long.df)) wide(long.df, "TRIBE")
  genus <- if("GENUS" %in% names(long.df)) wide(long.df, "GENUS")
  species <- if("SPECIES" %in% names(long.df)) wide(long.df, "SPECIES")
  
  pct_phylum <- if(length(phylum) > 0) calc_pct_taxa(phylum, taxa.rank = "PHYLUM")
  pct_subphylum <- if(length(subphylum) > 0) calc_pct_taxa(subphylum, "SUBPHYLUM")
  pct_class <- if(length(class) > 0) calc_pct_taxa(class, "CLASS")
  pct_subclass <- if(length(subclass) > 0) calc_pct_taxa(subclass, "SUBCLASS")
  pct_order <- if(length(order) > 0) calc_pct_taxa(order, "ORDER")
  pct_suborder <- if(length(suborder) > 0) calc_pct_taxa(suborder, "SUBORDER")
  pct_family <- if(length(family) > 0) calc_pct_taxa(family, "FAMILY")
  pct_subfamily <- if(length(subfamily) > 0) calc_pct_taxa(subfamily, "SUBFAMILY")
  pct_tribe <- if(length(tribe) > 0) calc_pct_taxa(tribe, "TRIBE")
  pct_genus <- if(length(genus) > 0) calc_pct_taxa(genus, "GENUS")
  pct_species <- if(length(species) > 0) calc_pct_taxa(species, "SPECIES")
  
  check_exists <- function(pct_taxa){
    pct_taxa <- if(length(pct_taxa) > 0){
      pct_taxa <- pct_taxa
    } else{
      pct_taxa <- pct_taxa[, 1:5]
    }
  }
  
  checked_pct_phylum <- check_exists(pct_phylum)
  checked_pct_subphylum <- check_exists(pct_subphylum)
  checked_pct_class <- check_exists(pct_class)
  checked_pct_subclass <- check_exists(pct_subclass)
  checked_pct_order <- check_exists(pct_order)
  checked_pct_suborder <- check_exists(pct_suborder)
  checked_pct_family <- check_exists(pct_family)
  checked_pct_subfamily <- check_exists(pct_subfamily)
  checked_pct_tribe <- check_exists(pct_tribe)
  checked_pct_genus <- check_exists(pct_genus)
  checked_pct_species <- check_exists(pct_species)
  
  #comb_all <- cbind(checked_pct_phylum, checked_pct_subphylum,
  #                  checked_pct_class, checked_pct_subclass,
  #                  checked_pct_order, checked_pct_suborder,
  #                  checked_pct_family, checked_pct_subfamily,
  #                  checked_pct_tribe, checked_pct_genus,
  #                  checked_pct_species)
  
  comb_all <- plyr::join_all(list(checked_pct_phylum, checked_pct_subphylum,
                                  checked_pct_class, checked_pct_subclass,
                                  checked_pct_order, checked_pct_suborder,
                                  checked_pct_family, checked_pct_subfamily,
                                  checked_pct_tribe, checked_pct_genus,
                                  checked_pct_species), by =  c("EVENT_ID", "STATION_ID",
                                                                "DATE", "SAMPLE_NUMBER",
                                                                "AGENCY_CODE"),
                             type = "full")
  comb_taxa <- comb_all[, 6:ncol(comb_all)]
  rm.cols <- names(comb_taxa[, colSums(comb_taxa) == 0])
  final.df <- comb_all[, !(names(comb_all) %in% rm.cols)]
  
  return(final.df)
}
