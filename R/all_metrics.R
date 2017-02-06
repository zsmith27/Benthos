#==============================================================================
# Calculate All Metrics
#==============================================================================
#'Run All Metrics
#'
#'@param master.df = The master taxa list contains taxonomic ranks from Phylum
#'to species and known taxonomic attributes.
#'@param long.df = Taxonomic counts in a long data format.
#'@param taxa.rank = The lowest taxonomic rank ("ORDER", "FAMILY", or "GENUS")
#' used to calculate the metrics.  If the majority of your taxa are identified
#' to the family level, then it would be inappropriate to perform metric
#' calculations at the genus level.
#'@return Calculates all of applicable and available metrics in the package.
#'@export
#'
all_metrics <- function(long.df, master.df, taxa.rank) {
  #Prep==========================================================================
  # Wide format data frames for necessary taxonomic levels.
  names(long.df) <- trimws(toupper(names(long.df)))

  class.wide <- wide(long.df, "CLASS")
  order.wide <- wide(long.df, "ORDER")
  family.wide <- wide(long.df, "FAMILY")
  if(taxa.rank %in% c("TRIBE", "GENUS", "SPECIES")){
    tribe.wide <- wide(long.df, "TRIBE")
  }else{
    tribe.wide <- NULL
  }

  if(taxa.rank %in% c("GENUS")){
    genus.wide <- wide(long.df, "GENUS")
  }else{
    genus.wide <- NULL
  }
  # Specified Taxonomic level
  if(taxa.rank == "FAMILY"){
    taxa.rank.df <- family.wide
  }
  if(taxa.rank == "GENUS"){
    taxa.rank.df <- genus.wide
  }
  if(!(taxa.rank %in% c("FAMILY", "GENUS"))){
    taxa.rank.df <- wide(long.df, taxa.rank)
  }




  #============================================================================

  metrics <- data.frame(taxa.rank.df[, c("EVENT_ID", "STATION_ID", "DATE",
                                         "SAMPLE_NUMBER", "AGENCY_CODE")])
  colnames(metrics) <- c("EVENT_ID", "STATION_ID", "DATE",
                         "SAMPLE_NUMBER", "AGENCY_CODE")
  metrics <- metrics[order(metrics$EVENT_ID), ]
  #============================================================================

  # Calculate diversity metrics
  print("Calculating Diversity Metrics:")
  print("...Richness")
  metrics$RICH <- vegan::specnumber(taxa.rank.df[, 6:ncol(taxa.rank.df)])
  print("...Shannon Diversity")
  metrics$SHANNON <- shannon(taxa.rank.df)
  #print("  Effecitiv Richness (Shannon)")
  #metrics$EFFECTIVE_RICH_SHANNON <- effective_richness(taxa.rank.df, "shannon")
  print("...Simpson Diversity")
  metrics$SIMPSONS <- simpsons(taxa.rank.df)
  #print("  Effective Richness (Simpson)")
  #metrics$EFFECTIVE_RICH_SIMPSON <- effective_richness(taxa.rank.df, "invsimpson")
  print("...Hurlberts PIE")
  metrics$HURLBERTS_PIE <- hurlberts_pie(taxa.rank.df)
  print("...Margalefs Diversity")
  metrics$MARGALEFS <- margalefs(taxa.rank.df)
  print("...Menhinicks Diversity")
  metrics$MENHINICKS <- menhinicks(taxa.rank.df)
  print("...Pielou Evenness")
  metrics$PIELOU <-  pielou(taxa.rank.df)

  if(taxa.rank %in% c("FAMILY", "GENUS")){
    print("...Ephemeroptera Richness")
    metrics$RICH_EPHEMEROPTERA <- rich_ephemeroptera(long.df, taxa.rank)
    print("...Plecoptera Richness")
    metrics$RICH_PLECOPTERA <- rich_plecoptera(long.df, taxa.rank)
    print("...Trichoptera Richness")
    metrics$RICH_TRICHOPTERA <- rich_trichoptera(long.df, taxa.rank)
    print("...EPT Richness")
    metrics$RICH_EPT <- rich_ept(long.df, taxa.rank)
    print("...%EPT Richness")
    metrics$PCT_EPT_RICH <- pct_ept_rich(long.df, taxa.rank)

    # FFGs
    print("...Collector Richness")
    metrics$RICH_COLLECT <- rich_attribute(taxa.rank.df, master.df,
                                           attribute.column = "BIBI_FFG",
                                           attribute.interest = c("CG", "CF"), taxa.rank)
    print("...Gatherer Richness")
    metrics$RICH_GATHER <- rich_attribute(taxa.rank.df, master.df,
                                          attribute.column = "BIBI_FFG",
                                          attribute.interest = "CG", taxa.rank)
    print("...Filter Feeder Richness")
    metrics$RICH_FILTER <- rich_attribute(taxa.rank.df, master.df,
                                          attribute.column = "BIBI_FFG",
                                          attribute.interest = "CF", taxa.rank)
    print("...Shredder Richness")
    metrics$RICH_SHRED <- rich_attribute(taxa.rank.df, master.df,
                                         attribute.column = "BIBI_FFG",
                                         attribute.interest = "SH", taxa.rank)
    print("...Scraper Richness")
    metrics$RICH_SCRAPE <- rich_attribute(taxa.rank.df, master.df,
                                          attribute.column = "BIBI_FFG",
                                          attribute.interest = "SC", taxa.rank)
    print("...Predator Richness")
    metrics$RICH_PREDATOR <- rich_attribute(taxa.rank.df, master.df,
                                            attribute.column = "BIBI_FFG",
                                            attribute.interest = "PR", taxa.rank)
    # Habits
    print("...Climber Richness")
    metrics$RICH_CLIMB <- rich_attribute(taxa.rank.df, master.df,
                                         attribute.column = "BIBI_HABIT",
                                         attribute.interest = "CB", taxa.rank)
    print("...Swimmer Richness")
    metrics$RICH_SWIM <- rich_attribute(taxa.rank.df, master.df,
                                        attribute.column = "BIBI_HABIT",
                                        attribute.interest = "SW", taxa.rank)
    print("...Clinger Richness")
    metrics$RICH_CLING <- rich_attribute(taxa.rank.df, master.df,
                                         attribute.column = "BIBI_HABIT",
                                         attribute.interest = "CN", taxa.rank)
    print("...Burrower Richness")
    metrics$RICH_BURROW <- rich_attribute(taxa.rank.df, master.df,
                                          attribute.column = "BIBI_HABIT",
                                          attribute.interest = "BU", taxa.rank)
    print("...Sprawler Richness")
    metrics$RICH_SPRAWL <- rich_attribute(taxa.rank.df, master.df,
                                          attribute.column = "BIBI_HABIT",
                                          attribute.interest = "SP", taxa.rank)
    print("...Skater Richness")
    metrics$RICH_SKATE <- rich_attribute(taxa.rank.df, master.df,
                                         attribute.column = "BIBI_HABIT",
                                         attribute.interest = "SK", taxa.rank)
  }

  if(taxa.rank %in% "GENUS"){
    print("...Non-Chironomidae/Oligochaeta (NCO) Richness")
    metrics$RICH_NCO <- rich_nco(long.df, taxa.rank)
  }
  #============================================================================
  # Tolerance metrics need to be updated once more taxonomic levels are
  # added to the taxa attributes table

  if(taxa.rank %in% c("FAMILY", "GENUS")){
    print("...Intolerant Richness")
    metrics$RICH_INTOL <- rich_tolerance(taxa.rank.df, master.df, "BIBI_TV", 0, 3)
    print("...Moderately-Tolerant Richness")
    metrics$RICH_MODTOL <- rich_tolerance(taxa.rank.df, master.df, "BIBI_TV", 4, 6)
    print("...Tolerant Richness")
    metrics$RICH_TOL <- rich_tolerance(taxa.rank.df, master.df, "BIBI_TV", 7, 10)
    print("...Intolerant EPT Richness")
    metrics$RICH_EPT_NO_TOL <- ept_rich_no_tol(long.df, "FAMILY", master.df)
    print("...Intolerant %EPT Richness")
    metrics$PCT_EPT_RICH_NO_TOL <- pct_ept_rich_no_tol(long.df, "FAMILY", master.df)
  }

  if(taxa.rank %in% c("FAMILY")){
    print("...Becks Version 1")
    metrics$BECKS_V1 <- becks(taxa.rank.df, taxa.rank, master.df, beck.version = 1)
    print("...Becks Version 3")
    metrics$BECKS_V3 <- becks(taxa.rank.df, taxa.rank, master.df, beck.version = 3)
  }

  if(taxa.rank %in% c("GENUS", "SPECIES")){
    print("...Epeorus Richness")
    metrics$RICH_EPHEM_EPEORUS <- rich_ephem_epeorus(long.df, genus.wide)
  }

  taxa <- c("PHYLUM", "SUBPHYLUM", "CLASS",
            "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE",
            "GENUS", "SPECIES")

  master.fill <- fill_taxa(master.df)
  master.fill <- unique(master.fill[, c("FINAL_ID", taxa)])
  #test <- (master.fill[duplicated(master.fill$TSN_R), ])
  long.fill <- long.df
  long.fill <- long.fill[, !names(long.fill) %in% taxa]
  long.fill <- merge(long.fill, master.fill, by = "FINAL_ID", all.x = T)

  ord.fill <- wide(long.fill, "ORDER")
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    fam.fill <- wide(long.fill, "FAMILY")
    if(taxa.rank %in% "GENUS") gen.fill <- wide(long.fill, "GENUS")
  }

  if(taxa.rank %in% "ORDER") taxa.rank.fill <- ord.fill
  if(taxa.rank %in% "FAMILY") taxa.rank.fill <- fam.fill
  if(taxa.rank %in% "GENUS") taxa.rank.fill <- gen.fill
  if(!taxa.rank %in% c("ORDER", "FAMILY", "GENUS")) taxa.rank.fill <- wide(long.fill, taxa.rank)
  #============================================================================

  print("...%Dominant 1")
  metrics$PCT_DOM1 <- pct_dom(taxa.rank.df, 1)
  print("...%Dominant 2")
  metrics$PCT_DOM2 <- pct_dom(taxa.rank.df, 2)
  print("...%Dominant 3")
  metrics$PCT_DOM3 <- pct_dom(taxa.rank.df, 3)
  print("...%Dominant 4")
  metrics$PCT_DOM4 <- pct_dom(taxa.rank.df, 4)
  print("...%Dominant 5")
  metrics$PCT_DOM5 <- pct_dom(taxa.rank.df, 5)
  print("...Abundance")
  metrics$ABUNDANCE <- abundance(taxa.rank.df)

  #Tolerance Metrics============================================================
  # Tolerance metrics need to be updated once more tolerance values are added
  # to the taxa attributes table
  print("Calculating Tolerance Metrics:")
  if(taxa.rank %in% "FAMILY") {
    print("...ASPT")
    metrics$ASPT_MOD <- tol_index(long.fill, master.df, "ASPT", taxa.rank)
  }
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    print("...HBI")
    metrics$HBI <- tol_index(long.fill, master.df, "BIBI_TV", taxa.rank)
    print("...% Urban Intolerant")
    metrics$PCT_URBAN_INTOL <- pct_urban_intol(long.df, master.df)
    print("...% Intolerant (0-3)")
    metrics$PCT_INTOL_0_3 <- pct_tol_val(taxa.rank.fill, master.df, "BIBI_TV", 0, 3)
    print("...% Intolerant (0-4)")
    metrics$PCT_INTOL_0_4 <- pct_tol_val(taxa.rank.fill, master.df, "BIBI_TV", 0, 4)
    print("...% Moderately Tolerant (4-6)")
    metrics$PCT_MOD_TOL_4_6 <- pct_tol_val(taxa.rank.fill, master.df, "BIBI_TV", 4, 6)
    print("...% Tolerant (7-10)")
    metrics$PCT_TOLERANT_7_10 <- pct_tol_val(taxa.rank.fill, master.df, "BIBI_TV", 7, 10)
    print("...% Tolerant (5-10)")
    metrics$PCT_TOLERANT_5_10 <- pct_tol_val(taxa.rank.fill, master.df, "BIBI_TV", 5, 10)
  }

  #Composition Metrics==========================================================
  print("Calculating Composition Metrics:")
  if(taxa.rank %in% c("FAMILY", "GENUS")){
    print("...%EPT Minus Hydropsychidae and Baetidae")
    metrics$PCT_EPT_HYDRO_BAETID <- pct_ept_hydro_baetid(order.wide, family.wide)
    print("...%EPT Minus Hydropsychidae")
    metrics$PCT_EPT_NO_HYDRO <- pct_ept_no_hydro(order.wide, family.wide)
    print("...%Trichoptera Minus Hydropsychidae")
    metrics$PCT_NON_HYDROP_TRICHOPTERA <- pct_non_hydrop_trichoptera(order.wide, family.wide)
    #metrics$PCT_HYDRO_TRICHOPTERA <- pct_hydro_trichoptera(order.wide, family.wide)
    print("...%Ephemeroptera Minus Baetidae")
    metrics$PCT_EPHEMEROPTERA_NO_BAETID <- pct_epmeroptera_no_baetid(order.wide, family.wide)
    print("...%EPT Composed of Hydropsychidae")
    metrics$PCT_HYDRO_EPT <- pct_hydro_ept(order.wide, family.wide)
    #print("...%Retreat Caddisfly")
    suborder.wide <- wide(long.df, "SUBORDER")
    #metrics$PCT_RETREAT_CADDISFLY <- pct_taxon(suborder.wide, "ANNULIPALPIA")
    #print("...%Corbiculidae")
    #metrics$PCT_CORBICULIDAE <- pct_corbiculidae(family.wide)
    #metrics$PCT_CORBICULIDA_TEST <- pct_taxon(family.wide, "CORBICULIDAE")
    #print("...%Simuliidae")
    #metrics$PCT_SIMULIIDAE <- pct_simuliidae(family.wide)
    #metrics$PCT_SIMULIIDAE_TEST <- pct_taxon(family.wide, "SIMULIIDAE")
    #print("...%Chironomidae")
    #metrics$PCT_CHIRONOMIDAE <- pct_chironomidae(family.wide)
    #metrics$PCT_CHIRONOMIDAE_TEST <- pct_taxon(family.wide, "CHIRONOMIDAE")
    print("...%Oligochaeta and Chironomidae")
    metrics$PCT_OLIGO_CHIRO <- pct_oligo_chiro(class.wide, family.wide)
    phylum.wide <- wide(long.df, "PHYLUM")
    print("...%Annelida and Chironomidae")
    metrics$PCT_ANNELID_CHIRO <- pct_annelid_chiro(phylum.wide, family.wide)
    print("...%Limeston Taxa")
    metrics$PCT_LIMESTONE <- pct_limestone(order.wide, family.wide)
  }

  if(taxa.rank %in% c("GENUS")){
    print("...%Chironomidae Composed of Chricotopus and Chironomis")
    metrics$PCT_CC_CHIRO <- pct_cc_chironomidae(family.wide, genus.wide)
    print("...%EPT Composed of Cheumatopsyche")
    metrics$PCT_EPT_CHEUMATOPSYCHE <- pct_ept_cheumatopsyche(order.wide, genus.wide)
    print("...%EPT Composed of Hydropsyche")
    metrics$PCT_EPT_HYDROPSYCHE <- pct_ept_hydropsyche(order.wide, genus.wide)
    tribe.wide <- wide(long.df, "TRIBE")
    #print("...%Tanytarsini")
    #metrics$PCT_TANYTARSINI <- pct_tanytarsini(tribe.wide)
    #print("...%Orthocladiinae")
    #metrics$PCT_ORTHOCLADIINAE <- pct_orthocladiinae(long.df)
  }
  print("...%EPT")
  metrics$PCT_EPT <- pct_ept(order.wide)
  #print("...%Ephemeroptera")
  #metrics$PCT_EPHEMEROPTERA <- pct_ephemeroptera(order.wide)
  #print("...%Plecoptera")
  #metrics$PCT_PLECOPTERA <- pct_plecoptera(order.wide)
  #print("...%Trichoptera")
  #metrics$PCT_TRICHOPTERA <- pct_trichoptera(order.wide)

  #print("...%Coleoptera")
  #metrics$PCT_COLEOPTERA <- pct_coleoptera(order.wide)
  #print("...%Odonata")
  #metrics$PCT_ODONATA <- pct_odonata(order.wide)
  print("...%COTE")
  metrics$PCT_COTE <- pct_cote(order.wide)
  print("...%POTEC")
  metrics$PCT_POTEC <- pct_potec(order.wide)
  #print("...%Amphipoda")
  #metrics$PCT_AMPHIPODA <- pct_amphipoda(order.wide)
  #print("...%Bivalvia")
  #metrics$PCT_BIVALVIA <- pct_bivalvia(class.wide)
  #print("...%Unionoida")
  #metrics$PCT_UNIONOIDA <- pct_unionoida(order.wide)
  #print("...%Diptera")
  #metrics$PCT_DIPTERA <- pct_diptera(order.wide)
  print("...GOLD")
  metrics$GOLD <- gold(class.wide, order.wide)
  #print("...%Oligochaeta")
  #metrics$PCT_OLIGOCHAETA <- pct_oligochaeta(class.wide)
  #print("...%Non-Insect")
  #metrics$PCT_NON_INSECT <- pct_non_insect(class.wide)


  #Functional Feeding Group (FFG) Metrics========================================
  print("Calculating FFG Metrics:")
  #if(!(taxa.rank %in% "FAMILY") | taxa.rank %in% "FAMILY") taxa.rank.att <- "FAMILY"
  if(taxa.rank %in% c("ORDER")) taxa.rank.att <- "ORDER"
  if(taxa.rank %in% c("FAMILY", "SUBFAMILY", "TRIBE")) taxa.rank.att <- "FAMILY"
  # Update once genus attributes added
  if(taxa.rank %in% c("GENUS", "SPECIES")) taxa.rank.att <- "GENUS"
  #metrics$PCT_FFG_UNASSIGNED <- pct_attribute(taxa.rank.df, master.df, "BIBI_FFG", "UA", taxa.rank.att)
  #metrics$PCT_DOM_FFG <- pct_dom1_group(long.df, master.df, "BIBI_FFG", taxa.rank.att)
  print("...%Collector")
  metrics$PCT_COLLECT <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG", c("CG", "CF"), taxa.rank.att)
  print("...%Filter Feeder")
  metrics$PCT_FILTER <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG", "CF", taxa.rank.att)
  print("...%Gatherer")
  metrics$PCT_GATHER <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG", "CG", taxa.rank.att)
  print("...%Predator")
  metrics$PCT_PREDATOR <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG", "PR", taxa.rank.att)
  print("...%Scraper")
  metrics$PCT_SCRAPE <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG", "SC", taxa.rank.att)
  print("...%Shredder")
  metrics$PCT_SHRED <- pct_attribute(taxa.rank.fill, master.df, "BIBI_FFG", "SH", taxa.rank.att)

  #Habit Metrics=================================================================
  print("Calculating Habit Metrics:")
  #metrics$PCT_DOM_Habit <- pct_dom1_group(long.df, master.df, "BIBI_HABIT", taxa.rank.att)
  #metrics$PCT_HABIT_UNASSIGNED <- pct_attribute(taxa.rank.df, master.df, "BIBI_HABIT", "UA", taxa.rank.att)
  print("...%Burrower")
  metrics$PCT_BURROW <- pct_attribute(taxa.rank.fill, master.df, "BIBI_HABIT", "BU", taxa.rank.att)
  print("...%Climber")
  metrics$PCT_CLIMB <- pct_attribute(taxa.rank.fill, master.df, "BIBI_HABIT", "CB", taxa.rank.att)
  print("...%Clinger")
  metrics$PCT_CLING <- pct_attribute(taxa.rank.fill, master.df, "BIBI_HABIT", "CN", taxa.rank.att)
  print("...%Skater")
  metrics$PCT_SKATE <- pct_attribute(taxa.rank.fill, master.df, "BIBI_HABIT", "SK", taxa.rank.att)
  print("...%Sprawler")
  metrics$PCT_SPRAWL <- pct_attribute(taxa.rank.fill, master.df, "BIBI_HABIT", "SP", taxa.rank.att)
  print("...%Swimmer")
  metrics$PCT_SWIM <- pct_attribute(taxa.rank.fill, master.df, "BIBI_HABIT", "SW", taxa.rank.att)

  #print("...%Unidentified Taxa")
  #metrics$PCT_UNIDENTIFIED <- pct_unidentified(taxa.rank.df)
  #============================================================================
  print("...Sequence % Taxa")
  seq.pct <- seq_pct_taxa(long.df, master.df)
  merge.cols <- c("EVENT_ID", "STATION_ID", "AGENCY_CODE", "DATE", "SAMPLE_NUMBER")
  final.df <- merge(metrics, seq.pct, by = merge.cols)
  
  return(final.df)

}

