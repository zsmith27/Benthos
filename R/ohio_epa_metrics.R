#'Ohio Environmental Protection Agency Metrics
#'
#'@param long.df Taxonomic counts in a long data format.
#'@param master.df The master taxa list contains taxonomic ranks from Phylum
#'to species and known taxonomic attributes.
#'@return Ten metrics were selected by the Ohio Environmental Protection Agency
#'to create the Ohio River macorinverterbate IBI. The index was developed
#'with macorinvertebrate data collected using an artificial substrate 
#'sampling methodology. Therefore, these metrics may not work well
#'with other sampling methodologies, even if the samples were collected
#'from the Ohio River. 
#'Ohio Environmental Protection Agency. 1988. 
#'Biological criteria for the protection of aquatic life: Volume II: 
#'users manual for biological field assessment of Ohio surface waters. 
#'Ohio EPA, Division of Water Quality Monitoring and Assessment, 
#'Surface Water Section, Columbus.
#'@export

ohio_epa_metrics <- function(long.df, master.df){
  # Reformate tables to make it easier to calculate metrics.
  class.wide <- Benthos::wide(long.df, "CLASS")
  ord.wide <- Benthos::wide(long.df, "ORDER")
  tribe.wide <- Benthos::wide(long.df, "TRIBE")
  gen.wide <- Benthos::wide(long.df, "GENUS")
  #============================================================================
  # Create a new data frame to store metric outputs.
  metrics.df <- gen.wide[, 1:7]
  metrics.df$SAMPLE_COUNT <- as.numeric(gsub("\\..*","", metrics.df$UNIQUE_ID))
  #============================================================================
  # Calculate Metrics
  #============================================================================
  # 1. Richness
  metrics.df$RICH <- vegan::specnumber(gen.wide[, 8:ncol(gen.wide)])
  # 2. Ephemeroptera Richness
  metrics.df$RICH_EPHEMEROPTERA <- Benthos::taxon_richness(long.df, "EPHEMEROPTERA", "ORDER", "GENUS")
  # 3. Trichoptera Richness
  metrics.df$RICH_TRICHOPTERA <- Benthos::taxon_richness(long.df, "TRICHOPTERA", "ORDER", "GENUS")
  # 4. Diptera Richness
  metrics.df$RICH_DIPTERA <- Benthos::taxon_richness(long.df, "DIPTERA", "ORDER", "GENUS")
  # 5. % Ephemeroptera
  metrics.df$PCT_EPHEMEROPTERA <- Benthos::pct_taxon(ord.wide, "EPHEMEROPTERA")
  # 6. % Trichoptera
  metrics.df$PCT_TRICHOPTERA <- Benthos::pct_taxon(ord.wide, "TRICHOPTERA")
  # 7. % Tanytarsini
  metrics.df$PCT_TANYTARSINI <- Benthos::pct_taxon(tribe.wide, "TANYTARSINI")
  # 8. % Non-Tanytarsini Diptera and Non-Insecta
  metrics.df$PCT_NON_TANYTARSINI_NON_INSECTA <- Benthos::pct_non_tanytarsini_non_insecta(class.wide, ord.wide, tribe.wide)
  # 9. % Tolerant
  # NEEDS WORK
  #metrics.df$PCT_TOLERANT_OHIO <- Benthos::ohio_pct_tol()
  
  # 10. Qualitative EPT Taxa
  # Not included because this sample was based on a qualitative 
  # sample collected in addition to the artificial susbstrate sample.
  # This metric is simply richness of EPT taxa but it does not seem appropriate
  # to include this metric becuase it was developed for a specific sampling
  # methodology.
  #============================================================================
  return(metrics.df)
}