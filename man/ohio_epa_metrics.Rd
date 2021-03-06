% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ohio_epa_metrics.R
\name{ohio_epa_metrics}
\alias{ohio_epa_metrics}
\title{Ohio Environmental Protection Agency Metrics}
\usage{
ohio_epa_metrics(long.df, master.df)
}
\arguments{
\item{long.df}{Taxonomic counts in a long data format.}

\item{master.df}{The master taxa list contains taxonomic ranks from Phylum
to species and known taxonomic attributes.}
}
\value{
Ten metrics were selected by the Ohio Environmental Protection Agency
to create the Ohio River macorinverterbate IBI. The index was developed
with macorinvertebrate data collected using an artificial substrate 
sampling methodology. Therefore, these metrics may not work well
with other sampling methodologies, even if the samples were collected
from the Ohio River. 
Ohio Environmental Protection Agency. 1988. 
Biological criteria for the protection of aquatic life: Volume II: 
users manual for biological field assessment of Ohio surface waters. 
Ohio EPA, Division of Water Quality Monitoring and Assessment, 
Surface Water Section, Columbus.
}
\description{
Ohio Environmental Protection Agency Metrics
}
