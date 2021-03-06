% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composition_metrics.R
\name{pct_ept_hydro_baetid}
\alias{pct_ept_hydro_baetid}
\title{Percentage of EPT Individuals Minus Hydropsychid and Baetid Individuals}
\usage{
pct_ept_hydro_baetid(order.wide, family.wide)
}
\arguments{
\item{order.wide}{Taxonomic counts aggregated at the order level
classification in a wide data format. Use the wide function to
prepare the data.}

\item{family.wide}{Taxonomic counts aggregated at the family level
classification in a wide data format. Use the wide function to
prepare the data.}
}
\value{
The percentage of individuals identified as ephemeropterans
(Order: Ephemeroptera), plecopterans (Order: Plecoptera),
and trichopterans (Order: Trichoptera) (EPT) minus the percentage of
hydropsychids (Family: Hydropsychidae) and baetids (Family: Baetidae).
These two families are typically tolerant and can be hyperdominant,
resulting in elevated percentages of EPT that may not accurately
describe the degradation gradient of interest. This metric will
typically decrease with degradation.
}
\description{
Percentage of EPT Individuals Minus Hydropsychid and Baetid Individuals
}
