% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/domain_enrichment.R
\name{calculate_tbt_negative_ptms}
\alias{calculate_tbt_negative_ptms}
\title{Analyse negative cross talk events in domains.}
\usage{
calculate_tbt_negative_ptms(anchor_mod, cross_mod, anchor_mapped_df_Rds,
  cross_mapped_df_Rds, output_label)
}
\arguments{
\item{anchor_mod}{A string indicating the anchor modification.}

\item{cross_mod}{A string indicating the other modification.}

\item{anchor_mapped_df_Rds}{An Rds file containing the window score file of anchor PTM mapped domain.}

\item{cross_mapped_df_Rds}{An Rds file containing the window score file of the other PTM with mapped domain.}

\item{output_label}{The string to tag the output files.}
}
\description{
Test for enrichment of two user provided PTM types inside of domains.
}
\examples{
calculate_tbt_negative_ptms(anchor_mod = "ubi",
                            cross_mod = "acety",
                            anchor_mapped_df_Rds = "ubi_0103_mapped_df.Rds",
                            cross_mapped_df_Rds = "acety_0103_mapped_df.Rds",
                            output_label = "ubi_acety_negtive")
}
