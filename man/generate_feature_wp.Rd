% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_whole_proteom.R
\name{generate_feature_wp}
\alias{generate_feature_wp}
\title{Generate feature data for whole proteome}
\usage{
generate_feature_wp(ptm_site, flanking_size = 12, SPIDER = T,
  positive_info_file, protein_fasta_file, output_label)
}
\arguments{
\item{ptm_site}{The amino acid this PTM involves, in upper-case single letter representation.}

\item{flanking_size}{The number of residues surrounding each side of the center residue, the total window size will be 2*flanking_size+1, default to 12.}

\item{SPIDER}{A boolean variable indicating the usage of SPIDER3 features, default set to TRUE.}

\item{positive_info_file}{A text file containing the positive PTM sites info in required format.}

\item{protein_fasta_file}{A text file containing the proteins sequences of interest in Fasta format.}

\item{output_label}{The string to tag the output files.}
}
\description{
This function generates feature data for whole proteome by mapping the positive PTM info to protein sequences, constructing windows and extracting 3 sets of features.
}
\details{
This function outputs the features generated from input files.
}
\examples{
generate_feature_wp(ptm_site = "S",
             flanking_size = 12,
             SPIDER = T,
             positive_info_file = "ps_PSP.tsv",
             protein_fasta_file = "S_sp_fasta.tsv",
             output_label = "ps_0103")
}
