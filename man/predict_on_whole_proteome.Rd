% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_whole_proteom.R
\name{predict_on_whole_proteome}
\alias{predict_on_whole_proteome}
\title{A function to predict and annotate whole proteom PTM sites}
\usage{
predict_on_whole_proteome(ptm_site, flanking_size = 12, SPIDER = T,
  positive_info_file, protein_fasta_file, n_fold = 2, lower_bound = -1,
  upper_bound = 1, liblinear_dir, feature_file_path, cvlog_path_name,
  specificity_level = 0.99, output_label)
}
\arguments{
\item{ptm_site}{The target amino acid of the given PTM type, in upper-case single letter representation.}

\item{flanking_size}{The number of residues surrounding each side of the center residue, The total window size will be 2*flanking_size+1 (default to 12).}

\item{SPIDER}{A boolean variable indicating whether to use SPIDER3 features (default set to TRUE.)}

\item{positive_info_file}{A text file containing the positive PTM sites in the required format.}

\item{protein_fasta_file}{A text file containing the protein sequences of interest in fasta format.}

\item{n_fold}{The number of folds used for training and prediction in cross validation stage (default set to 2).}

\item{lower_bound}{The lower bound of the scaled data range (default to -1).}

\item{upper_bound}{The upper bound of the scaled data range (default to 1).}

\item{liblinear_dir}{The path for the Liblinear tool.}

\item{feature_file_path}{The path for the feature files.}

\item{cvlog_path_name}{The path and name of the log files, which hold the details of Liblinear procedures.}

\item{specificity_level}{A number ranges from 0 to 1 indicating the specificity user requires the classifier to achieve (default to 0.99).}

\item{output_label}{The string to tag the output files.}
}
\description{
A function to predict and annotate whole proteom PTM sites
}
\details{
This function outputs the features generated from input files.
}
\examples{
predict_on_whole_proteome(ptm_site = "S",
                         flanking_size = 12,
                         SPIDER = T,
                         positive_info_file = "ps_PSP.tsv",
                         protein_fasta_file = "S_sp_fasta.tsv",
                         n_fold = 2,
                         lower_bound = -1,
                         upper_bound = 1,
                         liblinear_dir = "/data/ginny/liblinear-2.11/",
                         feature_file_path = "/data/ginny/test_package/",
                         cvlog_path_name = "/data/ginny/test_package/cvlog.txt",
                         specificity_level = 0.99,
                         output_label = "ps_0103")
}
