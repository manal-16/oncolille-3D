% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iflamingo.main_func.R
\name{iflamingo.main_func}
\alias{iflamingo.main_func}
\title{iflamingo.main_func}
\usage{
iflamingo.main_func(
  hic_data_low,
  file_format,
  epigenomic_file,
  domain_res,
  frag_res,
  chr_size,
  chr_name,
  normalization,
  downsampling_rates,
  lambda,
  lambda_epi,
  max_dist,
  nThread,
  alpha = -0.25,
  max_iter = 500,
  hic_data_high = NULL,
  norm_low = NULL,
  norm_high = NULL
)
}
\arguments{
\item{hic_data_low}{Input Hi-C data in .hic format, .mcool format or sparse matrix format.}

\item{file_format}{Foramt of the input hic data. Could be .hic, .mcool or sparse matrix format.}

\item{epigenomic_file}{Path to the epigenimic data in bedgraph format.}

\item{domain_res}{Size of the domains in bps, e.g. 1e6.}

\item{frag_res}{Size of the fragment in bps, e.g. 5e3.}

\item{chr_size}{Size of the chromosome in bps.}

\item{chr_name}{Name of the chromosome, e.g. chr1.}

\item{normalization}{Normalization method.}

\item{downsampling_rates}{Fraction of contacts to be used during the reconstruction.}

\item{lambda}{Lagrigian coefficient.}

\item{lambda_epi}{Lagrigian coefficient for the epigenomic term.}

\item{max_dist}{Maximum allowed distance betwee two consecutive points.}

\item{nThread}{Number of thread avalable for the reconstruction.}

\item{alpha}{Convertion factor between interaction frequency and pairwise distance. default -0.25.}

\item{max_iter}{Maximum iteration for the assembling algorithm. default 500.}

\item{hic_data_high}{Optional. The high resolution HiC data in sparse matrix format. Only required if the file_format is 'sparse matrix'.}

\item{norm_low}{Optional. The normalization vector for the low resolution Hi-C data. Only required if the file_format is 'sparse matrix'.}

\item{norm_high}{Optional. The normalization vector for the high resolution Hi-C data. Only required if the file_format is 'sparse matrix'.}
}
\value{
A data.frame containing the FLAMINGO predicted 3D structure.
}
\description{
Main function of iFLAMINGO, a wraper for all steps
}
\keyword{FLAMINGO}
