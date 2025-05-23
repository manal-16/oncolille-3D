% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flamingo.reconstruct_structure_worker.R
\name{flamingo.reconstruct_structure_worker}
\alias{flamingo.reconstruct_structure_worker}
\title{flamingo.reconstruct_structure_worker}
\usage{
flamingo.reconstruct_structure_worker(
  input_if,
  pd,
  sw,
  lambda,
  max_dist,
  nThread = 2
)
}
\arguments{
\item{input_if}{Input interaction frequency matrix.}

\item{pd}{Input pairwise distance matrix.}

\item{sw}{Downsample rates, suggesting the fraction of the observed data to be used in the model.}

\item{lambda}{The penalty term for the diagnal entries.}

\item{max_dist}{The maximum allowed distance for the diagnal entries.}

\item{nThread}{Number of thread to be used in the model, default is 28.}
}
\value{
A flamingo_prediction object containing the fragment id and 3D coordinates
}
\description{
Reconstruct the 3D genome structure using low-rank tensor completion
}
\keyword{FLAMINGO}
