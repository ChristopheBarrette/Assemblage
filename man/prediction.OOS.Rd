% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assemblage_v240228.R
\name{prediction.OOS}
\alias{prediction.OOS}
\title{prediction.OOS}
\usage{
prediction.OOS(
  assemblage.info,
  rank.OOS,
  comp.OOS,
  Bench.OOS,
  pred.Y.OOS,
  model.select
)
}
\arguments{
\item{assemblage.info}{Information from the previous steps}

\item{rank.OOS}{Rank matrix}

\item{comp.OOS}{Component matrix}

\item{Bench.OOS}{Benchmark matrix}

\item{pred.Y.OOS}{Target for the OOS forecast}

\item{model.select}{Select the model}
}
\description{
Function to do the OOS prediction using the coefficients that used all the data set for training.
}
