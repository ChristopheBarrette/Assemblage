% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assemblage_v240228.R
\name{bench.transformation}
\alias{bench.transformation}
\title{bench.transformation}
\usage{
bench.transformation(bench, moving.average)
}
\arguments{
\item{bench}{The benchmark used for model comparison. It can be left empty. If provided, it should follow the month-over-month format and will undergo a transformation similar to x.}

\item{moving.average}{This parameter governs the transformation of the x matrix. If set to 1, the x matrix will be retained for Month over Month analysis. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 signifies Year over Year analysis.}
}
\description{
Function used to transform the bench matrix into the right moving average unit
}
