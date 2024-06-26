% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assemblage_v240228.R
\name{assemblage}
\alias{assemblage}
\title{assemblage}
\usage{
assemblage(
  y,
  x,
  bench = c(),
  x.weight = c(),
  x.rank = c(),
  train.size = 1,
  moving.average = 1,
  horizon.gap = 0,
  window.size = c(),
  step.estim.window = 1,
  standardize.values = TRUE,
  coef.sum.1 = TRUE,
  cores = 1,
  rank.OOS = c(),
  comp.OOS = c(),
  Bench.OOS = c(),
  pred.Y.OOS = c(),
  lambda.grid.C = c(),
  model.select = c("ALL"),
  all.sample.fit = FALSE,
  Progression = 0,
  Packages = TRUE
)
}
\arguments{
\item{y}{Target vector (no transformation are done)}

\item{x}{The matrix we need to transform and rank (Should be in month-over-month format)}

\item{bench}{The benchmark used for model comparison. It can be left empty. (default is empty) If provided, it should follow the month-over-month format and will undergo the transformation of x.}

\item{x.weight}{The raw weight matrix can be left empty, it can be a vector or a matrix. (default is empty) If used, the components should have the same order as those in 'x'.}

\item{x.rank}{The rank matrix already pre-sorted and follows the correct moving.average format (default is empty). If utilized, the y, x, bench, x.weight, and x.rank inputs must already be in the correct format, and the moving.average must be set as c().}

\item{train.size}{The endpoint of the training sample (default is 1) defines the beginning of the out-of-sample period for calculation of the RMSE. If train.size is less than or equal to 1 (i.e., a fraction), it represents the proportion of observations to use in the in-sample estimation. Alternatively, if it is greater than 1 (e.g., 240), it signifies the index of the last row used in the in-sample training set. For instance, if the dataset comprises 300 rows and train.size is set to 240, 80\% of the dataset would be used for the in-sample estimation.}

\item{moving.average}{This parameter governs the transformation of the x matrix. (default is 1) If set to 1, the x matrix will keep the Month over Month format. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 corresponds Year over Year analysis. If left empty (using c()), the y, x, bench, x.weight, and x.rank inputs are expected to be in the correct format beforehand. Opting for this choice reduces estimation time by bypassing the transformation step.}

\item{horizon.gap}{The number of observations to be excluded to ensure that there is no overlapping information within the out-of-sample dataset (default is 0)}

\item{window.size}{The size of the rolling window (default is all the in-sample)}

\item{step.estim.window}{Re estimating the rolling window model every ... observations (default is 1)}

\item{standardize.values}{Standardize the response variables (default is TRUE)}

\item{coef.sum.1}{Does the coefficients need to sum to one, only affects comp and Bench models. (default is TRUE)}

\item{cores}{Number of cores to run the code, helps for glmnet function and cross-validation. (default is 1)}

\item{rank.OOS}{Xrank matrix for out-of-sample forecast using the coefficients that used all the data set for training. (Matrix will undergo no transformation) (default is c())}

\item{comp.OOS}{Xcomp matrix for out-of-sample forecast using the coefficients that used all the data set for training. (Matrix will undergo no transformation) (default is c())}

\item{Bench.OOS}{Benchmark matrix for out-of-sample forecast using the coefficients that used all the data set for training. (Matrix will undergo no transformation) (default is c())}

\item{pred.Y.OOS}{Target for the OOS forecast. (Used to have an RMSE) (default is c())}

\item{lambda.grid.C}{Personalize lambda grid (default is a grid of 20 lambda)}

\item{model.select}{Select the model to run 'rank' = rank space ,'comp' = component space, 'blend' = blend rank/comp/bench ,'Bench'= benchmark without intercept, 'BInt' = benchmark with Intercept,'ALL' (default is 'ALL')}

\item{all.sample.fit}{Do an estimation using all the data set as training (default is FALSE)}

\item{Progression}{Give feedback on the estimation. (default 0) 0 = No feedback, 1 = Print each iteration, 2 = Print after every model estimation & each iteration}

\item{Packages}{Automatically update the required packages. (default TRUE)}
}
\value{
Return a list containing the loadings, features, fitted values and more. With the predictions and fitted values at the moving.average frequency specified by the user.
}
\description{
Assemblage Regression represents a specialized form of generalized nonnegative ridge regression, designed to optimize the weights of subcomponents to maximize the predictive capability of the aggregate.
}
\details{
When running the code, it's possible to utilize the bypassing option, ensuring all available observations are considered within the desired moving average unit. This entails saving various crucial variables: assemblage.info[\link{'x.info'}][\link{'x.ranks.ma'}] for the x.rank matrix, assemblage.info[\link{'x.info'}][\link{'x.comps.ma'}] for the x matrix, assemblage.info[\link{'x.info'}][\link{'x.weights'}] for the x.weight matrix, and assemblage.info[\link{'bench.info'}][\link{'x.bench.ma'}] for the bench matrix. It's important to note that y will retain its original matrix. Subsequently, executing a loop using the 'assemblage' function becomes feasible, utilizing the aforementioned variables and setting the 'moving.average' option as 'c()'.
}
\examples{
Check out the GitHub page to see example using real inflation data.
}
