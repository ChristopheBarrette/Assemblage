

#' @title assemblage
#' @description Assemblage Regression represents a specialized form of generalized nonnegative ridge regression, designed to optimize the weights of subcomponents to maximize the predictive capability of the aggregate.
#' @param y Target vector (no transformation are done)
#' @param x The matrix we need to transform and rank (Should be in month-over-month format)
#' @param bench The benchmark used for model comparison. It can be left empty. (default is empty) If provided, it should follow the month-over-month format and will undergo the transformation of x.
#' @param x.weight The raw weight matrix can be left empty, it can be a vector or a matrix. (default is empty) If used, the components should have the same order as those in 'x'.
#' @param x.rank The rank matrix already pre-sorted and follows the correct moving.average format (default is empty). If utilized, the y, x, bench, x.weight, and x.rank inputs must already be in the correct format, and the moving.average must be set as c().
#' @param train.size The endpoint of the training sample (default is 1), which defines the beginning of the out-of-sample period for calculation of the RMSE.If train.size is less than or equal to 1, it represents the proportion of observations to use. If greater than 1, it signifies the index of the last row used in the in-sample training set.
#' @param moving.average This parameter governs the transformation of the x matrix. (default is 1) If set to 1, the x matrix will keep the Month over Month format. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 corresponds Year over Year analysis. If left empty (using c()), the y, x, bench, x.weight, and x.rank inputs are expected to be in the correct format beforehand. Opting for this choice reduces estimation time by bypassing the transformation step.
#' @param horizon.gap The number of observations to be excluded to ensure that there is no overlapping information within the out-of-sample dataset (default is 0)
#' @param window.size The size of the rolling window (default is all the in-sample)
#' @param step.estim.window Re estimating the rolling window model every ... observations (default is 1)
#' @param standardize.values Standardize the response variables (default is TRUE)
#' @param coef.sum.1 Does the coefficients need to sum to one (default is TRUE)
#' @param cores Number of cores to run the code, helps for glmnet function. (default is 1)
#' @param rank.OOS Xrank matrix for out-of-sample forecast using the coefficients that used all the data set for training. (Matrix will undergo no transformation) (default is c())
#' @param comp.OOS Xcomp matrix for out-of-sample forecast using the coefficients that used all the data set for training. (Matrix will undergo no transformation) (default is c())
#' @param Bench.OOS Benchmark matrix for out-of-sample forecast using the coefficients that used all the data set for training. (Matrix will undergo no transformation) (default is c())
#' @param pred.Y.OOS Target for the OOS forecast. (Used to have an RMSE) (default is c())
#' @param lambda.grid.C Personalize lambda grid (default is a grid of 20 lambda)
#' @param model.select Select the model to run 'rank' = rank space ,'comp' = component space, 'blend' = blend rank/comp/bench ,'Bench'= benchmark without intercept, 'BInt' = benchmark with Intercept,'ALL' (default is 'ALL')
#' @param all.sample.fit Do an estimation using all the data set as training (default is FALSE)
#' @param Progression Give feedback on the estimation. (default 0) 0 = No feedback, 1 = Print each iteration, 2 = Print after every model estimation & each iteration
#' @return Return a list containing the loadings, the features, the fitted values and more.
#' @details When running the code, it's possible to utilize the bypassing option, ensuring all available observations are considered within the desired moving average unit. This entails saving various crucial variables: assemblage.info[['x.info']][['x.ranks.ma']] for the x.rank matrix, assemblage.info[['x.info']][['x.comps.ma']] for the x matrix, assemblage.info[['x.info']][['x.weights']] for the x.weight matrix, and assemblage.info[['bench.info']][['x.bench.ma']] for the bench matrix. It's important to note that y will retain its original matrix. Subsequently, executing a loop using the 'assemblage' function becomes feasible, utilizing the aforementioned variables and setting the 'moving.average' option as 'c()'.
#' @rdname assemblage
#' @examples Check out the GitHub page to see example using real inflation data.
#' @export 
assemblage = function( y, x, bench=c(), x.weight=c(), x.rank=c(), train.size=1, moving.average=1,
                     horizon.gap=0, window.size=c(), step.estim.window = 1, standardize.values = TRUE, 
                     coef.sum.1 = TRUE, cores=1, rank.OOS=c(),comp.OOS=c(),Bench.OOS=c(), 
                     pred.Y.OOS=c(), lambda.grid.C=c(), model.select=c('ALL'),
                     all.sample.fit=FALSE, Progression=0){
 
# assemblage main function version 29-02-2024
  
################################################################
############################ Legend ############################
#
# y = Target vector (no transformation are done)
#
# x = The matrix we need to transform and rank (Should be in month-over-month format)
#
# bench = The benchmark used for model comparison. It can be left empty. (default is empty)
#         If provided, it should follow the month-over-month format and will undergo the transformation of x.
#
# x.weight = The raw weight matrix can be left empty, it can be a vector or a matrix. (default is empty)
#            If used, the components should have the same order as those in 'x'.
#
# x.rank = The rank matrix already pre-sorted and follows the correct moving.average format (default is empty). 
#          If utilized, the y, x, bench, x.weight, and x.rank inputs must already be in the correct format, 
#          and the moving.average must be set as c().
#
# train.size = The endpoint of the training sample (default is 1), which defines 
#              the beginning of the out-of-sample period for calculation of the RMSE. 
#              If train.size is less than or equal to 1, it represents the proportion 
#              of observations to use. If greater than 1, it signifies the index of 
#              the last row used in the in-sample training set.
#
# moving.average = This parameter governs the transformation of the x matrix. (default is 1)
#                  If set to 1, the x matrix will keep the Month over Month format. 
#                  Setting it to 3 corresponds to Quarter over Quarter analysis, 
#                  while 12 corresponds Year over Year analysis.
#                  If left empty (using c()), the y, x, bench, x.weight, and x.rank 
#                  inputs are expected to be in the correct format beforehand. 
#                  Opting for this choice reduces estimation time by bypassing the transformation step.
#
# horizon.gap = The number of observations to be excluded to ensure that there is no 
#               overlapping information within the out-of-sample dataset (default is 0)
#
# window.size = The size of the rolling window (default is all the in-sample)
#  
# step.estim.window = Re estimating the rolling window model every ... observations (default is 1)
#
# standardize.values = Standardize the response variables (default is TRUE)
#
# coef.sum.1 = Does the coefficients need to sum to one (default is TRUE)
#
# cores = Number of cores to run the code, helps for glmnet function. (default is 1)
#
# rank.OOS = Xrank matrix for out-of-sample forecast using the coefficients that 
#                  used all the data set for training. (Matrix will undergo no transformation) (default is c())
#
# comp.OOS = Xcomp matrix for out-of-sample forecast using the coefficients that 
#                  used all the data set for training. (Matrix will undergo no transformation) (default is c())
#
# Bench.OOS = Benchmark matrix for out-of-sample forecast using the coefficients that 
#                  used all the data set for training. (Matrix will undergo no transformation) (default is c())
#
# pred.Y.OOS = Target for the OOS forecast. (Used to have an RMSE) (default is c())
#
# lambda.grid.C = Personalize lambda grid (default is a grid of 20 lambda)
#
# model.select = Select the model to run 'rank' = rank space ,'comp' = component space,
#                'blend' = blend rank/comp/bench ,'Bench'= benchmark without intercept,
#                'BInt' = benchmark with Intercept,'ALL' (default is 'ALL')
#
# all.sample.fit = Do an estimation using all the data set as training (default is FALSE)
#
# Progression = Give feedback on the estimation. (default 0)
#                 0 = No feedback
#                 1 = Print each iteration
#                 2 = Print after every model estimation & each iteration
#
################################################################

############ Suggestion to use the bypassing option ############
# 
# Run the code with all the available observation in the wanted moving.average unit. 
# Save the following variable :
#   assemblage.info[['x.info']][['x.ranks.ma']] = will the x.rank matrix
#   assemblage.info[['x.info']][['x.comps.ma']] = will the x matrix
#   assemblage.info[['x.info']][['x.weights']] = will the x.weight matrix
#   assemblage.info[['bench.info']][['x.bench.ma']] = will the bench matrix
#   y will be the same matrix
# Then, you can execute a loop by utilizing the 'assemblage' function, employing the 
# recently listed variables, and by setting the 'moving.average' option as 'c()'.
#
################################################################
################################################################
 
  
  ############# See if the packages are installed ###############

  # Get the names of installed packages
  InstalledPKGs <- names(installed.packages()[,'Package'])
  
  # Define the packages you want to install
  myPKGs <- c("glmnet", "pracma", "CVXR", "foreach", "doParallel",'stats','methods','Matrix'
              ,'iterators','datasets','base')
  
  # Identify packages that are not installed
  InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
  
  # If there are packages to install, install them
  if (length(InstallThesePKGs) > 0) {
    install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")
  }

  # load packages
  library(glmnet)
  library(pracma)
  library(CVXR)
  library(foreach)
  library(doParallel)
  library(stats)
  library(methods)
  library(Matrix)
  library(iterators)
  library(datasets)
  library(base)
  
# --- Set the cores for glmnet
ncores=cores

# --- See if the target is in Month over Month units
meanY= mean(y,na.rm = TRUE)
meanX = mean(rowMeans(x,na.rm = TRUE))
if(meanY >2*meanX | meanY <.5*meanX){
  warning("It seems that the target unit is not month over month. This could affect the model performance. \n")
}


if(sum(model.select=='ALL')==1){
  if(isempty(bench)) {
    model.select=c('blend','rank','comp')
  }else{
    model.select=c('blend','rank','comp','BInt','Bench')
  }
}else{
  # Define the desired order
  desired_order <- c('blend', 'rank', 'comp', 'BInt', 'Bench')

  # Reorder model.select based on the desired order
  model.select <- factor(model.select, levels = desired_order)

  # Reorder model.select according to the desired order
  model.select <- as.character(sort(model.select))
  
  if(isempty(bench)) {
    model.select=model.select[!model.select%in%c('BInt','Bench')]
  }
}



# --- Find where to start the pseudo out of sample
if(train.size <= 1){
  end.in.sample = floor(nrow(x)*train.size)
  
  if(train.size == 1){
    all.sample.fit = T
  }
  
}else{
  if(train.size>=nrow(x)){
    print(paste0('The train.size need to be smaller or equal than nrow(x), new train.size=nrow(x) \n'))
    train.size = nrow(x) 
    all.sample.fit = T
  }
  end.in.sample = train.size
  
}

if(!isempty(x.rank) & isempty(moving.average)){ # --- Bypassing the transformation
  # --- Insert the inputs matrix in x.info
  x.info = list('x.ranks.ma' = x.rank, 'x.comps.ma' = x, 'x.weights' = x.weight)

  # --- Are there benchmarks
  if(!is.null(bench)){ # --- YES
    # --- Prepare the benchmark matrix
    bench.info=list('x.bench.ma'= bench)
    
    # --- Save the x info and Bench info
    assemblage.info=list('x.info'=x.info, 'bench.info'=bench.info)
    
  }else{ # --- NO
    # --- Save the x info
    assemblage.info=list('x.info'=x.info[])
  }
  
}else if(isempty(x.rank) & !isempty(moving.average)){ # --- Doing the transformation
  
  # --- Prepare the component and rank matrix
  x.info = x.transformation(x,moving.average)
  
  # --- Take the ranks and components matrix and put it in MA units
  x.info[['x.ranks.ma']] = x.info[['x.ranks']]/moving.average
  x.info[['x.comps.ma']] = x.info[['x.comps']]/moving.average
  
  # --- Calculate the weights 
  x.info[['x.weights']] = weight.transformation(x.weight,x.info[['x.comps.ma']],moving.average)  
  
  # --- Are there benchmarks
  if(!is.null(bench)){ # --- YES
    # --- Prepare the benchmark matrix
    bench.info = bench.transformation(bench,moving.average)
    bench.info[['x.bench.ma']] = bench.info[['x.bench']]/moving.average
    
    # --- Save the x info and Bench info
    assemblage.info=list('x.info'=x.info[], 'bench.info'=bench.info[])
    
  }else{ # --- NO
    # --- Save the x info
    assemblage.info=list('x.info'=x.info[])
  }

}else if(!isempty(x.rank) & !isempty(moving.average)){ # --- Does not respect the conditions
  print(paste0('If both the x.rank and moving.average options are not empty, the function will not work as intended. For transformation to take place, x.rank should be empty while moving.average should not be empty. Conversely, if the input matrices are already in moving average format, x.rank should not be empty, and moving.average should be empty. \n'))
  return()
  
}else if(isempty(x.rank) & isempty(moving.average)){ # --- Does not respect the conditions
  print(paste0('If both the x.rank and moving.average options are empty, the function will not work as intended. For transformation to take place, x.rank should be empty while moving.average should not be empty. Conversely, if the input matrices are already in moving average format, x.rank should not be empty, and moving.average should be empty. \n'))
  return()
    
}

if(!is.null(rank.OOS) | !is.null(comp.OOS)){
  #Do an estimation using all the data set as training for the OOS forecast
  all.sample.fit = T
}

# --- Run the expanding window estimation
assemblage.info = assemblage.estimation.RW( y, x, bench, assemblage.info, moving.average, end.in.sample, window.size, step.estim.window, standardize.values, horizon.gap, coef.sum.1, lambda.grid.C, ncores, model.select, all.sample.fit, Progression)

# --- Calculate RMSE for the POOS section
if ( "Pseudo.OOS.info" %in% names(assemblage.info)) {
  
  assemblage.info[["Pseudo.OOS.info"]] = c(assemblage.info[["Pseudo.OOS.info"]], list('RMSE.POOS' = RMSE.function(assemblage.info[["Pseudo.OOS.info"]][["Prediction"]],assemblage.info[["Pseudo.OOS.info"]][["Start.POOS"]]) ))
  
  # --- Are there benchmarks
  if(sum(model.select=='BInt')>0){ # --- YES
    assemblage.info[["Pseudo.OOS.info"]] = c(assemblage.info[["Pseudo.OOS.info"]], list('RMSE.POOS.norm' = RMSE.function(assemblage.info[["Pseudo.OOS.info"]][["Prediction"]],assemblage.info[["Pseudo.OOS.info"]][["Start.POOS"]],bench) ))
    }
  
}
  
 # --- Calculate RMSE for all the sample
if(all.sample.fit){
  assemblage.info[["All.in.sample.info"]] = c(assemblage.info[["All.in.sample.info"]], list('RMSE.in.sample' = RMSE.function(assemblage.info[["All.in.sample.info"]][["Prediction"]])))
  
  # --- Are there benchmarks
  if(sum(model.select=='BInt')>0){ # --- YES
    assemblage.info[["All.in.sample.info"]] = c(assemblage.info[["All.in.sample.info"]], list('RMSE.in.sample.norm' = RMSE.function(assemblage.info[["All.in.sample.info"]][["Prediction"]],c(),bench)))
  }
}

  # --- Are there still prediction to do
  if(!is.null(rank.OOS) | !is.null(comp.OOS)){
    assemblage.info = prediction.OOS(assemblage.info, rank.OOS, comp.OOS, Bench.OOS, pred.Y.OOS, model.select)
  }

if(!isempty(x.rank) & isempty(moving.average)){ # --- Bypassing the transformation
  cat(sprintf("No errors occurred, but the function utilized the input as the transformed data due \n to its utilization of the x.rank option alongside the moving.average=c() option.  \n"))
}
  return(assemblage.info)
} 




#' @title x.transformation
#' @description Function used to rank and transform the x matrix into the right moving average unit
#' @param x The matrix we need to transform and rank (Should be in month-over-month format)
#' @param moving.average This parameter governs the transformation of the x matrix. If set to 1, the x matrix will be retained for Month over Month analysis. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 signifies Year over Year analysis.
#' @rdname x.transformation
#' @export 
x.transformation = function(x , moving.average){
  
  # Function used to rank and transform the x matrix into the right moving average unit
  
  ################### Legend #####################
  #  
  # x = The matrix we need to transform and rank (Should be in month-over-month format)
  # moving.average = This parameter governs the transformation of the x matrix. 
  #                  If set to 1, the x matrix will be retained for Month over Month analysis. 
  #                  Setting it to 3 corresponds to Quarter over Quarter analysis, 
  #                  while 12 signifies Year over Year analysis.
  #
  ################################################ 
  
  # --- Assign colnames and rownames if there is none
  if(is.null(rownames(x))){
    rownames(x) = paste('t.',1:nrow(x),sep='')
  }
  if(is.null(colnames(x))){
    colnames(x) = paste('V.',1:nrow(x),sep='')
  }
  
  # --- Create output variable
  x.info=list()
  
  # --- Save the x in ma 1
  x.ma = x
  x.info[['x.raw']]=x.ma
  
  # --- Find the rank order 
  x.order=x.ma
  x.rank.ma=x.ma
  for (i in 1:nrow(x.order)) {
    x.order[i,] = order(x.ma[i,]) # Find the order
    x.rank.ma[i,] = x.ma[i,x.order[i,]] # Sort the components
  }
  
  # --- Save the rank order 
  x.info[['x.ma.order']]=x.order
  
  # --- Put the components and the ranks in the right ma
  if(moving.average>1){
    # --- Create the matrix
    x.ranks=array(NA,c(dim(x.rank.ma)[1],dim(x.rank.ma)[2]))
    x.comps=array(NA,c(dim(x.rank.ma)[1],dim(x.rank.ma)[2]))
    # --- Do the loop
    for(i in 1:(dim(x.rank.ma)[1]-moving.average+1)){
      x.ranks[(i+moving.average-1),]=(apply(((as.matrix(x.rank.ma[i:(i+moving.average-1),])/100)+1),2, function(x) prod(x, na.rm = TRUE))-1)*100
      x.comps[(i+moving.average-1),]=(apply(((as.matrix(x.ma[i:(i+moving.average-1),])/100)+1),2, function(x) prod(x, na.rm = TRUE))-1)*100
    }
    
  }else{ # If moving.average == 1, x is already in the right ma
    x.ranks=x.rank.ma
    x.comps=x.ma
  }
  
  # --- Be sure the new matrix have names
  rownames(x.comps) = rownames(x.ma)
  colnames(x.comps) = colnames(x.ma)
  
  rownames(x.ranks) = rownames(x.ma)
  colnames(x.ranks) = colnames(x.ma)
  for (i in 1:ncol(x.ranks)){
    if (i <= 9){
      colnames(x.ranks)[i]=paste("R.00",i,sep="")
    }
    if (i > 9){
      colnames(x.ranks)[i]=paste("R.0",i,sep="")
    }
    if (i > 99){
      colnames(x.ranks)[i]=paste("R.",i,sep="")
    }
  }
  
  # --- Save the matrix
  x.info[['x.ranks']]=x.ranks
  x.info[['x.comps']]=x.comps
  
  return(x.info)
  
}

#' @title bench.transformation
#' @description Function used to transform the bench matrix into the right moving average unit
#' @param bench The benchmark used for model comparison. It can be left empty. If provided, it should follow the month-over-month format and will undergo a transformation similar to x.
#' @param moving.average This parameter governs the transformation of the x matrix. If set to 1, the x matrix will be retained for Month over Month analysis. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 signifies Year over Year analysis.
#' @rdname bench.transformation
#' @export 
bench.transformation = function( bench, moving.average){
  
  # Function used to transform the bench matrix into the right moving average unit
  
  ################### Legend #####################
  #  
  # bench = The benchmark used for model comparison. It can be left empty. 
  #         If provided, it should follow the month-over-month format and will undergo a transformation similar to x.
  # moving.average = This parameter governs the transformation of the x matrix.
  #                  If set to 1, the x matrix will be retained for Month over Month analysis. 
  #                  Setting it to 3 corresponds to Quarter over Quarter analysis, 
  #                  while 12 signifies Year over Year analysis.
  #
  ################################################ 
  
  # --- Assign colnames and rownames if there is none
  if(is.null(rownames(bench))){
    rownames(bench) = paste('t.',1:nrow(bench),sep='')
  }
  if(is.null(colnames(bench))){
    colnames(bench) = paste('B.',1:nrow(bench),sep='')
  }
  
  # --- Create the output variable
  bench.info=list()
  
  # --- Save the bench in ma 1
  bench.ma = bench
  bench.info[['bench.raw']]=bench.ma
  
  # --- Put the benchmarks in the right ma
  if(moving.average>1){
    # --- Create the matrix
    x.bench=array(NA,c(dim(bench.ma)[1],dim(bench.ma)[2]))
    # --- Do the loop
    for(i in 1:(dim(bench.ma)[1]-moving.average+1)){
      x.bench[(i+moving.average-1),]=(apply(((as.matrix(bench.ma[i:(i+moving.average-1),])/100)+1),2, function(x) prod(x, na.rm = TRUE))-1)*100
    }
    
  }else{ # If moving.average == 1, Bench is already in the right ma
    x.bench=bench.ma
  }
  
  # --- Be sure the new matrix have names
  rownames(x.bench) = rownames(bench.ma)
  colnames(x.bench) = colnames(bench.ma)
  
  # --- Save the matrix
  bench.info[['x.bench']]=x.bench
  
  return(bench.info)
  
}


#' @title weight.transformation
#' @description Function used to find and calculate the components weights in the right moving average unit
#' @param x.weight The raw weight matrix can be left empty, or it can be a vector or a matrix. If used, its components should align in the same order as those in 'x'.
#' @param x.comps The components matrix
#' @param moving.average This parameter governs the transformation of the x matrix. If set to 1, the x matrix will be retained for Month over Month analysis. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 signifies Year over Year analysis.
#' @rdname weight.transformation
#' @export 
weight.transformation = function( x.weight, x.comps, moving.average){
  
  # Function used to find and calculate the components weights in the right moving average unit
  
  ################### Legend #####################
  #  
  # x.weight = The raw weight matrix can be left empty, or it can be a vector or a matrix. 
  #            If used, its components should align in the same order as those in 'x'.
  # x.comps = The components matrix
  # moving.average = This parameter governs the transformation of the x matrix. 
  #                  If set to 1, the x matrix will be retained for Month over Month analysis. 
  #                  Setting it to 3 corresponds to Quarter over Quarter analysis, 
  #                  while 12 signifies Year over Year analysis.
  #
  ################################################ 
  
  if(is.null(x.weight)){
    # --- x.weight is empty, we will suppose the components all have the same weight
    
    x.Weight =array((1/ncol(x.comps)), dim(x.comps))
    
  } else if(is.null(dim(x.weight))){
    # --- x.weight is a vector, we will suppose that the weights are time invariant
    
    if (ncol(x.comps)!=length(x.weight)){
      print('The x.weight vector must be the same length as the number of components')
      return(break)
      
    }else{
      # ---  Replicate the vector and do (x.weight/sum(x.weight)) to be sure they sum to 1
      x.Weight= matrix(rep((x.weight/sum(x.weight)), each = nrow(x.comps)), ncol = length(x.weight), byrow = FALSE)
    }
    
  }else{
    # --- x.weight is a matrix, the weights are time variant
    
    if (ncol(x.comps)!=ncol(x.weight) | nrow(x.comps)!=nrow(x.weight)){
      print('The x.weight matrix must be the same size as the x matrix')
      return(break)
      
    }else{
      
      x.Weight=x.weight
      
      # --- fill the NAs
      x.Weight[is.na(x.Weight)]= 1/ncol(x.weight)
      
      # --- Calculate the mean value during the moving.average
      if(moving.average>1) {
        for(i in 1:(dim(x.weight)[1]-moving.average+1)){
          x.Weight[(i+moving.average-1),]=apply(as.matrix(x.Weight[i:(i+moving.average-1),]),2, function(x) mean(x, na.rm = TRUE))
        }
      }
      
      
      # --- Be sure the weights sum to 1
      for (i in (moving.average):nrow(x.Weight)) {
        x.Weight[i,]=x.Weight[i,]/sum(x.Weight[i,],na.rm = TRUE)
      }
      
      # --- Put NAs at the rows used for the higher ma
      if(moving.average>1){
        x.Weight[1:(moving.average-1),] = NA 
      }
    }
  }
  
  # --- Be sure the new matrix have names
  rownames(x.Weight) = rownames(x.comps)
  colnames(x.Weight) = colnames(x.comps)
  
  return(x.Weight)
}



#' @title nonneg.ridge
#' @description Function used to find the coefficients in the nonnegative ridge regression
#' @param y.in The target variable
#' @param x.in The features matrix
#' @param standardize.values Standardize the response variables 
#' @rdname nonneg.ridge
#' @export 
nonneg.ridge = function( y.in, x.in, standardize.values){
  
  # Function used to find the coefficients in the nonnegative ridge regression
  
  ################### Legend #####################
  #  
  # y.in = The target variable
  # x.in = The features matrix
  # standardize.values = Standardize the response variables 
  #
  ################################################   
  
  # --- Create folds for cross-validation
  fd=c(rep(1,nrow(x.in)/10),rep(2,nrow(x.in)/10),rep(3,nrow(x.in)/10),
       rep(4,nrow(x.in)/10),rep(5,nrow(x.in)/10),rep(6,nrow(x.in)/10),
       rep(7,nrow(x.in)/10),rep(8,nrow(x.in)/10),rep(9,nrow(x.in)/10))
  
  fd = c(fd,rep(10,nrow(x.in)-length(fd)))
  
  # --- Cross-Validation
  cvresults = cv.glmnet(x.in, y.in, intercept = F,  alpha=0, lower.limits = 0, 
                        type.measure = "mse",
                        standardize=standardize.values,standardize.response=standardize.values,
                        penalty.factor= apply(x.in,2,sd)^(2*!standardize.values),
                        nfolds = length(unique(fd)),foldid = fd)
  
  # --- Fit the model
  r2 = glmnet(x.in, y.in, intercept = F, lambda=cvresults$lambda.min, alpha=0, lower.limits = 0, 
              type.measure = "mse",penalty.factor=apply(x.in,2,sd)^(2*!standardize.values),
              standardize=standardize.values,standardize.response=standardize.values)
  
  
  return(r2)
  
}


#' @title nonneg.ridge.sum1
#' @description Function used to find the coefficients of the components in the nonnegative ridge regression
#' @param y.in The target variable
#' @param x.in The features matrix
#' @param x.weights The coefficients are adjusted or shrunk closer to the specific values standardize.values = Standardize the response variables 
#' @param standardize.values Standardize the response variables 
#' @param lambda.grid.C Personalize lambda grid
#' @param ncores Number of cores to run the code, helps for glmnet function. (default is 1)
#' @rdname nonneg.ridge.sum1
#' @export 
nonneg.ridge.sum1 = function( y.in, x.in, x.weights, standardize.values, lambda.grid.C=c(),ncores=1){
  
  # Function used to find the coefficients of the components in the 
  # nonnegative ridge regression
  
  ################### Legend #####################
  #  
  # y.in = The target variable
  # x.in = The features matrix
  # x.weights = The coefficients are adjusted or shrunk closer to the specific values
  # standardize.values = Standardize the response variables 
  # lambda.grid.C = Personalize lambda grid
  # ncores = Number of cores to run the code, helps for glmnet function. (default is 1)
  #
  ################################################    
  
  # --- Create folds for cross-validation
  fd=c(rep(1,nrow(x.in)/10),rep(2,nrow(x.in)/10),rep(3,nrow(x.in)/10),
       rep(4,nrow(x.in)/10),rep(5,nrow(x.in)/10),rep(6,nrow(x.in)/10),
       rep(7,nrow(x.in)/10),rep(8,nrow(x.in)/10),rep(9,nrow(x.in)/10))
  
  fd = c(fd,rep(10,nrow(x.in)-length(fd)))
  
  # --- 
  shrinkw = apply(x.in,2,sd)
  # --- Lambda grid
  if(is.null(lambda.grid.C)){  # --- No grid given
    
    # --- Cross-Validation 
    cvresults = cv.glmnet(x.in, y.in, intercept = F,  alpha=0, lower.limits = 0, 
                          type.measure = "mse",
                          standardize=F,standardize.response=standardize.values,
                          penalty.factor= apply(x.in,2,sd),
                          nfolds = length(unique(fd)),foldid = fd)
    lam.seq= c(0,exp(c(0,0.5,seq(1,9,length.out=17))))*cvresults$lambda.1se
    
  }else{ # --- Use grid given
    
    lam.seq = lambda.grid.C
  }
  
  # --- With Cross-Validation
  mse.stack = matrix(NA,length(lam.seq),10)
  if(!length(lam.seq)==1){
    # --- Do it in parallel to speed up the process
    cl <- makeCluster(min(ncores,length(lam.seq)))
    registerDoParallel(cl)
    
    mse.stack <- foreach(lam = 1:length(lam.seq), .combine = 'cbind', .packages =c('CVXR','pracma','glmnet')) %dopar% {
      for(ff in 1:10) {
        # --- Sample Period:
        d.in <- c(which(fd != ff))
        # --- Hold Out Period:
        d.out <- c(which(fd == ff))
        # --- Initialize the Coefficients
        coeffs <- Variable(ncol(x.in))
        # --- Define the Loss-Function
        loss = Minimize(sum((y.in[d.in]-x.in[d.in,]%*%coeffs)^2)+lam.seq[lam]*sum(shrinkw*(coeffs-x.weights)^2))
        # --- Set the constraints
        constr = list(coeffs >=0,sum(coeffs)==1)
        # --- Set the Problem
        prob <- Problem(loss, constr)
        # --- Solve the Problem
        sol <- solve(prob)
        # --- Get the betas
        beta <- sol$getValue(coeffs)
        # --- MSE on the hold out set
        mse.stack[lam, ff] <- as.numeric(mean((y.in[d.out] - x.in[d.out,] %*% beta)^2))
      }
      return(as.numeric(mse.stack[lam,]))
    }
    # Stop the parallel backend
    stopCluster(cl)
    
    mse.stack=t(mse.stack)
    
    # --- Get the Optimal lambda & run estimation on the whole training set
    lbd = lam.seq[order(apply(mse.stack,1,mean))[1]]
  }else{lbd=lam.seq}
  # --- Initialize the Coefficients
  coeffs = Variable(ncol(x.in))
  # --- Define the Loss-Function
  loss = Minimize(sum((y.in-x.in%*%coeffs)^2)+lbd*sum(shrinkw*(coeffs-x.weights)^2))
  # --- Set the constraints
  constr = list(coeffs >=0,sum(coeffs)==1)
  # --- Set the Problem
  prob = Problem(loss,constr)
  # --- Solve the Problem
  sol = solve(prob)
  # --- Get the betas
  beta = sol$getValue(coeffs)
  
  return(as.numeric(beta))
  
}

#' @title nonneg.ridge.meanD
#' @description Function used to find the coefficients of the ranks in the nonnegative ridge regression
#' @param y.in The target variable
#' @param x.in The features matrix
#' @param standardize.values Standardize the response variables 
#' @param lambda.grid.C Personalize lambda grid
#' @param ncores Number of cores to run the code, helps for glmnet function. (default is 1)
#' @rdname nonneg.ridge.meanD
#' @export 
nonneg.ridge.meanD = function( y.in, x.in, standardize.values, lambda.grid.C=c(),ncores=1){
  
  # Function used to find the coefficients of the ranks in the 
  # nonnegative ridge regression
  
  ################### Legend #####################
  #  
  # y.in = The target variable
  # x.in = The features matrix
  # standardize.values = Standardize the response variables 
  # lambda.grid.C = Personalize lambda grid
  # ncores = Number of cores to run the code, helps for glmnet function. (default is 1)
  #
  ################################################    
  library(glmnet)
  library(pracma)
  library(CVXR)
  library(foreach)
  library(doParallel)
  library(stats)
  library(methods)
  library(Matrix)
  library(iterators)
  library(datasets)
  library(base)
  # --- Create folds for cross-validation
  fd=c(rep(1,nrow(x.in)/10),rep(2,nrow(x.in)/10),rep(3,nrow(x.in)/10),
       rep(4,nrow(x.in)/10),rep(5,nrow(x.in)/10),rep(6,nrow(x.in)/10),
       rep(7,nrow(x.in)/10),rep(8,nrow(x.in)/10),rep(9,nrow(x.in)/10))
  
  fd = c(fd,rep(10,nrow(x.in)-length(fd)))
  
  # --- 
  shrinkw = apply(x.in,2,sd)
  # --- Lambda grid
  if(is.null(lambda.grid.C)){  # --- No grid given
    lam.seq = exp(seq(-5,12,length.out=20))
  }else{  # --- Use grid given
    lam.seq = lambda.grid.C
  }
  # --- With Cross-Validation
  mse.stack = matrix(NA,length(lam.seq),10)
  
  if(!length(lam.seq)==1){
    # --- Do it in parallel to speed up the process
    cl <- makeCluster(min(ncores,length(lam.seq)))
    registerDoParallel(cl)
    
    mse.stack <- foreach(lam = 1:length(lam.seq), .combine = 'cbind', .packages =c('CVXR','pracma','glmnet')) %dopar% {
      for(ff in 1:10) {
        # --- Sample Period:
        d.in <- c(which(fd != ff))
        # --- Hold Out Period:
        d.out <- c(which(fd == ff))
        # --- Initialize the Coefficients
        coeffs <- Variable(ncol(x.in))
        # --- Define the Loss-Function
        loss = Minimize(sum((y.in[d.in]-x.in[d.in,]%*%coeffs)^2)+lam.seq[lam]*sum(diff(shrinkw*coeffs)^2))
        # --- Set the constraints
        constr = list(coeffs >=0,t(coeffs)%*%apply(x.in,2,mean) ==mean(y.in))
        # --- Set the Problem
        prob <- Problem(loss, constr)
        # --- Solve the Problem
        sol <- solve(prob)
        # --- Get the betas
        beta <- sol$getValue(coeffs)
        # --- MSE on the hold out set
        mse.stack[lam, ff] <- as.numeric(mean((y.in[d.out] - x.in[d.out,] %*% beta)^2))
      }
      return(as.numeric(mse.stack[lam,]))
    }
    # Stop the parallel backend
    stopCluster(cl)
    
    mse.stack=t(mse.stack)
    
    # --- Get the Optimal lambda & run estimation on the whole training set
    lbd = lam.seq[order(apply(mse.stack,1,mean))[1]]
  }else{lbd=lam.seq}
  # --- Initialize the Coefficients
  coeffs = CVXR::Variable(ncol(x.in),1)
  # --- Define the Loss-Function
  loss = CVXR::Minimize(sum((y.in-x.in%*%coeffs)^2)+lbd*sum(diff(shrinkw*coeffs)^2))
  # --- Set the constraints
  constr = list(coeffs >=0,t(coeffs)%*%apply(x.in,2,mean) ==mean(y.in))
  # --- Set the Problem
  prob = CVXR::Problem(loss,constr)
  # --- Solve the Problem
  sol = CVXR::solve(prob,solver="ECOS")
  # --- Get the betas
  beta = sol$getValue(coeffs)
  
  return(as.numeric(beta))
  
}

#' @title nonneg.ridge.mean
#' @description Function used to find the coefficients of the blend in the nonnegative ridge regression
#' @param y.in The target variable
#' @param x.in The features matrix
#' @param standardize.values Standardize the response variables 
#' @param lambda.grid.C Personalize lambda grid
#' @param ncores Number of cores to run the code, helps for glmnet function. (default is 1)
#' @rdname nonneg.ridge.mean
#' @export 
nonneg.ridge.mean = function( y.in, x.in, standardize.values, lambda.grid.C=c(), ncores=1){
  
  # Function used to find the coefficients of the blend in the 
  # nonnegative ridge regression
  
  ################### Legend #####################
  #  
  # y.in = The target variable
  # x.in = The features matrix
  # standardize.values = Standardize the response variables 
  # lambda.grid.C = Personalize lambda grid
  # ncores = Number of cores to run the code, helps for glmnet function. (default is 1)
  #
  ################################################    
  
  # --- Create folds for cross-validation
  fd=c(rep(1,nrow(x.in)/10),rep(2,nrow(x.in)/10),rep(3,nrow(x.in)/10),
       rep(4,nrow(x.in)/10),rep(5,nrow(x.in)/10),rep(6,nrow(x.in)/10),
       rep(7,nrow(x.in)/10),rep(8,nrow(x.in)/10),rep(9,nrow(x.in)/10))
  
  fd = c(fd,rep(10,nrow(x.in)-length(fd)))
  
  # --- 
  shrinkw = apply(x.in,2,sd)
  # --- Lambda grid
  if(is.null(lambda.grid.C)){  # --- No grid given
    lam.seq = exp(seq(-5,12,length.out=20))
  }else{  # --- Use grid given
    lam.seq = lambda.grid.C
  }
  
  # --- With Cross-Validation
  mse.stack = matrix(NA,length(lam.seq),10)
  if(!length(lam.seq)==1){
    # --- Do it in parallel to speed up the process
    cl <- makeCluster(min(ncores,length(lam.seq)))
    registerDoParallel(cl)
    
    mse.stack <- foreach(lam = 1:length(lam.seq), .combine = 'cbind', .packages =c('CVXR','pracma','glmnet')) %dopar% {
      for(ff in 1:10) {
        # --- Sample Period:
        d.in <- c(which(fd != ff))
        # --- Hold Out Period:
        d.out <- c(which(fd == ff))
        # --- Initialize the Coefficients
        coeffs <- Variable(ncol(x.in))
        # --- Define the Loss-Function
        loss = Minimize(sum((y.in[d.in]-x.in[d.in,]%*%coeffs)^2)+lam.seq[lam]*sum(shrinkw*(coeffs)^2))
        # --- Set the constraints
        constr <- list(coeffs >= 0, t(coeffs) %*% apply(x.in, 2, mean) == mean(y.in))
        # --- Set the Problem
        prob <- Problem(loss, constr)
        # --- Solve the Problem
        sol <- solve(prob)
        # --- Get the betas
        beta <- sol$getValue(coeffs)
        # --- MSE on the hold out set
        mse.stack[lam, ff] <- as.numeric(mean((y.in[d.out] - x.in[d.out,] %*% beta)^2))
      }
      return(as.numeric(mse.stack[lam,]))
    }
    # Stop the parallel backend
    stopCluster(cl)
    
    mse.stack=t(mse.stack)
    # --- Get the Optimal lambda & run estimation on the whole training set
    lbd = lam.seq[order(apply(mse.stack,1,mean))[1]]
  }else{lbd=lam.seq}
  # --- Initialize the Coefficients
  coeffs = Variable(ncol(x.in))
  # --- Define the Loss-Function
  loss = Minimize(sum((y.in-x.in%*%coeffs)^2)+lbd*sum(shrinkw*(coeffs)^2))
  # --- Set the constraints
  constr = list(coeffs >=0,t(coeffs)%*%apply(x.in,2,mean) ==mean(y.in))
  # --- Set the Problem
  prob = Problem(loss,constr)
  # --- Solve the Problem
  sol = solve(prob)
  # --- Get the betas
  beta = sol$getValue(coeffs)
  
  return(as.numeric(beta))
  
}

#' @title Time.function
#' @description Function to have a time estimation for the expanding window
#' @param ooo Iteration of the expanding window 
#' @param end.ooo Number of iterations
#' @param start.time When the expanding window started
#' @rdname Time.function
#' @export 
Time.function = function(ooo, end.ooo,start.time ){
  
  # Function to have a time estimation for the expanding window
  
  ################### Legend #####################
  #  
  # ooo = Iteration of the expanding window 
  # end.ooo = How many iteration there are
  # start.time = When the expanding window started
  #
  ################################################    
  
  # --- Print
  if(ooo>1){
    # --- Calculate the elapsed time
    elapsed.time = difftime(Sys.time(), start.time, units = "secs")
    # --- Estimate the remaining time
    remaining.time = (elapsed.time / ((ooo-1)/(end.ooo))) * (1 - ((ooo-1)/(end.ooo)))
    cat(sprintf('Progress: starting iteration %.0f out of %.0f, Elapsed Time: %.2f seconds, Estimation of the Remaining Time: %.2f seconds \n',
                ooo, end.ooo, elapsed.time, remaining.time))
  }else{
    cat(sprintf('Progress: starting iteration %.0f out of %.0f \n',
                ooo, end.ooo))
  }
}

#' @title assemblage.estimation.RW
#' @description Function estimating the rolling window assemblage
#' @param y Target vector 
#' @param x The raw x matrix
#' @param bench The raw Benchmarks matrix
#' @param assemblage.info  List containing the current assemblage information
#' @param moving.average This parameter governs the transformation of the x matrix. If set to 1, the x matrix will be retained for Month over Month analysis. Setting it to 3 corresponds to Quarter over Quarter analysis, while 12 signifies Year over Year analysis. 
#' @param end.in.sample The observation when the in-sample estimation end
#' @param window.size Size of the rolling window
#' @param step.estim.window Estimating the rolling window every ... observations
#' @param standardize.values Standardize the response variables
#' @param horizon.gap The number of observations to be excluded to ensure that there is no overlapping information within the out-of-sample dataset
#' @param coef.sum.1 Does the coefficients need to sum to one
#' @param lambda.grid.C Personalize lambda grid
#' @param ncores Number of cores to run the code, helps for glmnet function. (default is 1)
#' @param model.select Select the model to run
#' @param all.sample.fit Do an estimation using all the data set as training (default is FALSE)
#' @param Progression Give feedback on the estimation. (default 0)
#' @rdname assemblage.estimation.RW
#' @export 
assemblage.estimation.RW = function( y, x, bench, assemblage.info, moving.average, 
                                     end.in.sample, window.size, step.estim.window, standardize.values, 
                                     horizon.gap, coef.sum.1, lambda.grid.C=c(), ncores, model.select, all.sample.fit, Progression) {
  
  # Function estimating the rolling window
  
  ################### Legend #####################
  #  
  # y = Target vector 
  # x = The raw x matrix
  # bench = The raw Benchmarks matrix
  # assemblage.info = List containing the current assemblage information
  # moving.average = This parameter governs the transformation of the x matrix. 
  #                  If set to 1, the x matrix will be retained for Month over Month analysis. 
  #                  Setting it to 3 corresponds to Quarter over Quarter analysis, 
  #                  while 12 signifies Year over Year analysis.
  # end.in.sample = The observation when the in-sample estimation end
  # window.size = Size of the rolling window
  # step.estim.window = Estimating the rolling window every ... observations
  # standardize.values = Standardize the response variables
  # horizon.gap = The number of observations to be excluded to ensure that there is no 
  #           overlapping information within the out-of-sample dataset
  # coef.sum.1 = Does the coefficients need to sum to one
  # lambda.grid.C = Personalize lambda grid
  # ncores = Number of cores to run the code, helps for glmnet function. (default is 1)
  # model.select = Select the model to run
  # all.sample.fit = Do an estimation using all the data set as training (default is FALSE)
  # Progression = Give feedback on the estimation. (default 0)
  #
  ################################################  
  
  # --- Assign the feature matrix
  if(sum(model.select=='rank')>0){
    x.ranks =assemblage.info[['x.info']][['x.ranks.ma']]
  }
  if(sum(model.select=='comp')>0){
    x.comps =assemblage.info[['x.info']][['x.comps.ma']]
    x.weights =assemblage.info[['x.info']][['x.weights']]
  }
  
  # --- Are there benchmarks
  if(!is.null(bench)) { # --- YES 
    # --- Assign the benchmarks
    if(sum(model.select%in%c('Bench','BInt'))>0){
      x.bench =assemblage.info[['bench.info']][['x.bench.ma']]
    }
    
    # --- Assign the feature for the blend model and save the matrix
    if(sum(model.select=='blend')>0){
      x.bench =assemblage.info[['bench.info']][['x.bench.ma']]
      x.blend = cbind(x.bench,x.comps,x.ranks)
      assemblage.info[['x.info']][['x.blend.ma']] =x.blend
    }
    
    
    # --- Create the prediction matrix
    fitted.values.all.in = array(NA,c(nrow(x),6))  
    colnames(fitted.values.all.in)= c('Target', 'Mod.blend', 'Mod.rank', 'Mod.comp','Mod.BInt','Mod.Bench')
    
  }else{ # --- NO
    # --- Assign the feature for the blend model and save the matrix
    if(sum(model.select=='blend')>0){
      x.blend = cbind(x.comps,x.ranks)
      assemblage.info[['x.info']][['x.blend.ma']] =x.blend
    }
    
    # --- Create the prediction matrix
    fitted.values.all.in = array(NA,c(nrow(x.comps),4))  
    colnames(fitted.values.all.in)= c('Target', 'Mod.blend', 'Mod.rank', 'Mod.comp')
  }
  fitted.values.all.in[1:length(y),1]= y
  rownames(fitted.values.all.in)=rownames(x)
  
  
  # --- Find how many iteration we are doing
  if(end.in.sample==nrow(x) ){
    horizon = 0 # --- Fit on all the available data
    end.ooo = 0 # --- Only doing the in-sample estimation
    all.sample.fit = TRUE
  }else{
    end.ooo = ifelse((nrow(x)-end.in.sample)<step.estim.window, 1, floor( (nrow(x)-end.in.sample)/step.estim.window )) 
  }
  
  if(is.null(window.size)){ # --- Use all the data
    window.size = end.in.sample 
  }
  
  # --- Are we doing the pseudo out of sample expanding window
  if(end.ooo>0){ # --- YES
    
    # --- Create the prediction matrix for the POOS and coefficients
    poos.prediction = fitted.values.all.in
    
    if(sum(model.select=='blend')>0){
      coefficients.blend =array(NA,c(end.ooo,(ncol(x.blend))))
      rownames(coefficients.blend)=paste(1:end.ooo)
      colnames(coefficients.blend)=c(colnames(x.blend))
    }
    if(sum(model.select=='rank')>0){
      coefficients.rank =array(NA,c(end.ooo,(ncol(x))))
      rownames(coefficients.rank)=paste(1:end.ooo)
      colnames(coefficients.rank)=c(colnames(x.ranks))
    }
    
    if(sum(model.select=='comp')>0){
      coefficients.comp =array(NA,c(end.ooo,(ncol(x))))
      rownames(coefficients.comp)=paste(1:end.ooo)
      colnames(coefficients.comp)=c(colnames(x.comps))
    }
    
    # --- Are there benchmarks
    if(!is.null(bench)) { # --- YES 
      if(sum(model.select=='BInt')>0){
        coefficients.BInt =array(NA,c(end.ooo,(ncol(x.bench)+1)))
        rownames(coefficients.BInt)=paste(1:end.ooo)
        colnames(coefficients.BInt)=c('intercept',colnames(x.bench))
      }
      
      if(sum(model.select=='Bench')>0){
        coefficients.Bench =array(NA,c(end.ooo,(ncol(x.bench))))
        rownames(coefficients.Bench)=paste(1:end.ooo)
        colnames(coefficients.Bench)=c(colnames(x.bench))
      }
    }
    
    # --- Give a feedback on the time it will take
    if(Progression>0){
      start.time = Sys.time()
      cat(sprintf("Starting assemblage rolling estimation \n"))
    }
    
    # --- Start the iteration loop
    for (ooo in 1:end.ooo){
      
      # --- Give a feedback on the time it will take
      if(Progression>0){
        Time.function(ooo,end.ooo,start.time)
      }
      
      end.in = end.in.sample - horizon.gap + (ooo-1) * step.estim.window   # --- Current in-sample end
      start.in = max(moving.average, end.in - window.size) # --- Take into account the observation used for the ma
      start.oo = min(nrow(x), end.in.sample + 1 + (ooo-1) * step.estim.window)   # --- Current POOS start
      end.oo = min(nrow(x), end.in.sample  + (ooo) * step.estim.window )  # --- Current POOS end
      
      y.in = y[start.in:end.in] # --- Training target
      
      
      # ------------------------------------ Model Blend ------------------------------------------ #
      
      if(sum(model.select=='blend')>0){
        if(Progression>1){
          cat(sprintf("Estimating mod.blend \n"))
        }
        x.in = cbind(x.blend)[start.in:end.in,] # --- Training feature
        betas = nonneg.ridge.mean(y.in,x.in,standardize.values,lambda.grid.C,ncores) # --- Estimation
        
        # --- Saving the coefficients
        coefficients.blend[ooo,]=betas
        rownames(coefficients.blend)[ooo] =rownames(x)[end.in]
        # --- Saving the predictions
        if(ooo == 1){
          poos.prediction[start.in:(start.oo-1),2] = cbind(x.blend)[start.in:(start.oo-1),] %*% betas
        }
        poos.prediction[start.oo:end.oo,2] = cbind(x.blend)[start.oo:end.oo,] %*% betas
      }
      
      # ------------------------------------ Model Rank ------------------------------------------ #
      
      if(sum(model.select=='rank')>0){
        if(Progression>1){
          cat(sprintf("Estimating mod.rank \n"))
        }
        x.in = cbind(x.ranks)[start.in:end.in,] # --- Training feature
        betas = nonneg.ridge.meanD(y.in,x.in,standardize.values,lambda.grid.C,ncores) # --- Estimation
        
        # --- Saving the coefficients
        coefficients.rank[ooo,]=betas
        rownames(coefficients.rank)[ooo] =rownames(x)[end.in]
        # --- Saving the predictions
        if(ooo == 1){
          poos.prediction[start.in:(start.oo-1),3] = cbind(x.ranks)[start.in:(start.oo-1),] %*% betas
        }
        poos.prediction[start.oo:end.oo,3] = cbind(x.ranks)[start.oo:end.oo,] %*% betas
      }
      
      # ------------------------------------ Model Component ------------------------------------------ #
      
      if(sum(model.select=='comp')>0){
        if(Progression>1){
          cat(sprintf("Estimating mod.comp \n"))
        }
        x.in = cbind(x.comps)[start.in:end.in,] # --- Training feature
        weights = x.weights[end.in,] # --- Training weights
        
        # --- Should the weights sum to one
        if(coef.sum.1) { # --- YES
          betas = nonneg.ridge.sum1( y.in, x.in, weights, standardize.values,lambda.grid.C,ncores) # --- Estimation
          
          # --- Saving the coefficients
          coefficients.comp[ooo,]= betas
          rownames(coefficients.comp)[ooo] =rownames(x)[end.in]
          # --- Saving the predictions
          if(ooo == 1){
            poos.prediction[start.in:(start.oo-1),4] = as.vector(x.comps[start.in:(start.oo-1),] %*% betas)
          }
          poos.prediction[start.oo:end.oo,4] = as.vector(x.comps[start.oo:end.oo,] %*% betas)
          
        }else{ # --- NO
          r2 = nonneg.ridge(y.in,x.in,standardize.values) # --- Estimation
          
          # --- Saving the coefficients
          coefficients.comp[ooo,]=as.vector(r2$beta)
          rownames(coefficients.comp)[ooo] =rownames(x)[end.in]
          # --- Saving the predictions
          if(ooo == 1){
            poos.prediction[start.in:(start.oo-1),4] = predict(r2,newx=cbind(x.comps)[start.in:(start.oo-1),])
          }
          poos.prediction[start.oo:end.oo,4] = predict(r2,newx=cbind(x.comps)[start.oo:end.oo,])
        }
      }
      
      # --- Are there benchmarks
      if(!is.null(bench)){ # --- YES
        
        # ------------------------------------ Model Bench (WITH Intercept) ------------------------------------------ #
        
        if(sum(model.select=='BInt')>0){
          if(Progression>1){
            cat(sprintf("Estimating mod.BInt \n"))
          }
          x.in = cbind(x.bench)[start.in:end.in,] # --- Training feature
          # --- Estimation
          r2 = glmnet(x.in, y.in, intercept = T, lambda=0, alpha=0, lower.limits = 0, 
                      type.measure = "mse",
                      standardize=T,standardize.response=F)
          
          # --- Saving the coefficients
          coefficients.BInt[ooo,]=c(as.vector(r2$a0),as.vector(r2$beta))
          rownames(coefficients.BInt)[ooo] =rownames(x)[end.in]
          # --- Saving the predictions
          if(ooo == 1){
            poos.prediction[start.in:(start.oo-1),5] = predict(r2,newx=cbind(x.bench)[start.in:(start.oo-1),])
          }
          poos.prediction[start.oo:end.oo,5] = predict(r2,newx=cbind(x.bench)[start.oo:end.oo,])
          
        }
        
        # ------------------------------------ Model Bench (NO Intercept) ------------------------------------------ #
        
        if(sum(model.select=='Bench')>0){
          if(Progression>1){
            cat(sprintf("Estimating mod.Bench \n"))
          }
          x.in = cbind(x.bench)[start.in:end.in,] # --- Training feature
          
          # --- Should the weights sum to one
          if(coef.sum.1) {# --- YES
            betas = nonneg.ridge.sum1( y.in, x.in, rep(1/ncol(x.bench),times=ncol(x.bench)), standardize.values,lambda.grid.C,ncores) # --- Estimation
            
            # --- Saving the coefficients
            coefficients.Bench[ooo,]= betas
            rownames(coefficients.Bench)[ooo] =rownames(x)[end.in]
            # --- Saving the predictions
            if(ooo == 1){
              poos.prediction[start.in:(start.oo-1),6] = as.vector(x.bench[start.in:(start.oo-1),] %*% betas)
            }
            poos.prediction[start.oo:end.oo,6] = as.vector(x.bench[start.oo:end.oo,] %*% betas)
            
          }else{# --- NO
            r2 = nonneg.ridge(y.in,x.in,standardize.values) # --- Estimation
            
            # --- Saving the coefficients
            coefficients.Bench[ooo,]=t(as.matrix(r2$beta))
            rownames(coefficients.Bench)[ooo] =rownames(x)[end.in]
            # --- Saving the predictions
            if(ooo == 1){
              poos.prediction[start.in:(start.oo-1),6] = predict(r2,newx=cbind(x.bench)[start.in:(start.oo-1),])
            }
            poos.prediction[start.oo:end.oo,6] = predict(r2,newx=cbind(x.bench)[start.oo:end.oo,])
          }
        }
        
      }
      
    }
    
    
    # --- Saving the coefficients inside assemblage.info
    Coefficients =list()
    
    if (sum(model.select == 'blend') > 0) {
      Coefficients[['Mod.blend']] <- coefficients.blend
    }
    
    if (sum(model.select == 'rank') > 0) {
      Coefficients[['Mod.rank']] <- coefficients.rank
    }
    
    if (sum(model.select == 'comp') > 0) {
      Coefficients[['Mod.comp']] <- coefficients.comp
    }
    
    if (sum(model.select == 'BInt') > 0) {
      Coefficients[['Mod.BInt']] <- coefficients.BInt
    }
    
    if (sum(model.select == 'Bench') > 0) {
      Coefficients[['Mod.Bench']] <- coefficients.Bench
    }
    
    poos.prediction=poos.prediction[,c(TRUE,c('blend','rank','comp','BInt','Bench')%in%model.select)]
    
    assemblage.info[['Pseudo.OOS.info']]=list('Prediction'=poos.prediction,'Coefficients'=Coefficients,
                                              'Start.POOS'= rownames(x)[min(nrow(x), end.in.sample + 1)] )
    
  }
  
  if(all.sample.fit==TRUE){ # --- Do an estimation using all the observation 
    if(Progression>0){
      cat(sprintf("Starting assemblage all.sample.fit estimation \n"))
    }
    
    start.in = ifelse(is.null(moving.average),1,moving.average)  # --- Take into account the observation used for the ma
    end.in = nrow(x) - horizon.gap  # --- Current in-sample end
    
    y.in = y[start.in:end.in] # --- Training target
    
    
    # ------------------------------------ Model Blend ------------------------------------------ #
    
    if(sum(model.select=='blend')>0){
      if(Progression>1){
        cat(sprintf("Estimating mod.blend \n"))
      }
      x.in = cbind(x.blend)[start.in:end.in,] # --- Training feature
      betas = nonneg.ridge.mean(y.in,x.in,standardize.values,lambda.grid.C,ncores) # --- Estimation
      
      # --- Saving the coefficients
      coefficients.blend.all=t(betas)
      rownames(coefficients.blend.all) =rownames(x)[end.in]
      colnames(coefficients.blend.all) =colnames(x.blend)
      # --- Saving the predictions
      fitted.values.all.in[start.in:nrow(x),2] = cbind(x.blend)[start.in:nrow(x),] %*% betas
    }
    
    # ------------------------------------ Model Rank ------------------------------------------ #
    
    if(sum(model.select=='rank')>0){
      if(Progression>1){
        cat(sprintf("Estimating mod.rank \n"))
      }
      x.in = cbind(x.ranks)[start.in:end.in,] # --- Training feature
      betas = nonneg.ridge.meanD(y.in,x.in,standardize.values,lambda.grid.C,ncores) # --- Estimation
      
      # --- Saving the coefficients
      coefficients.rank.all=t(betas)
      rownames(coefficients.rank.all) =rownames(x)[end.in]
      colnames(coefficients.rank.all) =colnames(x.ranks)
      # --- Saving the predictions
      fitted.values.all.in[start.in:nrow(x),3] = cbind(x.ranks)[start.in:nrow(x),] %*% betas
    }
    
    # ------------------------------------ Model Component ------------------------------------------ #
    
    if(sum(model.select=='comp')>0){
      if(Progression>1){
        cat(sprintf("Estimating mod.comp \n"))
      }
      x.in = cbind(x.comps)[start.in:end.in,] # --- Training feature
      weights = x.weights[end.in,] # --- Training weights
      
      # --- Should the weights sum to one
      if(coef.sum.1) { # --- YES
        betas = nonneg.ridge.sum1( y.in, x.in, weights, standardize.values,lambda.grid.C,ncores) # --- Estimation
        
        # --- Saving the coefficients
        coefficients.comp.all=t(betas)
        rownames(coefficients.comp.all) =rownames(x)[end.in]
        colnames(coefficients.comp.all) = colnames(x)
        # --- Saving the predictions
        fitted.values.all.in[start.in:nrow(x),4] = as.vector(x.comps[start.in:nrow(x),] %*% betas)
        
      }else{ # --- NO
        r2 = nonneg.ridge(y.in,x.in,standardize.values) # --- Estimation
        
        # --- Saving the coefficients
        coefficients.comp.all=t(as.matrix(r2$beta))
        rownames(coefficients.comp.all) =rownames(x)[end.in]
        colnames(coefficients.comp.all) = colnames(x)
        # --- Saving the predictions
        fitted.values.all.in[start.in:nrow(x),4] = predict(r2,newx=cbind(x.comps)[start.in:nrow(x),])
      }
    }
    
    # --- Are there benchmarks
    if(!is.null(bench)){ # --- YES
      
      # ------------------------------------ Model Bench (WITH Intercept) ------------------------------------------ #
      
      if(sum(model.select=='BInt')>0){
        if(Progression>1){
          cat(sprintf("Estimating mod.BInt \n"))
        }
        x.in = cbind(x.bench)[start.in:end.in,] # --- Training feature
        # --- Estimation
        r2 = glmnet(x.in, y.in, intercept = T, lambda=0, alpha=0, lower.limits = 0, 
                    type.measure = "mse",
                    standardize=T,standardize.response=F)
        
        # --- Saving the coefficients
        coefficients.BInt.all=t(c(as.vector(r2$a0),as.vector(r2$beta)))
        rownames(coefficients.BInt.all) =rownames(x)[end.in]
        colnames(coefficients.BInt.all) = c('Intercept',colnames(x.bench))
        # --- Saving the predictions
        fitted.values.all.in[start.in:nrow(x),5] = predict(r2,newx=cbind(x.bench)[start.in:nrow(x),])
      }
      
      # ------------------------------------ Model Bench (NO Intercept) ------------------------------------------ #
      
      if(sum(model.select=='Bench')>0){
        if(Progression>1){
          cat(sprintf("Estimating mod.Bench \n"))
        }
        x.in = cbind(x.bench)[start.in:end.in,] # --- Training feature
        
        # --- Should the weights sum to one
        if(coef.sum.1) {# --- YES
          betas = nonneg.ridge.sum1( y.in, x.in, rep(1/ncol(x.bench),times=ncol(x.bench)), standardize.values,lambda.grid.C,ncores) # --- Estimation
          
          # --- Saving the coefficients
          coefficients.Bench.all=t(betas)
          rownames(coefficients.Bench.all) =rownames(x)[end.in]
          colnames(coefficients.Bench.all) = colnames(x.bench)
          # --- Saving the predictions
          fitted.values.all.in[start.in:nrow(x),6] = as.vector(x.bench[start.in:nrow(x),] %*% betas)
          
        }else{# --- NO
          r2 = nonneg.ridge(y.in,x.in,standardize.values) # --- Estimation
          
          # --- Saving the coefficients
          coefficients.Bench.all=t(as.matrix(r2$beta))
          rownames(coefficients.Bench.all) =rownames(x)[end.in]
          colnames(coefficients.Bench.all) = colnames(x.bench)
          # --- Saving the predictions
          fitted.values.all.in[start.in:nrow(x),6] = predict(r2,newx=cbind(x.bench)[start.in:nrow(x),])
        }
      }
      
    } 
    
    
    # --- Saving the coefficients inside assemblage.info
    Coefficients <- list()
    
    if (sum(model.select == 'blend') > 0) {
      Coefficients[['Mod.blend']] <- coefficients.blend.all
    }
    
    if (sum(model.select == 'rank') > 0) {
      Coefficients[['Mod.rank']] <- coefficients.rank.all
    }
    
    if (sum(model.select == 'comp') > 0) {
      Coefficients[['Mod.comp']] <- coefficients.comp.all
    }
    
    if (sum(model.select == 'BInt') > 0) {
      Coefficients[['Mod.BInt']] <- coefficients.BInt.all
    }
    
    if (sum(model.select == 'Bench') > 0) {
      Coefficients[['Mod.Bench']] <- coefficients.Bench.all
    }
    
    if(!is.null(bench)) {
      fitted.values.all.in=fitted.values.all.in[,c(TRUE,c('blend','rank','comp','BInt','Bench')%in%model.select)]
    }else{
      fitted.values.all.in=fitted.values.all.in[,c(TRUE,c('blend','rank','comp')%in%model.select)]
    }
    assemblage.info[['All.in.sample.info']]=list('Prediction'=fitted.values.all.in,'Coefficients'=Coefficients)
  }
  
  return(assemblage.info)
}

#' @title RMSE.function
#' @description Function to calculate the RMSE
#' @param Prediction The prediction matrix
#' @param Start.POOS Where the POOS start
#' @param bench The benchmark matrix
#' @rdname RMSE.function
#' @export 
RMSE.function = function(Prediction,Start.POOS=c(), bench=c()) {
  
  # Function to calculate the RMSE
  
  ################### Legend #####################
  #  
  # Prediction = The prediction matrix
  # Start.POOS = Where the POOS start
  # bench = The benchmark matrix
  #
  ################################################   
  
  # --- Is there a POOS section ?
  if(is.null(Start.POOS)){ # --- NO
    # --- Find the first non NA value
    start.index=which(!is.na(Prediction[,2]))[1]
  }else{ # --- YES
    # --- Find the row index where the POOS started
    start.index =which(rownames(Prediction)==Start.POOS)
  }
  
  # --- Find the dimension
  end.index = tail(which(!is.na(Prediction[,1])),n=1)
  col.index = ncol(Prediction)
  
  # --- Calculate the RMSE
  if(ncol(Prediction)==2){
    RMSE = sqrt(mean((Prediction[start.index:end.index,2:col.index]-Prediction[start.index:end.index,1])^2))
  }else{
    RMSE = t(as.matrix(sqrt(colMeans((Prediction[start.index:end.index,2:col.index]-Prediction[start.index:end.index,1])^2))))
  }
  
  # --- Are there benchmarks
  if(!is.null(bench)){ # --- YES
    # --- Normalize by mod BInt
    
    RMSE = RMSE/RMSE[colnames(Prediction)[2:ncol(Prediction)]=='Mod.BInt']
  }
  
  return(RMSE)
}


#' @title prediction.OOS
#' @description Function to do the OOS prediction using the coefficients that used all the data set for training.
#' @param assemblage.info Information from the previous steps
#' @param rank.OOS Rank matrix
#' @param comp.OOS Component matrix
#' @param Bench.OOS Benchmark matrix
#' @param pred.Y.OOS Target for the OOS forecast
#' @param model.select Select the model
#' @rdname prediction.OOS
#' @export 
prediction.OOS = function(assemblage.info, rank.OOS, comp.OOS, Bench.OOS, pred.Y.OOS, model.select) {
  
  # Function to do the OOS prediction using the coefficients that used all the data set for training.
  
  ################### Legend #####################
  #  
  # assemblage.info = Information from the previous steps
  # rank.OOS = Rank matrix
  # comp.OOS = Component matrix
  # Bench.OOS = Benchmark matrix
  # pred.Y.OOS = Target for the OOS forecast
  # model.select = Select the model
  #
  ################################################   
  
  # --- Are there benchmarks
  if(!is.null(Bench.OOS)){ # --- YES
    if (sum(model.select == 'BInt') > 0) {
      pred.BInt.OOS = assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.BInt"]][1] + Bench.OOS %*% assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.BInt"]][2:(ncol(Bench.OOS)+1)]
    }
    if (sum(model.select == 'Bench') > 0) {
      pred.Bench.OOS = Bench.OOS %*% t(assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.Bench"]])
    }
    if (sum(model.select == 'blend') > 0) {
      pred.blend.OOS = cbind(Bench.OOS,comp.OOS,rank.OOS) %*% t(assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.blend"]])
    }
  }else{ # --- NO
    # --- Be sure the estimation didn't had benchmarks 
    if (sum(model.select == 'blend') > 0) {
      ncolBlend = length(assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.blend"]])
      pred.blend.OOS = cbind(comp.OOS,rank.OOS) %*% assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.blend"]][max(1,(ncolBlend-ncol(cbind(comp.OOS,rank.OOS))+1)):ncolBlend]
    }
  }
  
  if (sum(model.select == 'rank') > 0) {
    pred.rank.OOS = rank.OOS %*% t(assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.rank"]])
  }
  if (sum(model.select == 'comp') > 0) {
    pred.comp.OOS = comp.OOS %*% t(assemblage.info[["All.in.sample.info"]][["Coefficients"]][["Mod.comp"]])
  }
  
  if(!is.null(pred.Y.OOS)){
    # --- Are there benchmarks
    if(!is.null(Bench.OOS)){ # --- YES
      
      # --- Be sure everything is the same length
      sample.size = sample.max = length(pred.Y.OOS)
      for (i in model.select) {
        sample.size=min(sample.size,length(get(paste0('pred.', i, '.OOS',sep=''))))
        sample.max=max(sample.max,length(get(paste0('pred.', i, '.OOS',sep=''))))
      }
      
      # --- Merge the prediction
      prediction = as.matrix(pred.Y.OOS[1:sample.size])
      for (i in model.select) {
        prediction = cbind(prediction,as.matrix(get(paste0('pred.', i, '.OOS',sep=''))[1:sample.size]))
      }
      colnames(prediction) = c('Target','Mod.blend','Mod.rank', 'Mod.comp','Mod.BInt','Mod.Bench')[c(T,c('blend','rank','comp','BInt','Bench')%in%model.select)]
      
      
      RMSE.OOS = as.matrix(RMSE.function(prediction,c(),Bench.OOS) )
      
      # --- Saving the Prediction inside assemblage.info
      Prediction =as.matrix(array(NA,c(sample.max,(1+length(model.select)))))
      Prediction[1:length(pred.Y.OOS),1] <- pred.Y.OOS
      rownames(Prediction)[1:length(pred.Y.OOS)] <- names(pred.Y.OOS)
      ii=1
      for (i in model.select) {
        ii = ii +1
        Prediction[1:length(get(paste0('pred.', i, '.OOS',sep=''))),ii] = get(paste0('pred.', i, '.OOS',sep=''))
      }
      colnames(Prediction) = c('Target','Mod.blend','Mod.rank', 'Mod.comp','Mod.BInt','Mod.Bench')[c(T,c('blend','rank','comp','BInt','Bench')%in%model.select)]
      
      # --- Save results
      assemblage.info[['OOS']]= list('Prediction'=Prediction,'RMSE'=RMSE.OOS)
      
    }else{
      
      # --- Be sure everything is the same length
      sample.size = sample.max = length(pred.Y.OOS)
      for (i in model.select) {
        sample.size=min(sample.size,length(get(paste0('pred.', i, '.OOS',sep=''))))
        sample.max=max(sample.max,length(get(paste0('pred.', i, '.OOS',sep=''))))
      }
      
      # --- Merge the prediction
      prediction = as.matrix(pred.Y.OOS[1:sample.size])
      for (i in model.select) {
        prediction = cbind(prediction,as.matrix(get(paste0('pred.', i, '.OOS',sep=''))[1:sample.size]))
      }
      colnames(prediction) = c('Target','Mod.blend','Mod.rank', 'Mod.comp','Mod.BInt','Mod.Bench')[c(T,c('blend','rank','comp','BInt','Bench')%in%model.select)]
      
      RMSE.OOS = RMSE.function(prediction) 
      
      
      # --- Saving the Prediction inside assemblage.info
      Prediction =as.matrix(array(NA,c(sample.max,(1+length(model.select)))))
      Prediction[1:length(pred.Y.OOS),1] <- pred.Y.OOS
      rownames(Prediction)[1:length(pred.Y.OOS)] <- names(pred.Y.OOS)
      ii=1
      for (i in model.select) {
        ii = ii +1
        Prediction[1:length(get(paste0('pred.', i, '.OOS',sep=''))),ii] = get(paste0('pred.', i, '.OOS',sep=''))
      }
      colnames(Prediction) = c('Target','Mod.blend','Mod.rank', 'Mod.comp','Mod.BInt','Mod.Bench')[c(T,c('blend','rank','comp','BInt','Bench')%in%model.select)]
      
      
      # --- Save results
      assemblage.info[['OOS']]= list('Prediction'=Prediction,'RMSE'=RMSE.OOS)
    }
  }else{
    
    # --- Be sure everything is the same length
    sample.max = 0
    for (i in model.select) {
      sample.max=max(sample.max,length(get(paste0('pred.', i, '.OOS',sep=''))))
    }
    
    # --- Saving the Prediction inside assemblage.info
    Prediction =as.matrix(array(NA,c(sample.max,(length(model.select)))))
    ii=0
    for (i in model.select) {
      ii = ii +1
      Prediction[1:length(get(paste0('pred.', i, '.OOS',sep=''))),ii] = get(paste0('pred.', i, '.OOS',sep=''))
      if(sample.max==length(get(paste0('pred.', i, '.OOS',sep=''))) & i!='BInt') {
        rownames(Prediction)=rownames(get(paste0(i, '.OOS',sep='')))
      }
    }
    colnames(Prediction) = c('Mod.blend','Mod.rank', 'Mod.comp','Mod.BInt','Mod.Bench')[c('blend','rank','comp','BInt','Bench')%in%model.select]
    
    # --- Save results
    assemblage.info[['OOS']]= list('Prediction'=Prediction)
    
  }
  
  return(assemblage.info)
  
}









