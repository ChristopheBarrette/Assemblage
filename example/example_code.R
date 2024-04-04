rm(list = ls())
# install.packages("devtools")
# devtools::install_github("ChristopheBarrette/Assemblage")

# Install the package and load it
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
library(Assemblage)

#####Paths#######
gen_path="~/"  #Specific to the computer
#Relative to gen path
wd_data =paste(gen_path,"data/",sep="")

cores = 1

# Doing forecasting of the next 12 months so need to remove the last 12 observations for the rmse, etc.
horizon.gap = 12

# Doing rolling window of 240
window.size = 240

# End training in 2019-12-01
end.dates = '2019-12-01'

# we do not want to re estimate the model past 2019-12-01
step.estim.window = 100000

# Select the models to run 'comp','rank','blend','BInt','Bench'
model.select = c('comp','rank','BInt','Bench')

# Import the USA data set 
load(paste(wd_data,"USA_PCE_lvl6_C.RData",sep=""))
load(paste(wd_data,"USA_PCE_lvl3_R.RData",sep=""))
load(paste(wd_data,"weight_USA_lvl6.RData",sep=""))

#select the in-sample and out of sample period
insample= (which(rownames(weight_USA_lvl6)==end.dates)-window.size-horizon.gap):which(rownames(weight_USA_lvl6)==end.dates)
outsample= (which(rownames(weight_USA_lvl6)==end.dates)+1) :nrow(USA_PCE_lvl6_C)

# Do the assemblage version with the inputs already transformed
assemblage_info = assemblage(y=USA_PCE_lvl6_C[insample,1], x=USA_PCE_lvl6_C[insample,5:ncol(USA_PCE_lvl6_C)], 
                         x.weight = weight_USA_lvl6[insample,], x.rank = USA_PCE_lvl3_R[insample,],
                         bench = USA_PCE_lvl6_C[insample,2:4], # in sample estimation part
                         
                          pred.Y.OOS = USA_PCE_lvl6_C[outsample,1],comp.OOS=USA_PCE_lvl6_C[outsample,5:ncol(USA_PCE_lvl6_C)], 
                         rank.OOS = USA_PCE_lvl3_R[outsample,], Bench.OOS=USA_PCE_lvl6_C[outsample,2:4], #out of sample part
                         
                         cores=cores, horizon.gap = horizon.gap, moving.average=c(),
                         model.select = model.select, Progression=2) #parameter part

print(assemblage_info[["OOS"]][["RMSE"]])


# Import the Canadian data set 
load(paste(wd_data,"Canada_PCE_lvl4.RData",sep=""))

#select the in-sample
Canada_PCE_lvl4=as.matrix(Canada_PCE_lvl4)
insample= (which(rownames(Canada_PCE_lvl4)==end.dates))

# Do the assemblage version with the transformation
assemblage_info = assemblage(y=Canada_PCE_lvl4[,1], x=Canada_PCE_lvl4[,5:ncol(Canada_PCE_lvl4)], bench = Canada_PCE_lvl4[,2:4],# in sample estimation part
                         
                         train.size = insample, step.estim.window = step.estim.window,horizon.gap = horizon.gap, #parameter part
                         cores=cores, moving.average=3, window.size = window.size, model.select = model.select, Progression=2)


print(assemblage_info[["Pseudo.OOS.info"]][["RMSE.POOS.norm"]])


