#' Optimisation of undersampling percentages for classification
#'
#' @description The function \code{\link{OptPrInDT}} applies an iterative technique for finding optimal undersampling percentages 
#'  'percl' for the larger class and 'percs' for the smaller class by a nested grid search for the use of the function \code{\link{PrInDT}} for
#' the relationship between the two-class factor variable 'classname' and all other factor and numerical variables
#' in the data frame 'data' by means of 'N' repetitions of undersampling. The optimization citerion is the balanced accuracy 
#' on the validation sample 'valdat' (default = full sample 'data'). The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The inputs plmax and psmax determine the maximal values of the percentages and the inputs distl and dists the
#' the distances to the next smaller percentage to be tried.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#' The parameter 'steps' controls, how many of the 3 possible optimization steps should be carried out; default=3.\cr
#'
#' @usage OptPrInDT(data,classname,ctestv=NA,N=99,plmax=0.09,psmax=0.9,
#'                distl=0.01,dists=0.1,conf.level=0.95,minsplit=NA,minbucket=NA,
#'                valdat=data,steps=3)
#'
#' @param data Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param N Number (> 7) of repetitions (integer)
#' @param plmax Maximal undersampling percentage of larger class (numerical, > 0 and <= 1);\cr
#'     default = 0.09
#' @param psmax Maximal undersampling percentage of smaller class (numerical, > 0 and <= 1);\cr
#'     default = 0.9
#' @param distl Distance to the next lower undersampling percentage of larger class (numerical, > 0 and < 1);\cr
#'     default = 0.01
#' @param dists Distance to the next lower undersampling percentage of smaller class (numerical, > 0 and < 1);\cr
#'     default = 0.1
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#' @param valdat validation data; default = data
#' @param steps number of optimization steps = 1, 2, 3; default = 3
#'
#' @return
#' \describe{
#' \item{besttree}{best tree on full sample}
#' \item{bestba}{balanced accuracy of best tree on full sample}
#' \item{percl}{undersampling percentage of large class of best tree on full sample}
#' \item{percs}{undersampling percentage of small class of best tree on full sample}
#' }
#'
#' @details
#' See help("RePrInDT") and help("PrInDT") for further information.
#'
#' Standard output can be produced by means of \code{print(name$besttree)} or just \code{name$besttree} as well as \code{plot(name$besttree)} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export OptPrInDT
#'
#' @examples
#' datastrat <- PrInDT::data_zero
#' data <- na.omit(datastrat) # cleaned full data: no NAs
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}','MLU == {1, 3}') # split exclusions
#' # call of OptPrInDT
#' out <- OptPrInDT(data,"real",ctestv,N=24,conf.level=0.99,steps=1) # unstratified
#' out # print best model and ensembles as well as all observations
#' plot(out)
#'
OptPrInDT <- function(data,classname,ctestv=NA,N=99,plmax=0.09,psmax=0.9,distl=0.01,dists=0.1,conf.level=0.95,minsplit=NA,minbucket=NA,valdat=data,steps=3){
  ## input check
  if (typeof(data) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || N <= 0 ||
      !(0 < plmax & plmax <= 1) || !(0 < psmax & psmax <= 1) || !(0 < distl & distl < 1) || !(0 < dists & dists < 1) ||
      !(0 < conf.level & conf.level <= 1) || !(typeof(minsplit) %in% c("logical","double")) || 
      !(typeof(minbucket) %in% c("logical", "double")) || typeof(valdat) != "list" || !(1 <= steps & steps <= 3) ) {
    stop("irregular input")
  }
  if (!(all(names(data) %in% names(valdat)))){
    stop("validation data variables unequal input data variables")
  }
  if ((is.na(minsplit) == TRUE) & (is.na(minbucket) == TRUE)){
    minsplit <- 20
    minbucket <- 7
  }
  if (!(is.na(minsplit) == TRUE) & (is.na(minbucket) == TRUE)){
    minbucket <- minsplit / 3
  }
  if ((is.na(minsplit) == TRUE) & !(is.na(minbucket) == TRUE)){
    minsplit <- minbucket * 3
  }
  if (N < 24){
      stop("Number of repetitions too low (must be > 23)")
    }
  names(data)[names(data)==classname] <- "class"
  names(valdat)[names(valdat)==classname] <- "class"
  if (!(identical(levels(data$class),levels(valdat$class)))){
    stop("levels of input class variable unequal levels of validation class variable")
  }
## parameter setting
incl <- c(distl,distl/2,distl/4) ## grid radius in 1st, 2nd, 3rd grid for large class
incs <- c(dists,dists/2,dists/4)    ## grid radius in 1st, 2nd, 3rd grid for small class
n <- c(N/8,N/4,N/2,N)  ## repetitions at the different stages (smaller version)
# 1st grid
N <- n[1]  # no. of repetitions 
psmall <- c(psmax,max(0.01,(psmax-dists)),max(0.01,(psmax-2*dists)))  ## check against 0
plarge <- c(plmax,max(0.01,(plmax-distl)),max(0.01,(plmax-2*distl)))
message("\n","1st grid")
out <- RePrInDT(data,"class",ctestv,N,plarge=plarge,psmall=psmall,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket,valdat=valdat)
#cat("Accuracies of best trees for different combinations of resampling percentages","\n")
#apply(out$acc1st,3,print)
im <- which.max(out$acc1st[,,5])
acc <- out$acc1st
dim(acc) <- c(dim(acc)[1] * dim(acc)[2],dim(acc)[3])
accbest <- acc[im,]
percl <- acc[im,1]  # undersampling percentage of the larger class
percs <- acc[im,2] # undersampling percentage of the smaller class
#message("Accuracies of best tree of resampling percentages","\n")
#message(acc[im,])
# 2nd grid
if (steps > 1){
N <- n[2]
psmall <- c(max(0.01,acc[im,2]-incs[2]),acc[im,2],min(1,acc[im,2]+incs[2])) ## check against 0 and 1
plarge <- c(max(0.01,acc[im,1]-incl[2]),acc[im,1],min(1,acc[im,1]+incl[2]))
message("\n","2nd grid")
out <- RePrInDT(data,"class",ctestv,N,plarge=plarge,psmall=psmall,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket,valdat=valdat)
#cat("Accuracies of best trees for different combinations of resampling percentages","\n")
#apply(out$acc1st,3,print)
im <- which.max(out$acc1st[,,5])
acc <- out$acc1st
dim(acc) <- c(dim(acc)[1] * dim(acc)[2],dim(acc)[3])
accbest <- acc[im,]
percl <- acc[im,1]  # undersampling percentage of the larger class
percs <- acc[im,2] # undersampling percentage of the smaller class
#message("Accuracies of best tree of resampling percentages","\n")
#message(acc[im,])
}
# 3rd grid
if (steps > 2){
N <- n[3]
psmall <- c(max(0.01,acc[im,2]-incs[3]),acc[im,2],min(1,acc[im,2]+incs[3]))  ## check against 0 and 1
plarge <- c(max(0.01,acc[im,1]-incl[3]),acc[im,1],min(1,acc[im,1]+incl[3]))
message("\n","3rd grid")
out <- RePrInDT(data,"class",ctestv,N,plarge=plarge,psmall=psmall,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket,valdat=valdat)
#cat("Accuracies of best trees for different combinations of resampling percentages","\n")
#apply(out$acc1st,3,print)
im <- which.max(out$acc1st[,,5])
acc <- out$acc1st
dim(acc) <- c(dim(acc)[1] * dim(acc)[2],dim(acc)[3])
accbest <- acc[im,]
percl <- acc[im,1]  # undersampling percentage of the larger class
percs <- acc[im,2] # undersampling percentage of the smaller class
#message("Accuracies of best tree of resampling percentages","\n")
#message(acc[im,])
}
##########################################
## 3. Best model
############################
N <- n[4]  # no. of repetitions 
## call of PrInDT
message("\n","Final evaluation")
out <- PrInDT(data,"class",ctestv,N,percl,percs,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket,valdat=valdat)
accbest <- out$ba1st
result <- list(tree1st=out$tree1st,bestba=accbest,percl=percl,percs=percs)
class(result) <- "OptPrInDT"
result
}
#' @export
print.OptPrInDT <- function(x,...){
#  NextMethod()
  print(x$tree1st)
  cat("\n")
  cat("Best balanced accuracy on full sample: ",x$bestba[3],"\n")
  cat("Best undersampling percentages: large class:",x$percl,", small class:",x$percs,"\n\n")
}
#' @export
plot.OptPrInDT <- function(x,...){
#  NextMethod()
  plot(x$tree1st,main="Best interpretable classification tree")
  invisible(x)
}