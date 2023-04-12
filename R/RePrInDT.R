#' Repeated \code{\link{PrInDT}} for specified percentage combinations
#'
#' @description
#' \code{\link{PrInDT}} is called repeatedly according to the percentages specified in the vectors 'plarge' and 
#' 'psmall'.\cr
#' The relationship between the two-class factor variable 'classname' and all other factor and numerical variables
#' in the data frame 'datain' is optimally modeled by means of 'N' repetitions of undersampling.\cr 
#' The trees generated from undersampling can be restricted by rejecting 
#' unacceptable trees which include split results specified in the character strings of the vector 'ctestv'.\cr
#' The probability threshold 'thres' for the prediction of the smaller class may be specified (default = 0.5).\cr
#' Undersampling may be stratified in two ways by the feature 'strat'.
#'
#' \strong{Reference}\cr Weihs, C., Buschfeld, S. 2021c. Repeated undersampling in PrInDT (RePrInDT): Variation in undersampling and prediction, 
#' and ranking of predictors in ensembles. arXiv:2108.05129
#'
#' @usage RePrInDT(datain, classname, ctestv=NA, N, plarge, psmall, conf.level=0.95,
#'        thres=0.5, stratvers=0, strat=NA, seedl=TRUE)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     {see function \code{\link{PrInDT}} for details.}\cr
#'     If no restrictions exist, the default = NA is used.
#' @param N Number of repetitions (integer > 0)
#' @param plarge Vector of undersampling percentages of larger class (numerical, between 0 and 1)
#' @param psmall Vector of undersampling percentages of smaller class (numerical, between 0 and 1)
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical between 0 and 1);\cr
#'     default = 0.95
#' @param thres Probability threshold for prediction of smaller class; default = 0.5
#' @param stratvers Version of stratification;\cr
#'     = 0: none (default),\cr
#'     = 1: stratification according to the percentages of the values of the factor variable 'strat',\cr
#'     > 1: stratification with minimum number 'stratvers' of observations per value of 'strat'
#' @param strat Name of one (!) stratification variable for undersampling (character);\cr
#'     default = NA (no stratification)
#' @param seedl Should the seed for random numbers be set (TRUE / FALSE)?\cr
#'     default = TRUE
#'
#' @return
#' \describe{
#' \item{treesb}{best trees for the different percentage combinations; refer to an individual tree as \code{treesb[[k]]}, k = 1, ..., length(plarge)*length(psmall)}
#' \item{acc1st}{accuracies of best trees on full sample}
#' \item{acc3en}{accuracies of ensemble of 3 best trees on full sample}
#' \item{simp_m}{mean of permutation losses for the predictors}
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.\cr
#' The plot function will produce a series of more than one plot. If you use R, you might want to specify \code{windows(record=TRUE)} before 
#' \code{plot(name)} to save the whole series of plots. In R-Studio this functionality is provided automatically.
#'
#' @exportS3Method print RePrInDT
#' @exportS3Method plot RePrInDT
#' @export RePrInDT
#' @importFrom graphics barplot

#' @examples
#' datastrat <- PrInDT::data_zero
#' data <- na.omit(datastrat) # cleaned full data: no NAs
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}', 'MLU == {1, 3}')
#' N <- 51  # no. of repetitions
#' conf.level <- 0.99 # 1 - significance level (mincriterion) in ctree
#' psmall <- c(0.95,1)     # percentages of the small class
#' plarge <- c(0.09,0.1)  # percentages of the large class
#' outRe <- RePrInDT(data,"real",ctestv,N,plarge,psmall,conf.level) # might take 5 minutes
#' outRe
#' plot(outRe)
#'
RePrInDT <- function(datain,classname,ctestv=NA,N,plarge,psmall,conf.level=0.95,thres=0.5,stratvers=0,strat=NA,seedl=TRUE){
  ## input check
  if (typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || N <= 0 ||
      !(all(0 <= plarge & plarge <= 1)) || !(all(0 <= psmall & psmall <= 1)) ||
      !(0 <= conf.level & conf.level <= 1) | !(0 <= thres & thres <= 1) |
      !(0 <= stratvers) || !(typeof(strat) %in% c("logical", "character")) || typeof(seedl) != "logical"){
    stop("irregular input")
  }
  if (seedl == TRUE){
    set.seed(7654321)  # set seed of random numbers
  }
  ###
  data <- datain
  ## initializations
  norep <- length(psmall) * length(plarge)
  acc1st <- array(0,dim=c(length(plarge),length(psmall),8))
  acc3en <- array(0,dim=c(length(plarge),length(psmall),6)) # !!!
  simp <- array(0,dim=c(length(plarge),length(psmall),(dim(data)[2])) )
  ## permuting columns
  data_per <- data
  for (j in 2:(dim(data)[2])) {
    data_per[,j] <- sample(data[,j])
  }
  k <- 0
  for (i in 1:length(plarge)) {
    for (j in 1:length(psmall)) {
      k <- k + 1
      message("\n")
#      message("\n")
#      message("sampling percentage for larger class: ",plarge[i],"\n")
#      message("sampling percentage for smaller class:",psmall[j],"\n")
      message("sampling percentage for larger class: ",plarge[i])
      message("sampling percentage for smaller class:",psmall[j])
      ## call of PrInDT
      out <- PrInDT(data,classname,ctestv,N,plarge[i],psmall[j],conf.level,thres,stratvers,strat,seedl=TRUE)
      if (any(is.na(out)) == TRUE){ # !!!
        return(NA)
      }
      acc1st[i,j,] <- c(plarge[i],psmall[j],out$ba1st)
      acc3en[i,j,] <- c(plarge[i],psmall[j],out$baen[2,2:5])
      if (k == 1){
        treesb <- out$tree1st
        n_class1 <- table(out$dataout$class)[1] # no. of elements of larger class 1
        n_class2 <- table(out$dataout$class)[2] # no. of elements of smaller class 2
      } else {
        treesb <- c(treesb,out$tree1st)
      }
      
      if (out$ba1st[3] > 0){
        for (l in 2:dim(data)[2]){
          data_imp <- out$dataout
          data_imp[,l] <- data_per[,l]
          ctpreds_imp <- predict(out$tree1st,newdata=data_imp)
          conf_imp <- table(ctpreds_imp, data_imp$class)
          simp[i,j,l] <- simp[i,j,l] + out$ba1st[3] - (conf_imp[1,1] / n_class1 + conf_imp[2,2] / n_class2)/2
          ctpreds_imp <- predict(out$tree2nd,newdata=data_imp)
          conf_imp <- table(ctpreds_imp, data_imp$class)
          simp[i,j,l] <- simp[i,j,l] + out$ba2nd[3] - (conf_imp[1,1] / n_class1 + conf_imp[2,2] / n_class2)/2
          ctpreds_imp <- predict(out$tree3rd,newdata=data_imp)
          conf_imp <- table(ctpreds_imp, data_imp$class)
          simp[i,j,l] <- simp[i,j,l] + out$ba3rd[3] - (conf_imp[1,1] / n_class1 + conf_imp[2,2] / n_class2)/2
        }
      }
    }
  }
## preparation of print
  dimnames(acc1st) <- list(1:length(plarge),1:length(psmall),c("plarge","psmall",paste0("full ",levels(out$dataout$class)[1]),
                                                               paste0("full ",levels(out$dataout$class)[2]),"full balanced",
                                                               paste0("test ",levels(out$dataout$class)[1]),paste0("test ",levels(out$dataout$class)[2]),"test balanced"))
  dimnames(acc3en) <- list(1:length(plarge),1:length(psmall),c("plarge","psmall",paste0("full ",levels(out$dataout$class)[1]),
                                                               paste0("full ",levels(out$dataout$class)[2]),"full balanced","mean test balanced")) # !!!
  simp_m <- sapply(2:7,function(x) mean(simp[,,x]))
  names(simp_m) <- colnames(out$dataout)[-(colnames(out$dataout) == "class")]
###
  result <- list(treesb=treesb,acc1st=acc1st,acc3en=acc3en,simp_m=simp_m)
  class(result) <- "RePrInDT"
  result
}
###
print.RePrInDT <- function(x, ...){
  ####
  ## print accuracies of best trees
  ####
  cat("\n","Accuracies of best trees for different combinations of resampling percentages: on full and test samples","\n")
  apply(x$acc1st,1,print)
  ## print overall best trees
  cat("\n\n","Best trees of percentage combinations","\n")
  maxba <- as.vector(t(x$acc1st[,,5] == max(x$acc1st[,,5])))
  mm <- which(maxba)
  k <- 0 
  for (i in 1:dim(x$acc1st)[1]) {
    for (j in 1:dim(x$acc1st)[2]) {
      k <- k + 1
      if (k %in% mm) {
        cat("\n","tree for  ( plarge  psmall  bal.acc. ) = (",x$acc1st[i,j,c(1,2,5)],")")
        print(x$treesb[[k]])
      }
    }
  }
  ## print accuracies of ensembles of 3 best trees
  cat("\n","Accuracies of ensembles of 3 best trees for different combinations of resampling percentages","\n")
  apply(x$acc3en,1,print)
  ####
  ## ranking of predictors
  ####
  cat("\n\n","Mean of permutation losses for the predictors","\n")
  print(x$simp_m)
}
###
plot.RePrInDT <- function(x, ...){
  ####
  ## plot best trees
  ####
  maxba <- as.vector(t(x$acc1st[,,5] == max(x$acc1st[,,5])))
  mm <- which(maxba)
  lmm <- length(mm)
  k <- 0 
  for (i in 1:dim(x$acc1st)[1]) {
    for (j in 1:dim(x$acc1st)[2]) {
      k <- k + 1
      if (k %in% mm) {
        v <- toString(x$acc1st[i,j,c(1,2,5)])
        plot(x$treesb[[k]],main=paste0("One of the ",lmm," best trees: ( plarge,  psmall,  bal.acc. ) = (",v,")"))
      }
    }
  }
  ####
  ## ranking of predictors
  ####
  cat("\n")
  barplot(sort(x$simp_m / max(x$simp_m)),horiz=TRUE,main="Normed means of permutation losses")
}