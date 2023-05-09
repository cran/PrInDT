#' Multiple label classification based on all observations
#'
#' @description Multiple label classification based on all observations. We consider two ways of modeling (Binary relevance modeling, 
#' dependent binary modeling) and three ways of model evaluation: single 
#' assessment, joint assessment, and true prediction (see the Value section for more information).\cr
#' Interpretability is checked (see ctestv).\cr
#' Variables should be arranged in 'datain' according to indices specified in 'indind', 'indaddind', and 'inddep'.\cr
#'
#' \strong{Reference}\cr Probst, P., Au, Q., Casalicchio, G., Stachl, C., and Bischl, B. 2017. Multilabel Classification with 
#' R Package mlr. arXiv:1703.08991v2
#'
#' @usage PrInDTMulabAll(datain, classnames, ctestv=NA, conf.level=0.95, indind, indaddind,
#'        inddep)
#' 
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classnames names of class variables (character vector)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     {see function \code{\link{PrInDT}} for details.}\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'   default = 0.95
#' @param indind indices of independent variables
#' @param indaddind indices of additional predictors used in the case of dependent binary relevance modeling
#' @param inddep indices of dependent variables
#'
#' @return
#' \describe{ 
#' \item{accabr}{model errors for Binary Relevance (single assessment) - only independent predictors are used for modeling one label at a time,  
#' the other labels are not used as predictors. The classification rules are trained on all observations. 
#' As the performance measure for the resulting classification rules, the balanced accuracy of the models for each individual label is employed.}
#' \item{errabin}{combined error for Binary Relevance (joint assessment) - the best prediction models for the different labels are combined to assess the 
#' combined prediction. The 01-accuracy counts a label combination as correct only if all labels are correctly predicted. 
#' The hamming accuracy corresponds to the proportion of labels whose value is correctly predicted.}
#' \item{accadbr}{model errors in Dependent Binary Relevance (Extended Model) (single assessment) - each label is trained by means of an extended model which 
#' not only includes the independent predictors but also the other labels. For these labels the truly observed values are used for estimation 
#' and prediction. In the extended model, further labels, which are not treated as dependent variables, can be used as additional predictors.}
#' \item{erraext}{combined errors for Dependent Binary Relevance (Extended Model) (joint assessment) }
#' \item{erratrue}{combined errors for Dependent Binary Relevance (True Prediction) - in the prediction phase, the values 
#' of all modeled labels are first predicted by the independent predictors only (see Binary Relevance) and then the predicted labels are used in  
#' the estimated extended model in a 2nd step to ultimately predict the labels.}
#' \item{coldata}{column names of input data}
#' \item{inddep}{indices of dependent variables (labels to be modeled)}
#' \item{treeabr}{list of trees from Binary Relevance modeling, one tree for each label; refer to an individual tree as \code{treeabr[[i]]}, 
#'                 i = 1, ..., no. of labels}
#' \item{treeadbr}{list of trees from Dependent Binary Relevance modeling, one for each label; refer to an individual tree as \code{treeadbr[[i]]}, 
#'                 i = 1, ..., no. of labels}
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)}  where 'name' is the output data 
#' frame of the function.\cr
#' The plot function will produce a series of more than one plot. If you use R, you might want to specify \code{windows(record=TRUE)} before 
#' \code{plot(name)} to save the whole series of plots. In R-Studio this functionality is provided automatically.
#'
#' @export PrInDTMulabAll
#' @exportS3Method print PrInDTMulabAll
#' @exportS3Method plot PrInDTMulabAll
#'
#' @examples
#' data <- PrInDT::data_land # load data
#' dataclean <- data[,c(1:7,23:24,11:13,22,8:10)]  # only relevant features
#' indind <- c(1:9) # original predictors
#' indaddind <- c(10:13) # additional predictors
#' inddep <- c(14:16) # dependent variables
#' dataclean <- na.omit(dataclean)
#' ctestv <- NA
#' ##
#' # Call PrInDTAll: language by language
#' ##
#' outmultAll <- PrInDTMulabAll(dataclean,colnames(dataclean)[inddep],ctestv,conf.level=0.95,
#'                      indind,indaddind,inddep)
#' outmultAll
#' plot(outmultAll)
#'
#' @import party
#' @import stats
#'
PrInDTMulabAll <- function(datain,classnames,ctestv=NA,conf.level=0.95,indind,indaddind,inddep){
  ## input check
  if ( typeof(datain) != "list" || typeof(classnames) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || 
      !(0 < conf.level & conf.level <= 1) || !(typeof(indind) %in% c("integer", "double")) || 
      !(typeof(indaddind) %in% c("integer", "double")) || !(typeof(inddep) %in% c("integer", "double")) ){
    stop("irregular input")
  }
####
####
## Modeling with no additional independent and no dependent variable as predictors
####
  data <- datain
#  message("\n")
  message("*** Binary relevance ***")
#  message("\n")
  message("Modeling the labels:")
  accabr <- matrix(0,nrow=length(inddep),ncol=1)
  rownames(accabr) <- vector(mode="character",length=length(inddep))
  colnames(accabr) <- "balanced"
  ctpreds <- matrix(0,ncol=length(inddep),nrow=dim(datain)[1])
  unequal <- matrix(TRUE,ncol=length(inddep),nrow=dim(datain)[1])
  hsum <- 0
  for (i in 1:length(inddep)){
    x <- colnames(datain)[inddep[i]]
    levels(data[,x]) <- c(levels(data[,x]),x)
    data[datain[,x] == 1,x] <- x
    data[,x] <- droplevels(data[,x])
#    message("\n")
    message(paste0("  ",colnames(datain)[inddep[i]]))
    out <- PrInDTAll(data[,c(indind,inddep[i])],x,ctestv,conf.level) # minimum call
    rownames(accabr)[i] <- x
    if (i == 1) {
      treeabr <- list(out$treeAll)
    } else {
      treeabr <- c(treeabr,out$treeAll)
    }
    accabr[i] <- out$baAll
    ctpreds[,i] <- stats::predict( treeabr[[i]],newdata=data[,c(indind,inddep[i])] )
    unequal[,i] <- as.integer(data[,inddep[i]]) != ctpreds[,i]
    if (table(data[,inddep[i]])[1] < table(data[,inddep[i]])[2] ){
      unequal[,i] <- as.integer(data[,inddep[i]]) == ctpreds[,i]
    }
    hsum <- hsum + sum(unequal[,i])
  }
# Joint assessment: all
  errabin <- c(0,0)
  # 01 loss
  err01 <- rowSums(unequal)
  err01[err01 > 1] <- 1
  errabin[1] <- sum(err01)
  errabin[1] <- errabin[1] / dim(data)[1]
  # hamming loss
  errabin[2] <- hsum / (3*dim(data)[1])
  names(errabin) <- c("01-accuracy", "hamming-accuracy")
####
## Modeling with additional independent and dependent variable as predictors
####
#  message("\n")
  message("\n","*** Dependent binary relevance: Extended model ***")
#  message("\n")
  message("Modeling the labels:")
  accadbr <- matrix(0,nrow=length(inddep),ncol=1)
  rownames(accadbr) <- vector(mode="character",length=length(inddep))
  colnames(accadbr) <- "balanced"
  hsum <- 0
  for (i in 1:length(inddep)) {
    x <- colnames(datain)[inddep[i]]
#   message("\n")
    message(paste0("  ",x))
    out <- PrInDTAll(data,x,ctestv,conf.level) # minimum call
    rownames(accadbr)[i] <- x
    if (i == 1) {
      treeadbr <- list(out$treeAll)
    } else {
      treeadbr <- c(treeadbr,out$treeAll)
    }
    accadbr[i] <- out$baAll
    ctpreds[,i] <- stats::predict( treeadbr[[i]],newdata=data )
    unequal[,i] <- as.integer(data[,inddep[i]]) != ctpreds[,i]
    if (table(data[,inddep[i]])[1] < table(data[,inddep[i]])[2] ){
      unequal[,i] <- as.integer(data[,inddep[i]]) == ctpreds[,i]
    }
    hsum <- hsum + sum(unequal[,i])
  }
#  message("\n")
# Joint assessment: all
  erraext <- c(0,0)
  # 01 loss
  err01 <- rowSums(unequal)
  err01[err01 > 1] <- 1
  erraext[1] <- sum(err01)
  erraext[1] <- erraext[1] / dim(data)[1]
  # hamming loss
  erraext[2] <- hsum / (3*dim(data)[1])
  names(erraext) <- c("01-accuracy", "hamming-accuracy")
  ##
  ## nested
  ##
# Dependent binary relevance: True prediction
  datain2 <- data
  ctpreds2 <- matrix(0,ncol=length(inddep),nrow=dim(data)[1])
  unequal2 <- matrix(TRUE,ncol=length(inddep),nrow=dim(data)[1])
  hsum2 <- 0
  for (i in 1:length(inddep)){
    #  p <- ctpreds[,i] - 1
    datain2[,inddep[i]] <- stats::predict(treeabr[[i]],newdata=data) # as.factor(p)
    if (table(data[,inddep[i]])[1] < table(data[,inddep[i]])[2] ){
      datain2[,inddep[i]] <- relevel(datain2[,inddep[i]],"0")
    }
  }
  for (i in 1:length(inddep)){
    ctpreds2[,i] <- stats::predict(treeadbr[[i]],newdata=datain2)
    if (table(data[,inddep[i]])[1] >= table(data[,inddep[i]])[2] ){
      unequal2[,i] <- as.integer(data[,inddep[i]]) != ctpreds2[,i]
    } else {
      unequal2[,i] <- as.integer(data[,inddep[i]]) == ctpreds2[,i]
    }
    hsum2 <- hsum2 + sum(unequal2[,i])
  }
  erratrue <- c(0,0)
  # 01 loss
  err01 <- rowSums(unequal2)
  err01[err01 > 1] <- 1
  erratrue[1] <- sum(err01)
  erratrue[1] <- erratrue[1] / dim(data)[1]
  # hamming loss
  erratrue[2] <- hsum2 / (3*dim(data)[1])
  names(erratrue) <- c("01-accuracy", "hamming-accuracy")
  result <- list(accabr = accabr, errabin = errabin, accadbr = accadbr, erraext = erraext, erratrue = erratrue, 
     coldata = colnames(datain), inddep = inddep, treeabr = treeabr, treeadbr = treeadbr)
  class(result) <- "PrInDTMulabAll"
  result
}

## print
print.PrInDTMulabAll <- function(x,...){
  cat("\n","Multi-label classification on full sample","\n")
  cat("\n")
  cat("*** Binary relevance ***")
  cat("\n")
  L <- length(x$inddep)
  for (i in 1:L){
    cat("\n")
    cat( "Binary Relevance: Best tree on subsamples for: ",x$coldata[x$inddep[i]] )
    print( x$treeabr[[i]] )
  }
  ### single assessment: all
  cat("\n\n")
  cat("Single assessment: all","\n")
  print(x$accabr)
  ### joint assessment all
  cat("\n")
  cat("Joint assessment: all")
  cat("\n","Assessment of models on full sample (labels are NOT predictors)","\n")
  print(1-x$errabin)
  cat("\n\n")
  cat("*** Dependent binary relevance: Extended Model ***")
  cat("\n")
  for (i in 1:L){
    cat("\n")
    cat( "Dependent Binary Relevance: Best tree on full sample for: ",x$coldata[x$inddep[i]] ) 
    print( x$treeadbr[[i]] )
  }
  ### single assessment: all
  cat("\n\n")
  cat("Single assessment: all","\n")
  print(x$accadbr)
  ### joint assessment all
  cat("\n")
  cat("Joint assessment: all")
  cat("\n","Assessment of models on full sample (labels are predictors)","\n")
  print(1-x$erraext)
  ##
  ## nested
  ##
  cat("\n\n")
  cat("*** Dependent binary relevance: True prediction ***")
  cat("\n")
  cat("\n","Assessment of PrInDT models (true prediction)","\n")
  print(1-x$erratrue)
}


## plot
plot.PrInDTMulabAll <- function(x,...){
L <- length(x$inddep)
  for (i in 1:L){
    plot( x$treeabr[[i]],main=paste0("Binary Relevance on full sample: Tree on all data for: ",x$coldata[x$inddep[i]]) )
  }
  for (i in 1:L){
    plot( x$treeadbr[[i]],main=paste0("Dependent Binary Relevance on full sample: Tree on all data for: ",x$coldata[x$inddep[i]]) )
  }
}
