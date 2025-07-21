#' Multiple label classification based on resampling by \code{\link{PrInDT}}
#'
#' @description Multiple label classification based on resampling by \code{\link{PrInDT}}. We consider two ways of modeling (Binary relevance modeling, 
#' dependent binary modeling) and three ways of model evaluation: single assessment, joint assessment, and true prediction  
#' (see the Value section for more information).\cr
#' Variables should be arranged in 'datain' according to indices specified in 'indind', 'indaddind', and 'inddep'.\cr
#' Please note that the dependent variables have to be specified as dummies, i.e. as 'absent' (value 0) or 'present' (value 1).\cr
#' Undersampling is repeated 'N' times.\cr
#' Undersampling percentages 'percl' for the larger class and 'percs' for the smaller class can be 
#' specified, one each per dependent class variable.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' \strong{References}\cr 
#' - Buschfeld, S., Weihs, C. and Ronan, P.. 2024. Modeling linguistic landscapes: A comparison of St Martin’s two capitals Philipsburg and Marigot Linguistic Landscape 10(3): 302–334. https://doi.org/10.1075/ll.23070.bus\cr
#' - Probst, P., Au, Q., Casalicchio, G., Stachl, C., and Bischl, B. 2017. Multilabel Classification with R Package mlr. arXiv:1703.08991v2
#'
#' @usage PrInDTMulab(datain,classnames=NA,ctestv=NA,conf.level=0.95,percl,percs=1,N,
#'                indind=NA,indaddind=NA,inddep,minsplit=NA,minbucket=NA)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classnames names of class variables (character vector); default = NA: from old version: superfluous now
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     (see function \code{\link{PrInDT}} for details.)\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#' default = 0.95
#' @param percl list of undersampling percentages of larger class (numerical, > 0 and <= 1): one per dependent class variable
#' @param percs list of undersampling percentage of smaller class (numerical, > 0 and <= 1); one per dependent class variable
#' @param N no. of repetitions (integer > 0)
#' @param indind indices of independent variables
#' @param indaddind indices of additional independent variables used in the case of dependent binary relevance modeling
#' @param inddep indices of dependent variables
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#'
#' @return
#' \describe{
#' \item{accbr}{model errors for Binary Relevance (single assessment) - only independent predictors are used for modeling one label at a time,
#' the other labels are not used as predictors. As the performance measure for the resulting classification rules,
#' the balanced accuracy of the best model from PrInDT is employed for each individual label.}
#' \item{errbin}{combined error for Binary Relevance (joint assessment) - the best prediction models for the different labels are combined  
#' to assess the combined prediction. The 01-accuracy counts a label combination as correct only if all labels are correctly predicted.
#' The hamming accuracy corresponds to the proportion of labels whose value is correctly predicted.}
#' \item{accdbr}{model errors for Dependent Binary Relevance (Extended Model) (single assessment) - each label is trained by means of an extended model which 
#' not only includes the independent predictors but also the other labels. For these labels, the truly observed values are used for 
#' estimation and prediction. In the extended model, other labels, which are not treated as dependent variables, can also be used as additional predictors.}
#' \item{errext}{combined errors for Dependent Binary Relevance (Extended Model) (joint assessment)}
#' \item{errtrue}{combined errors for Dependent Binary Relevance (True Prediction) - in the prediction phase, the values 
#' of all modeled labels are first predicted by the independent predictors only and then the predicted labels are used in the estimated  
#' extended model in a 2nd step to 
#' ultimately predict the labels.}
#' \item{coldata}{column names of input data}
#' \item{inddep}{indices of dependent variables (labels to be modeled)}
#' \item{treebr}{list of trees from Binary Relevance modeling, one tree for each label; refer to an individual tree as \code{treebr[[i]]}, 
#'                 i = 1, ..., no. of labels}
#' \item{treedbr}{list of trees from Dependent Binary Relevance modeling, one for each label; refer to an individual tree as \code{treedbr[[i]]}, 
#'                 i = 1, ..., no. of labels}
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)}  where 'name' is the output data 
#' frame of the function.\cr
#' The plot function will produce a series of more than one plot. If you use R, you might want to specify \code{windows(record=TRUE)} before 
#' \code{plot(name)} to save the whole series of plots. In R-Studio this functionality is provided automatically.
#'
#' @export PrInDTMulab
#'
#' @examples
#' data <- PrInDT::data_land # load data
#' dataclean <- data[,c(1:7,23:24,11:13,22,8:10)]  # only relevant features
#' indind <- c(1:9) # original predictors
#' indaddind <- c(10:13) # additional predictors
#' inddep <- c(14:16) # dependent variables
#' dataclean <- na.omit(dataclean)
#' ctestv <- NA
#' N <- 21  # no. of repetitions
#' perc <- c(0.45,0.05,0.25)   # percentages of observations of larger class, 
#' #                             1 per dependent class variable
#' perc2 <- c(0.75,0.95,0.75)  # percentages of observations of smaller class, 
#' #                             1 per dependent class variable
#' ##
#' # Call PrInDT: language by language
#' ##
#' outmult <- PrInDTMulab(dataclean,ctestv=NA,conf.level=0.95,percl=perc,percs=perc2,N=N,
#'                   indind=indind,indaddind=indaddind,inddep=inddep)
#' print(outmult)
#' plot(outmult)
#'
#' @import party
#' @import stats
#'
PrInDTMulab <- function(datain,classnames=NA,ctestv=NA,conf.level=0.95,percl,percs=1,N,indind=NA,indaddind=NA,inddep,minsplit=NA,minbucket=NA){
  ## input check
  if (typeof(datain) != "list" || !(typeof(classnames) %in% c("logical", "character")) || !(typeof(ctestv) %in% c("logical", "character")) || N <= 0 ||
      !all(0 < percl & percl <= 1) || !all(0 < percs & percs <= 1) || !(0 < conf.level & conf.level <= 1) || 
      !(typeof(indind) %in% c("integer", "double","logical")) || 
      !(typeof(indaddind) %in% c("integer", "double","logical")) || !(typeof(inddep) %in% c("integer", "double")) || !(typeof(minsplit) %in% c("logical","double")) || 
      !(typeof(minbucket) %in% c("logical", "double"))  ) {
    stop("irregular input")
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
####
## Modeling with no additional independent and no dependent variable as predictors
####
  data <- datain
  classnames <- colnames(data)[inddep]
  if (all(is.na(indind)) == TRUE & all(is.na(indaddind)) == TRUE){
    indind <- c(1:dim(data)[2])[-inddep]
  }  
  if (all(is.na(indind)) == TRUE & all(is.na(indaddind)) != TRUE){
    indind <- c(1:dim(data)[2])[-c(inddep,indaddind)]
  }  
#  print(indind)
#  message("\n\n")
  message("*** Binary relevance ***")
#  message("\n")
  message("Modeling the labels:")
  accbr <- matrix(0,nrow=length(inddep),ncol=2)
  rownames(accbr) <- vector(mode="character",length=length(inddep))
  colnames(accbr) <- c("full balanced","test balanced")
  ctpreds <- matrix(0,ncol=length(inddep),nrow=dim(datain)[1])
  unequal2 <- matrix(TRUE,ncol=length(inddep),nrow=dim(datain)[1])
  hsum2 <- 0
  for (i in 1:length(inddep)){
    x <- colnames(datain)[inddep[i]]
    levels(data[,x]) <- c(levels(data[,x]),x)
    data[datain[,x] == 1,x] <- x
    data[,x] <- droplevels(data[,x])
#    message("\n")
    message(paste0("  ",x))
    perc <- percl[i]  # undersampling percentage of the larger class
    perc2 <- percs[i] # undersampling percentage of the smaller class
    ## call of PrInDT
    out <- PrInDT(data[,c(indind,inddep[i])],x,ctestv,N,perc,perc2,conf.level,minsplit=minsplit,minbucket=minbucket)
    rownames(accbr)[i] <- colnames(datain)[inddep[i]]
    if (i == 1) {
      treebr <- list(out$tree1st)
      outtreet <- list(out$treet1st)
    } else {
      treebr <- c(treebr,out$tree1st)
      outtreet <- c(outtreet,out$treet1st)
    }
    accbr[i,] <- out$ba1st[c(3,6)]
    #  dev.off()
    ctpreds[,i] <- stats::predict( treebr[[i]],newdata=datain[,c(indind,inddep[i])] )
    unequal2[,i] <- as.integer(datain[,inddep[i]]) != ctpreds[,i]
    if (table(datain[,inddep[i]])[1] < table(datain[,inddep[i]])[2] ){
      unequal2[,i] <- as.integer(datain[,inddep[i]]) == ctpreds[,i]
    }
    hsum2 <- hsum2 + sum(unequal2[,i])
  }
  ### joint assessment: PrInDT
  errbin <- c(0,0)
  # 01 loss
  err01 <- rowSums(unequal2)
  err01[err01 > 1] <- 1
  errbin[1] <- sum(err01)
  errbin[1] <- errbin[1] / dim(datain)[1]
  # hamming loss
  errbin[2] <- hsum2 / (3*dim(datain)[1])
  names(errbin) <- c("01-accuracy", "hamming-accuracy")
####
## Modeling with additional independent and dependent variable as predictors
####
#  message("\n")
  message("\n","*** Dependent binary relevance: Extended model ***")
#  message("\n")
  message("Modeling the labels:")
  accdbr <- matrix(0,nrow=length(inddep),ncol=2)
  rownames(accdbr) <- vector(mode="character",length=length(inddep))
  colnames(accdbr) <- c("full balanced","test balanced")
  ctpreds2 <- matrix(0,ncol=length(inddep),nrow=dim(datain)[1])
  unequal2 <- matrix(TRUE,ncol=length(inddep),nrow=dim(datain)[1])
  hsum2 <- 0
  for (i in 1:length(inddep)) {
    x <- colnames(data)[inddep[i]]
#    message("\n")
    message(paste0("  ",x))
    perc <- percl[i]  # undersampling percentage of the larger class
    perc2 <- percs[i] # undersampling percentage of the smaller class
    ## call of PrInDT
    out <- PrInDT(data,x,ctestv,N,perc,perc2,conf.level,minsplit=minsplit,minbucket=minbucket)
    rownames(accdbr)[i] <- x
    if (i == 1) {
      treedbr <- list(out$tree1st)
      outtreetl <- list(out$treet1st)
    } else {
      treedbr <- c(treedbr,out$tree1st)
      outtreetl <- c(outtreetl,out$treet1st)
    }
    accdbr[i,] <- out$ba1st[c(3,6)]
    # dev.off()
    ctpreds2[,i] <- stats::predict( treedbr[[i]],newdata=data )
    unequal2[,i] <- as.integer(data[,inddep[i]]) != ctpreds2[,i]
    if (table(data[,inddep[i]])[1] < table(datain[,inddep[i]])[2] ){
      unequal2[,i] <- as.integer(data[,inddep[i]]) == ctpreds2[,i]
    }
    hsum2 <- hsum2 + sum(unequal2[,i])
  }
# Joint assessment: PrInDT
  errext <- c(0,0)
  # 01 loss
  err01 <- rowSums(unequal2)
  err01[err01 > 1] <- 1
  errext[1] <- sum(err01)
  errext[1] <- errext[1] / dim(data)[1]
  # hamming loss
  errext[2] <- hsum2 / (3*dim(data)[1])
  names(errext) <- c("01-accuracy", "hamming-accuracy")
  ##
  ## nested: True prediction
  ##
  datain2 <- data
  ctpreds2 <- matrix(0,ncol=length(inddep),nrow=dim(data)[1])
  unequal2 <- matrix(TRUE,ncol=length(inddep),nrow=dim(data)[1])
  accdbrt <- matrix(0,nrow=length(inddep),ncol=1)
  rownames(accdbrt) <- classnames
  hsum2 <- 0
  for (i in 1:length(inddep)){
    #  p <- ctpreds[,i] - 1
    datain2[,inddep[i]] <- stats::predict(treebr[[i]],newdata=data)
    if (table(data[,inddep[i]])[1] < table(data[,inddep[i]])[2] ){
      datain2[,inddep[i]] <- relevel(datain2[,inddep[i]],"0")
    }
  }
  for (i in 1:length(inddep)){
#    rownames(accdbrt)[i] <- colnames(data)[inddep[i]]   ### ?????
    ctpreds2[,i] <- stats::predict(treedbr[[i]],newdata=datain2)
    tabdbrt <- table(data[,inddep[i]],ctpreds2[,i])
#print(table(data[,inddep[i]]))
#print(tabdbrt)
    if (length(colnames(tabdbrt)) == 2){
      if (table(data[,inddep[i]])[1] >= table(data[,inddep[i]])[2] ){
        unequal2[,i] <- as.integer(data[,inddep[i]]) != ctpreds2[,i]
        accdbrt[i] <- ( tabdbrt[1,1] / sum(tabdbrt[1,]) + tabdbrt[2,2] / sum(tabdbrt[2,]) ) / 2
      } else {
        unequal2[,i] <- as.integer(data[,inddep[i]]) == ctpreds2[,i]
        accdbrt[i] <- (tabdbrt[1,2] / sum(tabdbrt[1,]) + tabdbrt[2,1] / sum(tabdbrt[2,]) ) / 2
      }
    } else {
    if (table(data[,inddep[i]])[1] >= table(data[,inddep[i]])[2] ){
        unequal2[,i] <- as.integer(data[,inddep[i]]) != ctpreds2[,i]
        accdbrt[i] <- 0.5  ## tabdbrt[2,1] / sum(tabdbrt[,1])    ## ???
      } else {
        unequal2[,i] <- as.integer(data[,inddep[i]]) == ctpreds2[,i]
        accdbrt[i] <- 0.5  ## tabdbrt[1,1] / sum(tabdbrt[,1])    ## ???
      }
    } 
#print(accdbrt[i]) 
    hsum2 <- hsum2 + sum(unequal2[,i])
  }
  errtrue <- c(0,0)
  # 01 loss
  err01 <- rowSums(unequal2)
  err01[err01 > 1] <- 1
  errtrue[1] <- sum(err01)
  errtrue[1] <- errtrue[1] / dim(data)[1]
  # hamming loss
  errtrue[2] <- hsum2 / (3*dim(data)[1])
  names(errtrue) <- c("01-accuracy", "hamming-accuracy")
#  message("\n")
  result <- list(accbr = accbr, errbin = errbin, accdbr = accdbr, errext = errext, errtrue = errtrue, 
    coldata = colnames(datain), inddep = inddep, treebr = treebr, treedbr = treedbr,accdbrt=accdbrt )
  class(result) <- "PrInDTMulab"
  result
}
#' @export
print.PrInDTMulab <- function(x,...){
  cat("\n","Optimal multi-label classification on subsamples","\n")
  cat("\n")
  cat("*** Binary relevance ***")
  cat("\n")
  L <- length(x$inddep)
  for (i in 1:L){
    cat("\n")
    cat( "Binary Relevance: Best tree on subsamples for label: ",x$coldata[x$inddep[i]] )
    print( x$treebr[[i]] )
  }
  ### single assessment: all
  cat("\n")
  cat("Single assessment: all","\n")
  print(x$accbr)
  ### joint assessment all
  cat("\n")
  cat("Joint assessment: all")
  cat("\n","Assessment of models on full sample (lables are NOT predictors)","\n")
  print(1-x$errbin)
  cat("\n\n")
  cat("*** Dependent binary relevance: Extended Model ***","\n")
  for (i in 1:L){
    cat("\n")
    cat( "Dependent Binary Relevance: Best tree on subsamples for label: ",x$coldata[x$inddep[i]] ) 
    print( x$treedbr[[i]] )
  }
  ### single assessment: all
  cat("\n\n")
  cat("Single assessment: all","\n")
  print(x$accdbr)
  ### joint assessment all
  cat("\n")
  cat("Joint assessment: all")
  cat("\n","Assessment of models on full sample (labels are predictors)","\n")
  print(1-x$errext)
  ##
  ## nested
  ##
  cat("\n\n")
  cat("*** Dependent binary relevance: True prediction ***")
  colnames(x$accdbrt) <- "balanced"
  cat("\n\n")
  cat("Single assessment: all","\n")
  print(x$accdbrt)
  cat("\n")
  cat("\n","Joint assessment all","\n")
  print(1-x$errtrue)
}
#' @export
plot.PrInDTMulab <- function(x,...){
L <- length(x$inddep)
  for (i in 1:L){
    plot( x$treebr[[i]],main=paste0("Binary Relevance: Best tree on subsamples for label: ",x$coldata[x$inddep[i]]) )
  }
  for (i in 1:L){
    plot( x$treedbr[[i]],main=paste0("Dependent Binary Relevance: Best tree on subsamples for label: ",x$coldata[x$inddep[i]]) )
  }
}
