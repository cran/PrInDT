#' PrInDT analysis for a classification problem with multiple classes.

#' @description PrInDT analysis for a classification problem with more than 2 classes. For each combination of one class vs. 
#' the other classes a 2-class \code{\link{PrInDT}} analysis is carried out.\cr
#' The percentages for undersampling of the larger class ('percl' in \code{\link{PrInDT}}) are chosen so that the resulting sizes
#' are comparable with the size of the smaller classes for which all their observations are used in undersampling ('percs' = 1 in \code{\link{PrInDT}}).\cr
#' The class with the highest probability in the K (= number of classes) analyses is chosen for prediction.\cr
#' Interpretability is checked (see 'ctestv').
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage PrInDTMulev(datain,classname,ctestv=NA,N,conf.level=0.95,seedl=FALSE,
#'                                     minsplit=NA,minbucket=NA)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     (see function \code{\link{PrInDT}} for details.)\cr
#'     If no restrictions exist, the default = NA is used.
#' @param N Number of repetitions (integer > 0)
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1)\cr
#'     (default = 0.95)
#' @param seedl Should the seed for random numbers be set (TRUE / FALSE)?\cr
#'     default = FALSE
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#'
#' @return 
#' \describe{
#' \item{class}{levels of class variable} 
#' \item{trees}{trees for the levels of the class variable; refer to an individual tree as \code{trees[[k]]}, k = 1, ..., no. of levels}
#' \item{ba}{balanced accuracy of combined predictions} 
#' \item{conf}{confusion matrix of combined predictions} 
#' \item{ninterp}{no. of non-interpretable trees} 
#' \item{acc}{balanced accuracies of best models for individual classes}
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)}  where 'name' is the output data 
#' frame of the function.\cr
#' The plot function will produce a series of more than one plot. If you use R, you might want to specify \code{windows(record=TRUE)} before 
#' \code{plot(name)} to save the whole series of plots. In R-Studio this functionality is provided automatically.
#'
#' @export PrInDTMulev
#'
#' @examples
#' datastrat <- PrInDT::data_zero
#' data <- na.omit(datastrat)
#' ctestv <- NA
#' data$rel[data$ETH %in% c("C1a","C1b","C1c") & data$real == "zero"] <- "zero1"
#' data$rel[data$ETH %in% c("C2a","C2b","C2c") & data$real == "zero"] <- "zero2"
#' data$rel[data$real == "realized"] <- "real"
#' data$rel <- as.factor(data$rel) # rel is new class variable
#' data$real <- NULL # remove old class variable
#' N <- 51
#' conf.level <- 0.99 # 1 - significance level (mincriterion) in ctree
#' out <- PrInDTMulev(data,"rel",ctestv,N,conf.level) 
#' out # print best models based on subsamples
#' plot(out) # corresponding plots
#'
#' @importFrom stats relevel predict
#' @importFrom party ctree ctree_control
#'
PrInDTMulev <- function(datain,classname,ctestv=NA,N,conf.level=0.95,seedl=FALSE,minsplit=NA,minbucket=NA){
  ## input check
  if (typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || N <= 0 ||
      !(0 < conf.level & conf.level <= 1) || typeof(seedl) != "logical" || !(typeof(minsplit) %in% c("logical","double")) || 
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
  if (seedl == TRUE){
    set.seed(7654321)  # set seed of random numbers
  }
  data <- datain
  names(data)[names(data)==classname] <- "class"
  K <- nlevels(data$class)
  n <- dim(data)[1]
  datastud <- data
  datastud$class <- NULL
  ctpredsk <- rep(0,n)
  ctpreds <- matrix(0,nrow=n,ncol=K)
  crit1 <- rep(FALSE,K)
  acci <- rep(0,K)
##
## loop over the classes (one vs. rest)
##
message("loop over all individual classes vs. rest")
 for (k in 1:K){  # loop over the classes
   message("individual class = ",levels(data$class)[k])
   n_class1 <- table(data$class)[k] # no. of elements of studied class
   n_class2 <- n - n_class1 # no. of elements in all other classes
   if (min(n_class1,n_class2) == 0){
     stop("irregular situation: only 1 class left")
   }
   class1 <- as.character(levels(data$class)[k])
   datastud$classstud[data$class == levels(data$class)[k]] <- class1
   datastud$classstud[!(data$class == levels(data$class)[k])] <- "rest"
   datastud$classstud <- as.factor(datastud$classstud)
   names(datastud$classstud) <- "classstud"
  ## determining the percentage for larger class
   if (n_class1 > n_class2){
     perclin <- n_class2 / (n - n_class2) # * 0.8
   } else {
     perclin <- n_class1 / (n - n_class1) # * 0.8
   }
  ######
  ## call of PrInDT model for one class vs. rest 
  ######
   outin <- PrInDT(datastud,"classstud",ctestv,N,percl=perclin,percs=1,conf.level,minsplit=minsplit,minbucket=minbucket)
#   outin <- PrInDT(datastud,"classstud",ctestv,N,percl=perclin,percs=0.8,conf.level)
   ct <- outin$tree1st
   if (k == 1){
     trees <- ct
   }
   else {
     trees <- c(trees,ct)
   }
   acci[k] <- outin$ba1st[3]
   if (is.na(ctestv[1]) == FALSE) {
     crit1[k] <- FindSubstr(ct,ctestv) # call of the above function for overall tree
   }
   x <- t(array(unlist(stats::predict(ct,type="prob",newdata=datastud)),dim=c(2,n)))
   if (n_class1 < n_class2){
     ctpreds[,k] <- x[,2]
   }
   else {
     ctpreds[,k] <- x[,1]
   }
    datastud$classstud <- NULL
 }
## deployment 
## comparison of probabilities of predicted classes: identify maximum
##
  predclassno <- rep(0,n)
  for (i in 1:n){
    predclassno[i] <- which.max(ctpreds[i,1:K])
  }
  predclass <- names(table(data$class))[predclassno]
  conf <- table(predclass,data$class) # confusion matrix
  ## balanced accuracy
  acc <- rep(0,K)
  for (i in 1:length(rownames(conf))){
    for (j in 1:length(colnames(conf))){
      if (rownames(conf)[i] == colnames(conf)[j]){
        acc[j] <- conf[i,j] / colSums(conf)[j]
      }
    }
  }
  crit2 <- sum(acc) / K
  crit1sum <- sum(crit1)
## results
  result <- list(class = levels(data$class), trees = trees, ba=crit2, conf = conf, ninterp=crit1sum, acc=acci)
  class(result) <- "PrInDTMulev"
  result
}
#' @export
print.PrInDTMulev <- function(x, ...){
  # output of models
  K <- length(x$class)
  for (k in 1:K){
    cat("\n","tree for class",x$class[k],"\n")
    cat("balanced accuracy ",x$acc[k],"\n")
    print(x$trees[[k]])
  }
  # output for model characteristics
  cat("\n")
  cat("****************","\n")
  cat("Confusion matrix of combined predictions","\n")
  cat("****************")
  print(x$conf)
  cat("***********************************************","\n")
  cat("Criteria: balanced accuracy","\n")
  cat("***********************************************","\n")
  cat("             ",unname(x$ba),"\n")
  cat("***********************************************","\n")
  cat("Criteria: no. of non-interpretable trees","\n")
  cat("***********************************************","\n")
  cat("             ",unname(x$ninterp),"\n")
}
#' @export
plot.PrInDTMulev <- function(x, ...){
  # plot of models
  K <- length(x$class)
  for (k in 1:K){
    plot(x$trees[[k]],main=paste0("Tree for class = ",x$class[k]))
  }
}
#' @export
predict.PrInDTMulev <- function(object,classname,newdata,...){
K <- nlevels(newdata[,classname])
ctpreds <- matrix(0,nrow=dim(newdata)[1],ncol=K)
  for (k in 1:K)  {
    y <- t(array(unlist(stats::predict(object$trees[[k]],type="prob",newdata=newdata)),dim=c(2,dim(newdata)[1])))
    if (table(newdata[,classname])[k] < dim(newdata)[1]/2){
      ctpreds[,k] <- y[,2]
   } else {
      ctpreds[,k] <- y[,1]
   }
  }
  predclassno <- rep(0,dim(newdata)[1])
  for (j in 1:dim(newdata)[1]){
    predclassno[j] <- which.max(ctpreds[j,1:K])
  }
  predclass <- names(table(newdata[,classname]))[predclassno]
  result <-  as.factor(predclass)
  result
}