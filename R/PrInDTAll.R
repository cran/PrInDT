#' Conditional inference tree (ctree) based on all observations
#'
#' @description ctree based on all observations in 'datain'. Interpretability is checked (see 'ctestv'); probability threshold can be specified.
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.
#'
#' \strong{Reference} \cr
#' Weihs, C., Buschfeld, S. 2021a. Combining Prediction and Interpretation in  Decision Trees (PrInDT) - 
#' a Linguistic Example. arXiv:2103.02336
#'
#' In the case of repeated measurements ('indrep=1'), the values of the substructure variable have to be given in 'repvar'.
#' Only one value of 'classname' is allowed for each value of 'repvar'.
#' If for a value of 'repvar' the percentage 'thr' of the observed occurence of a value of 'classname' is not reached by the number of predictions of the value of 'classname', a misclassification is detected.
#'
#' @usage PrInDTAll(datain, classname, ctestv=NA, conf.level=0.95, thres=0.5,
#'                  minsplit=NA,minbucket=NA,repvar=NA,indrep=0,thr=0.5)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     (see function \code{\link{PrInDT}} for details.)\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1); default = 0.95
#' @param thres Probability threshold for prediction of smaller class (numerical, >= 0 and < 1); default = 0.5
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#' @param repvar Values of variable defining the substructure in the case of repeated measurements; default=NA
#' @param indrep Indicator of repeated measurements ('indrep=1'); default = 0
#' @param thr threshold for element classification: minimum percentage of correct class entries; default = 0.5
#'
#' @return
#' \describe{
#' \item{treeall}{ctree based on all observations}
#' \item{baAll}{balanced accuracy of 'treeall'}
#' \item{interpAll}{criterion of interpretability of 'treeall' (TRUE / FALSE)}
#' \item{confAll}{confusion matrix of 'treeall'}
#' \item{acc1AE}{Accuracy of full sample tree on Elements of large class} 
#' \item{acc2AE}{Accuracy of full sample tree on Elements of small class}  
#' \item{bamaxAE}{balanced accuracy of full sample tree on Elements} 
#' \item{namA1}{Names of misclassified Elements by full sample tree of large class}
#' \item{namA2}{Names of misclassified Elements by full sample tree of small class}
#' \item{lablarge}{Label of large class} 
#' \item{labsmall}{Label of small class}
#' \item{thr}{Threshold for repeated measurements} 
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.\cr
#'
#' @export PrInDTAll
#'
#' @examples
#' datastrat <- PrInDT::data_zero
#' data <- na.omit(datastrat)
#' ctestv <- rbind('ETH == {C2a,C1a}','MLU == {1, 3}')
#' conf.level <- 0.99 # 1 - significance level (mincriterion) in ctree
#' outAll <- PrInDTAll(data,"real",ctestv,conf.level) 
#' print(outAll) # print model based on all observations
#' plot(outAll) # plot model based on all observations
#'
#' @importFrom stats relevel predict
#' @importFrom party ctree ctree_control
#'
#'
PrInDTAll <- function(datain,classname,ctestv=NA,conf.level=0.95,thres=0.5,minsplit=NA,minbucket=NA,repvar=NA,indrep=0,thr=0.5){
  ## input check
  if ( typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character"))
      || !(0 < conf.level & conf.level <= 1) || !(0 <= thres & thres < 1) || !(typeof(minsplit) %in% c("logical","double")) 
      || !(typeof(minbucket) %in% c("logical", "double")) || !(typeof(repvar) %in% c("logical","integer","character"))
      || !(typeof(indrep) %in% c("integer","double")) || typeof(thr) != "double" ){
    stop("irregular input")
  }  
  if (indrep == 1 & all(is.na(repvar)) == TRUE){
    message("\n","Note: repeated measurement analysis disabled: no repeated measurement variable specified","\n")
    indrep <- 0
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
  ##
  data <- datain
  repv <- repvar
  names(data)[names(data)==classname] <- "class"
  n_class1 <- table(data$class)[1] # no\item{of elements of larger class 1
  n_class2 <- table(data$class)[2] # no\item{of elements of smaller class 2
  if (min(n_class1,n_class2) == 0){
   stop("irregular input: only 1 class")
  }
  n <- n_class1 + n_class2
  if (n_class1 < n_class2){
    # relevel of classes if smaller class first
    data$class <- stats::relevel(data$class, levels(data$class)[2]) # larger class now first
    n_class1 <- table(data$class)[1] # no of elements of larger class
    n_class2 <- table(data$class)[2] # no of elements of smaller class
  }
  ## reordering of classes: smaller class first
  if (n_class1 > n_class2){
    order_class <- order(as.numeric(data$class),decreasing=TRUE)
  } else {
    order_class <- order(as.numeric(data$class))
  }
  data <- data[order_class,] # data now reordered: smaller class first
  if (indrep == 1) {repv <- repv[order_class]}
  ######
  ## model with all observations
  ######
#if (n < 60){
#print(n)
#      ct <- party::ctree(class ~ ., data = data,control = party::ctree_control(mincriterion=conf.level,minsplit=5,minbucket=3))  ## TEST TEST
#    } else {
      ct <- party::ctree(class ~ ., data = data,control = party::ctree_control(mincriterion=conf.level,minsplit=minsplit,minbucket=minbucket)) 
#    }
#  ct <- party::ctree(class ~ ., data = data, control = party::ctree_control(mincriterion=conf.level))
  if (thres != 0.5){
    predsprob <- stats::predict(ct,type="prob")
    ctpreds <- as.factor(sapply( 1:dim(data)[1],
                                 function(x)
                                   ifelse(predsprob[[x]][2] > thres,
                                          levels(data$class)[2], levels(data$class)[1] )))
  } else {
    ctpreds <- stats::predict(ct) # predictions for all observations
}
  conf <- table(ctpreds, data$class)  # confusion matrix calculation
  ## balanced accuracy
  acc1 <- conf[1,1] / n_class1
  acc2 <- conf[2,2] / n_class2
  crit2 <- (acc1 + acc2)/2
  crit1 <- "FALSE"
  if (is.na(ctestv[1]) == FALSE) {
    crit1 <- FindSubstr(ct,ctestv) # call of the above function for overall tree
  }
lablarge <- names(table(data$class))[1] # class 1 = large
labsmall <- names(table(data$class))[2] # class 2 = small 
acc1AE <- 0
acc2AE <- 0
bamaxAE <- 0
namA1 <- list()
namA2 <- list()
if (indrep == 1){
  pred <- predict(ct,newdata=data)
  ch <- table(repv,pred)
  conti <- cbind(table(repv,data$class),ch)
#print(conti)
  no1 <- sum(conti[,1] > 0)
  no2 <- sum(conti[,2] > 0)
  if ((no1 + no2) > dim(conti)[1]){
#    cat("repeated measurement variable has more than one class per value")
    acc2AE <- -1
    acc1AE <- -1 
  } else {
    for (i in 1:dim(conti)[1]){
      if (conti[i,2] > 0 & conti[i,4] <= (conti[i,2]*thr)){ 
        namA2 <- c(namA2,rownames(ch)[i])
      }
      if (conti[i,1] > 0 & conti[i,3] < (conti[i,1]*thr)){ 
        namA1 <- c(namA1,rownames(ch)[i])
      }
    }
    acc2AE <- 1 - length(namA2) / no2
    acc1AE <- 1 - length(namA1) / no1  
    bamaxAE <- (acc1AE + acc2AE)/2
    if (length(namA1) == 0) {namA1 <- "-"}
    if (length(namA2) == 0) {namA2 <- "-"} 
  }
}
###
  result <- list(treeAll = ct, baAll=crit2, interpAll=crit1, confAll = conf,acc1AE=acc1AE,acc2AE=acc2AE,bamaxAE=bamaxAE,
                 namA1=namA1,namA2=namA2,indrep=indrep,lablarge=lablarge,labsmall=labsmall,thr=thr)
  class(result) <- "PrInDTAll"
  result
}