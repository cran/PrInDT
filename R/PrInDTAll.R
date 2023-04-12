#' Conditional inference tree (ctree) based on all observations
#'
#' @description ctree based on all observations. Interpretability is checked (see 'ctestv'); probability threshold can be specified.
#'
#' \strong{Reference}: Weihs, C., Buschfeld, S. 2021a. Combining Prediction and Interpretation in  Decision Trees (PrInDT) - 
#' a Linguistic Example. arXiv:2103.02336
#'
#' @usage PrInDTAll(datain, classname, ctestv=NA, conf.level=0.95, thres=0.5)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     {see function \code{\link{PrInDT}} for details.}\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical between 0 and 1); default = 0.95
#' @param thres Probability threshold for prediction of smaller class; default = 0.5
#'
#' @return
#' \describe{
#' \item{treeall}{ctree based on all observations}
#' \item{baAll}{balanced accuracy of 'treeall'}
#' \item{interpAll}{criterion of interpretability of 'treeall' (TRUE / FALSE)}
#' \item{confAll}{confusion matrix of 'treeall'}
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
PrInDTAll <- function(datain,classname,ctestv=NA,conf.level=0.95,thres=0.5){
  ## input check
  if ( typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character"))
      || !(0 <= conf.level & conf.level <= 1) || !(0 <= thres & thres <= 1) ){
    stop("irregular input")
  }
  ##
  data <- datain
  names(data)[names(data)==classname] <- "class"
  n_class1 <- table(data$class)[1] # no\item{of elements of larger class 1
  n_class2 <- table(data$class)[2] # no\item{of elements of smaller class 2
  if (n_class1 < n_class2){
    # relevel of classes if smaller class first
    data$class <- stats::relevel(data$class, levels(data$class)[2]) # larger class now first
    n_class1 <- table(data$class)[1] # no\item{of elements of larger class
    n_class2 <- table(data$class)[2] # no\item{of elements of smaller class
  }
  ## reordering of classes: smaller class first
  if (n_class1 > n_class2){
    order_class <- order(as.numeric(data$class),decreasing=TRUE)
  } else {
    order_class <- order(as.numeric(data$class))
  }
  data <- data[order_class,] # data now reordered: smaller class first
  ######
  ## model with all observations
  ######
  ct <- party::ctree(class ~ ., data = data, control = party::ctree_control(mincriterion=conf.level))
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
###
  result <- list(treeAll = ct, baAll=crit2, interpAll=crit1, confAll = conf)
  class(result) <- "PrInDTAll"
  result
}