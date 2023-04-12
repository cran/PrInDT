#' Conditional inference tree (ctree) for multiple classes on all observations
#'
#' @description ctree for more than 2 classes on all observations. Interpretability is checked (see 'ctestv').
#'
#' @usage PrInDTMulevAll(datain, classname, ctestv=NA, conf.level=0.95)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     {see function \code{\link{PrInDT}} for details.}\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical between 0 and 1)\cr
#'     (default = 0.95)
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
#' frame of the function.
#'
#' @exportS3Method print PrInDTMulevAll
#' @exportS3Method plot PrInDTMulevAll
#' @export PrInDTMulevAll
#'
#' @examples
#' datastrat <- PrInDT::data_zero
#' data <- na.omit(datastrat)
#' ctestv <- rbind('ETH == {C2a,C1a}', 'MLU == {1, 3}')
#' data$rel[data$ETH %in% c("C1a","C1b","C1c") & data$real == "zero"] <- "zero1"
#' data$rel[data$ETH %in% c("C2a","C2b","C2c") & data$real == "zero"] <- "zero2"
#' data$rel[data$real == "realized"] <- "real"
#' data$rel <- as.factor(data$rel) # rel is new class variable
#' data$real <- NULL # remove old class variable
#' conf.level <- 0.99 # 1 - significance level (mincriterion) in ctree
#' outAll <- PrInDTMulevAll(data,"rel",ctestv,conf.level) 
#' outAll # print model based on all observations
#' plot(outAll)
#'
#' @importFrom stats relevel predict
#' @importFrom party ctree ctree_control
#'
PrInDTMulevAll <- function(datain,classname,ctestv=NA,conf.level=0.95){
  ## input check
  if (typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || 
      !(0 <= conf.level & conf.level <= 1)){
    stop("irregular input")
  }
  data <- datain
  names(data)[names(data)==classname] <- "class"
  K <- length(table(data$class))
  ######
  ## model with all observations
  ######
  ct <- party::ctree(class ~ ., data = data, control = party::ctree_control(mincriterion=conf.level))
  crit1 <- FALSE
  if (is.na(ctestv[1]) == FALSE) {
    crit1 <- FindSubstr(ct,ctestv) # call of the above function for overall tree
  }
  ctpreds <- stats::predict(ct) # predictive fit
  conf <- table(ctpreds,data$class)  # confusion matrix
  ## balanced accuracy
  acc <- rep(0,K)
  for (j in 1:K){
        acc[j] <- conf[j,j] / colSums(conf)[j]
  }
  crit2 <- sum(acc) / K
  result <- list(treeAll = ct, baAll=crit2, confAll = conf, interpAll=crit1)
  class(result) <- "PrInDTMulevAll"
  result
}

####
## print function
####
print.PrInDTMulevAll <- function(x, ...){
  # output for model on all observations
  cat("******************************","\n")
  cat("Tree built on all observations","\n")
  cat("******************************","\n")
  print(x$treeAll) # print of tree
  cat("\n")
  cat("****************","\n")
  cat("Confusion matrix","\n")
  cat("****************")
  print(x$confAll)
  cat("***********************************************","\n")
  cat("Criteria: balanced accuracy","\n")
  cat("***********************************************","\n")
  cat("             ",unname(x$baAll),"\n")
  cat("***********************************************","\n")
  cat("Criteria: interpretable tree","\n")
  cat("***********************************************","\n")
  cat("             ",unname(!(x$interpAll)),"\n")
}

####
## plot function
####
plot.PrInDTMulevAll <- function(x, ...){
  plot(x$treeAll,main="Multi-class problem: Tree on full sample") # plot of tree
}