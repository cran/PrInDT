#' Regression tree based on all observations
#'
#' @description Regression tree based on the full sample; interpretability is checked (see 'ctestv').\cr
#' The relationship between the target variable 'regname' and all other factor and numerical variables
#' in the data frame 'datain' is modeled based on all observations. \cr 
#'
#' @usage PrInDTregAll(datain, regname, ctestv=NA, conf.level=0.95)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param regname name of regressand variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     {see function \code{\link{PrInDT}} for details.}\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#'
#' @return
#' \describe{
#' \item{treeall}{tree based on all observations}
#' \item{R2All}{goodness of fit of 'treeall' based on all observations}
#' \item{interpAll}{criterion of interpretability of 'treeall' (TRUE / FALSE)}
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.
#'
#' @exportS3Method print PrInDTregAll
#' @exportS3Method plot PrInDTregAll
#' @export PrInDTregAll 
#'
#' @examples
#' data <- PrInDT::data_vowel
#' data <- na.omit(data)
#' ctestv <- 'vowel_maximum_pitch <= 320'
#' outreg <- PrInDTregAll(data,"target",ctestv)
#' outreg
#' plot(outreg)
#'
#' @importFrom party ctree ctree_control
#' @importFrom stats formula
#' @importFrom graphics plot
#'
PrInDTregAll <- function(datain,regname,ctestv=NA,conf.level=0.95){
  ## input check
  if (typeof(datain) != "list" || typeof(regname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) ||
      !(0 < conf.level & conf.level <= 1)){
    stop("irregular input")
  }
####
  data <- datain
  ct <- party::ctree(formula(paste(regname, "~ .")), data = data,control = party::ctree_control(mincriterion=conf.level))
  ctpreds <- predict(ct)
  R2 <- 1 - sum((data[,regname] - ctpreds)^2) / sum((data[,regname]-mean(data[,regname]))^2)
  crit1 <- "FALSE"
  if (all(is.na(ctestv)) == FALSE) {
    crit1 <- FindSubstr(ct,ctestv) # call of the above function for overall tree
  }
  result <- list(treeAll = ct, R2All = R2, interpAll = crit1)
  class(result) <- "PrInDTregAll"
  result
}
##
## print function
print.PrInDTregAll <- function(x,...){
  cat("\n","Regression tree on full data set","\n")
  print(x$treeAll)
  cat("\n","Interpretability:",unname(!as.logical(x$interpAll)))
  cat("\n","Goodness of Fit R2:",x$R2All,"\n\n")
}
##
## plot function
plot.PrInDTregAll <- function(x,...){
  plot(x$treeAll,main="Regression tree on full data set")
}