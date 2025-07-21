#' Conditional inference trees (ctrees) based on consecutive parts of the full sample
#'
#' @description ctrees based on the full sample of the smaller class and consecutive parts of the larger class of the nesting variable 'nesvar'. 
#' The variable 'nesvar' has to be part of the data frame 'datain'.\cr   
#' Interpretability is checked (see 'ctestv'); probability threshold can be specified.
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' \strong{Reference}\cr Weihs, C., Buschfeld, S. 2021b. NesPrInDT: Nested undersampling in PrInDT. 
#'	arXiv:2103.14931
#'
#' @usage PrInDTAllparts(datain, classname, ctestv=NA, conf.level=0.95, thres=0.5,
#'        nesvar, divt,minsplit=NA,minbucket=NA)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     (see function \code{\link{PrInDT}} for details.)\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1); default = 0.95
#' @param thres Probability threshold for prediction of smaller class (numerical, >= 0 and < 1); default = 0.5
#' @param nesvar Name of nesting variable (character)
#' @param divt Number of parts of nesting variable nesvar for which models should be determined individually
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#'
#' @return
#' \describe{
#' \item{baAll}{balanced accuracy of tree on full sample}
#' \item{nesvar}{name of nesting variable}
#' \item{divt}{number of consecutive parts of the sample}
#' \item{badiv}{balanced accuracy of trees on 'divt' consecutive parts of the sample}
#' }
#'
#' @details
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } where 'name' is the output data 
#' frame of the function.
#'
#' @export PrInDTAllparts
#'
#' @examples
#' data <- PrInDT::data_speaker
#' data <- na.omit(data)
#' nesvar <- "SPEAKER"
#' outNesAll <- PrInDTAllparts(data,"class",ctestv=NA,conf.level=0.95,thres=0.5,nesvar,divt=8)
#' outNesAll
#'
PrInDTAllparts <- function(datain,classname,ctestv=NA,conf.level=0.95,thres=0.5,nesvar,divt,minsplit=NA,minbucket=NA){
  ## input check
  if (typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || 
      !(0 < conf.level & conf.level <= 1) || !(0 <= thres & thres < 1) ||
      typeof(nesvar) != "character" || divt <= 0 || !(typeof(minsplit) %in% c("logical","double")) || 
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
  ##
  data <- datain
  names(data)[names(data)==classname] <- "class"
  n_class1 <- table(data$class)[1] # no. of elements of larger class 1
  n_class2 <- table(data$class)[2] # no. of elements of smaller class 2
  if (min(n_class1,n_class2) == 0){
   stop("irregular input: only 1 class")
  }
  if (n_class1 < n_class2){
    # relevel of classes if smaller class first
    data$class <- stats::relevel(data$class, levels(data$class)[2]) # larger class now first level
    n_class1 <- table(data$class)[1] # no. of elements of larger class
    n_class2 <- table(data$class)[2] # no. of elements of smaller class
  }
  ## full sample
  out <- PrInDTAll(data,"class",minsplit=minsplit,minbucket=minbucket)
  baAll <- out$baAll
  ## analyses of parts
  names(data)[names(data)==nesvar] <- "NEST"
  n_classind1 <- table(data$NEST)[1] # no. of elements of class 1 in nesvar
  n_classind2 <- table(data$NEST)[2] # no. of elements of class 2 in nesvar
  if (n_classind1 < n_classind2){
    # relevel of classes if smaller class first
    data$NEST <- stats::relevel(data$NEST, levels(data$NEST)[2]) # larger class now first level
    n_classind1 <- table(data$NEST)[1] # no. of elements of larger class
    n_classind2 <- table(data$NEST)[2] # no. of elements of smaller class
  }
  ## reordering of classes: smaller class first
  if (n_classind1 > n_classind2){
    order_class2 <- order(as.numeric(data$NEST),decreasing=TRUE)
  } else {
    order_class2 <- order(as.numeric(data$NEST))
  }
  x <- data[order_class2,] # data now reordered: smaller class of nesvar first
  badiv <- vector("numeric",length=divt)
  for (i in 0:(divt-1)){
    datadiv <- rbind (x[1:n_classind2,], x[(n_classind2 + (i*n_classind1/divt)+1):(n_classind2 + (i+1)*n_classind1/divt),])
    outdiv <- PrInDTAll(datadiv,"class",minsplit=minsplit,minbucket=minbucket)
    badiv[i+1] <- outdiv$baAll
  }
  result <- list(baAll = baAll, nesvar = nesvar, divt = divt, badiv = badiv)
  class(result) <- "PrInDTAllparts"
  result 
}
#' @export
print.PrInDTAllparts <- function(x,...){
  cat("\n","Balanced accuracy of tree on full sample: ",unname(x$baAll),"\n")
  cat("\n","Consecutive parts of variable ",x$nesvar,"\n")
  cat("Balanced accuracy of trees on", x$divt, "consecutive parts of the sample","\n")
  for (i in 1:x$divt){
    ba <- x$badiv[i]
    cat("\n","**** part",i," ",unname(ba))
  }
  cat("\n\n")
}
#' @export
plot.PrInDTAllparts <- function(x,...){
  cat("\n","No plots available","\n")
}

