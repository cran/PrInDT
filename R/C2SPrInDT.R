#' Two-stage estimation for classification
#'
#' @description The function C2SPrInDT applies two-stage estimation for finding an optimal model for relationships between the two-class factor variables
#' specified as column indices of 'datain' in the vector 'inddep' and all other factor and numerical variables in the data frame 'datain' by means of
#' 'N' repetitions of random subsampling with percentages 'percl' for the large classes and 'percs' for the small classes. One percentage of observations 
#' for each dependent variable has to be specified for the larger and the smaller class. For example, for three dependent variables, 'percl' consists of 
#' three percentages specified in the order in which the dependent variables appear in 'inddep'.\cr
#' The dependent variables have to be specified as dummies, i.e. as 0 for 'property absent' or 1 for 'property present' for a certain property the dependent variable is representing.\cr
#' The indices of the predictors relevant at 1st stage modeling can be specified in the vector 'indind'. If 'indind' is not specied, then all variables in 'datain' 
#' which are not specified in 'inddep' are used as predictors at 1st stage. At 2nd stage, all such variables are used as predictors anyway.\cr
#' The optimization citerion is the balanced accuracy on the full sample. The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'. \cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage C2SPrInDT(datain,ctestv=NA,conf.level=0.95,percl,percs,N=99,indind=NA,inddep,
#'                                  minsplit=NA,minbucket=NA)
#'
#' @param datain Input data frame with class factor variables and the influential variables,\cr
#'    which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param percl list of undersampling percentages of larger class (numerical, > 0 and <= 1): one per dependent class variable in the same order as in 'inddep'
#' @param percs list of undersampling percentage of smaller class (numerical, > 0 and <= 1); one per dependent class variable in the same order as in 'inddep'
#' @param N no. of repetitions (integer > 0); default = 99
#' @param indind indices of independent variables at stage 1; default = NA (means all independent variables used)
#' @param inddep indices of dependent variables
#' @param minsplit Minimum number of elements in a node to be splitted; default = 20
#' @param minbucket Minimum number of elements in a node; default = 7
#'
#' @return
#' \describe{
#'   \item{models1}{Best trees at stage 1} 
#'   \item{models2}{Best trees at stage 2} 
#'   \item{classnames}{names of classification variables}
#'   \item{baAll}{balanced accuracies of best trees at both stages} 
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.6.1, for further information.
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{name} as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export C2SPrInDT
#'
#' @examples
#' data <- PrInDT::data_land # load data
#' dataclean <- data[,c(1:7,23:24,11:13,22,8:10)]  # only relevant features
#' indind <- c(1:9) # original predictors
#' inddep <- c(14:16) # dependent variables
#' dataland <- na.omit(dataclean)
#' ctestv <- NA
#' perc <- c(0.45,0.05,0.25)   # percentages of observations of larger class, 
#' # 1 per dependent class variable
#' perc2 <- c(0.75,0.95,0.75)  # percentages of observations of smaller class, 
#' # 1 per dependent class variable
#' outland <- C2SPrInDT(dataland,percl=perc,percs=perc2,N=19,indind=indind,inddep=inddep)
#' outland
#' plot(outland)
#'
C2SPrInDT <- function(datain,ctestv=NA,conf.level=0.95,percl,percs,N=99,indind=NA,inddep,minsplit=NA,minbucket=NA){
## input check
if ( typeof(datain) != "list" || !(typeof(ctestv) %in% c("logical", "character"))
    || !all(0 < percl & percl <= 1) || !all(0 < percs & percs <= 1) 
    || !(typeof(N) %in% c("integer","double")) 
    || !(typeof(inddep) %in% c("integer","double")) || !(typeof(indind) %in% c("logical","integer","double"))
    || !(0 < conf.level & conf.level <= 1) || !(typeof(minsplit) %in% c("logical","double")) 
    || !(typeof(minbucket) %in% c("logical", "double"))  ) {
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
data <- datain
classnames <- colnames(data)[inddep]
anz <- length(inddep)
message ("1st stage")
models <- list()
equal <- matrix(TRUE,ncol=anz,nrow=dim(data)[1])
ctpreds <- matrix(0,ncol=anz,nrow=dim(data)[1])
ba <- c(1:(anz+2))*0
for (i in 1:anz){
    x <- colnames(datain)[inddep[i]]
    levels(data[,x]) <- c(levels(data[,x]),x)
    if (table(data[,x])[2] > dim(data)[1]/2){
      data[datain[,x] == 1,x] <- x
    } else {
      data[datain[,x] != 1,x] <- x
    } 
    data[,x] <- droplevels(data[,x])
#    message("\n")
    message(paste0("  ",x))
#  message(colnames(data)[inddep[i]])
  if (all(is.na(indind)) == TRUE){
    indall <- c(1:(dim(data)[2]))
    indall <- setdiff(indall,inddep)
    indall <- c(inddep[i],indall)
  } else {
    indall <- c(inddep[i],indind) # externals + ith dependent
  }
  perc <- percl[i]
  perc2 <- percs[i]
  outland <- PrInDT(data[,indall],classnames[i],ctestv=ctestv,N,perc,perc2,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
#  print(outland)  ##
#  plot(outland$tree1st)
  ctpredsVL <- predict(outland$tree1st,newdata=data)
  equal[,i] <- data[,inddep[i]] == ctpredsVL
#  equal[,i] <- (as.integer(data[,inddep[i]])-1) == ctpredsVL
  ctpreds[,i] <- ctpredsVL
  ba[i] <- outland$ba1st[3]
  models <- c(models,outland$tree1st)
}
  acc01 <- rowSums(equal)
  ba[anz+2] <- sum(acc01) / (anz*dim(data)[1])  ## indices generalized
  acc01[acc01 < anz] <- 0
  acc01[acc01 == anz] <- 1
  ba[anz+1] <- sum(acc01)
  ba[anz+1] <- ba[anz+1] / dim(data)[1]
#print(ba)
##
## 2nd stage
#
message("2nd stage")
dataP <- data
dataP <- cbind(dataP,ctpreds-1)  ## ????????????   change of coding: 2-1 -> 1-0
dim(dataP)
baP <- ba
modelsP <- models
ctpredsP <- ctpreds
models <- list()
equal <- matrix(TRUE,ncol=anz,nrow=dim(data)[1])
ctpreds <- matrix(0,ncol=anz,nrow=dim(data)[1])
ba <- c(1:(anz+2))*0
#dataP[,dim(data)[2]+anz] <- 1 - dataP[,dim(data)[2]+anz]  ## generalized above
for (i in 1:anz){
  dataP[,dim(data)[2]+i] <- as.factor(dataP[,dim(data)[2]+i])
  names(dataP)[dim(data)[2]+i] <- paste0(classnames[i],"_2S")
}
for (i in 1:anz){
  message(paste0("  ",classnames[i]))
  indall <- c(1:(dim(dataP)[2]))
  j <- dim(data)[2]+i
  indall <- indall[-j] # externals + predictions - endogenous
  j <- inddep[-i]
  indall <- indall[-j] # externals + predictions - endogenous
  perc <- percl[i]
  perc2 <- percs[i]
  outland <- PrInDT(dataP[,indall],classnames[i],ctestv=ctestv,N,perc,perc2,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
#  print(outland)  ##
#  plot(outland$tree1st)
  ctpredsVL <- predict(outland$tree1st,newdata=dataP)
  equal[,i] <- data[,inddep[i]] == ctpredsVL
#  equal[,i] <- (as.integer(dataP[,inddep[i]])-1) == ctpredsVL
  ctpreds[,i] <- ctpredsVL
  ba[i] <- outland$ba1st[3]
  models <- c(models,outland$tree1st)
}
  acc01 <- rowSums(equal)
  ba[anz+2] <- sum(acc01) / (anz*dim(data)[1])  ## indices generalized
  acc01[acc01 < anz] <- 0
  acc01[acc01 == anz] <- 1
  ba[anz+1] <- sum(acc01)
  ba[anz+1] <- ba[anz+1] / dim(data)[1]
baP <- rbind(baP,ba)
result <- list(models1=modelsP,models2=models,classnames=classnames,baAll=baP) 
class(result) <- "C2SPrInDT"
result
}
##
#' @export
print.C2SPrInDT <- function(x,...){
## output
cat("Accuracies of best models","\n")
colnames(x$baAll) <- c(x$classnames,"01-acc.","Hamming")
rownames(x$baAll) <- c("1st stage","2nd stage")
print(x$baAll)
cat("\n\n")
for (i in 1:length(x$classnames)){
  cat("Best model 2nd stage: ",x$classnames[i], "\n")
  print(x$models2[[i]])
  cat("\n")
}
}
#' @export
plot.C2SPrInDT <- function(x,...){
for (i in 1:length(x$classnames)){
  plot(x$models2[[i]],main=x$classnames[i])
}
}