#' Interdependent estimation for classification
#'
#' @description The function SimCPrInDT applies interdependent estimation (endogenous case) for finding an optimal model for relationships between the two-class factor variables
#' specified as column indices of 'data' in the vector 'inddep' and all other factor and numerical variables in the data frame 'data' by means of
#' 'N' repetitions of random subsampling with percentages 'percl' for the large classes and 'percs' for the small classes. One percentage of observations 
#' for each dependent variable has to be specified for the larger and the smaller class. For example, for three dependent variables, 'percl' consists of 
#' three percentages specified in the order in which the dependent variables appear in 'inddep'.\cr
#' The dependent variables have to be specified as dummies, i.e. as 'property absent' (value 0) or 'property present' (value 1).\cr
#' The optimization citerion is the balanced accuracy on the full sample. \cr
#' In an additional step, the mean balanced accuracy over all class variables is optimized (joint optimization).\cr
#' The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'. \cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage SimCPrInDT(data,ctestv=NA,inddep,percl,percs,N=99,M,psize,conf.level=0.95,
#'                                   minsplit=NA,minbucket=NA)
#'
#' @param data Input data frame with class factor variables and the influential variables,\cr
#'    which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param inddep indices of dependent variables
#' @param percl list of undersampling percentages of larger class (numerical, > 0 and <= 1): one per dependent class variable in the same order as in 'inddep'
#' @param percs list of undersampling percentage of smaller class (numerical, > 0 and <= 1); one per dependent class variable in the same order as in 'inddep'
#' @param N no. of repetitions of subsampling of observations (integer > 0); default = 99
#' @param M no. of repetitions of subsampling of predictors
#' @param psize no. of predictors in the subsamples of the predictors
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#'
#' @return
#' \describe{
#'   \item{models1}{Best trees at stage 1} 
#'   \item{models2}{Best trees at stage 2} 
#'   \item{models3}{Best trees from mean maximization} 
#'   \item{classnames}{names of classification variables}
#'   \item{baAll}{balanced accuracies of best trees at both stages} 
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.6.2, for further information.
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{name} as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export SimCPrInDT
#'
#' @examples
#' data <- PrInDT::data_land # load data
#' dataclean <- data[,c(1:7,23:24,11:13,22,8:10)]  # only relevant features
#' inddep <- c(14:16) # dependent variables
#' dataland <- na.omit(dataclean)
#' ctestv <- NA
#' perc <- c(0.45,0.05,0.25)   # percentages of observations of larger class, 
#' # 1 per dependent class variable
#' perc2 <- c(0.75,0.95,0.75)  # percentages of observations of smaller class, 
#' # 1 per dependent class variable
#' outland <- SimCPrInDT(dataland,percl=perc,percs=perc2,N=9,inddep=inddep,M=2,psize=10)
#' outland
#' plot(outland)
#'
SimCPrInDT <- function(data,ctestv=NA,inddep,percl,percs,N=99,M,psize,conf.level=0.95,minsplit=NA,minbucket=NA){
## input check
if ( typeof(data) != "list" || !(typeof(ctestv) %in% c("logical", "character")) 
    || !all(0 < percl & percl <= 1) || !all(0 < percs & percs <= 1) 
    || !(typeof(N) %in% c("integer","double")) || !(typeof(M) %in% c("integer","double")) 
    || !(typeof(inddep) %in% c("integer","double")) 
    || !(typeof(psize) %in% c("integer","double")) || !(0 < conf.level & conf.level <= 1) || !(typeof(minsplit) %in% c("logical","double")) 
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
#
## SimulPrInDT: classification
#
message ("1st stage: full sample")
depnames <- colnames(data)[inddep]
anz <- length(inddep)
D <- dim(data)[2]
models <- list()
modmax <- list()
equal <- matrix(TRUE,ncol=anz,nrow=dim(data)[1])
ctpreds <- matrix(0,ncol=anz,nrow=dim(data)[1])
ba <- c(1:anz)*0
for (i in 1:anz){
  outland <- PrInDTAll(data,depnames[i],ctestv=ctestv,conf.level=conf.level)
  ctpredsVL <- predict(outland$treeAll,newdata=data)
  equal[,i] <- (as.integer(data[,9])-1) == ctpredsVL
  ctpreds[,i] <- ctpredsVL
  ba[i] <- outland$baAll
  models[[i]] <- outland$treeAll
}
# Joint assessment: 01, Hamming
#  acc01 <- rowSums(equal)
#  ba[(anz+2)] <- sum(acc01) / (3*dim(data)[1])
#  acc01[acc01 < anz] <- 0
#  acc01[acc01 == anz] <- 1
#  ba[(anz+1)] <- sum(acc01)
#  ba[(anz+1)] <- ba[(anz+1)] / dim(data)[1]
##
#cat("Balanced accuracies for interdependent Models on full sample","\n")
max01 <- ba
#print(max01)
#print(mean(ba[1:anz]))
modmax <- models
#plot(unlist(modmax[[1]]))
#plot(unlist(modmax[[2]]))
#plot(unlist(modmax[[3]]))
baP <- c(ba[1:anz],mean(ba[1:anz]))
modelsP <- modmax
###############
##### 2nd stage: individual 
message("\n","2nd stage: individual optimization")
set.seed(7654321)
dataMend <- data[,-inddep]  ## all exogenous variables used
for (i in psize){
message("Number of predictors: ",i)
for (n in 1:M) {
  message("repetition ",n)
##  datas <- dataland[sample(rownames(dataland))[1:ssize],]  ## alternative: varying number of observations
  datas <- cbind(dataMend[,sample(colnames(dataMend))[1:psize]],data[,inddep])
  equal <- matrix(TRUE,ncol=anz,nrow=dim(data)[1])
  ba <- c(1:anz)*0
  for (j in 1:anz){
    perc <- percl[j]
    perc2 <- percs[j]
#    message(depnames[j])
#  outland <- PrInDT(dataland[,indall],"French",ctestv=NA,N,perc,perc2,conf.level=cf)
    outland <- PrInDT(datas,depnames[j],ctestv=ctestv,N=N,percl=perc,percs=perc2,conf.level=conf.level,seedl=FALSE,minsplit=minsplit,minbucket=minbucket)
#  print(outland)
#  plot(outland$tree1st)
    class.pred <- predict(outland$tree1st,newdata=data)
    class.pred <- relevel(class.pred,ref=levels(data[,inddep[j]])[1])
    conti <- table(class.pred,data[,inddep[j]])
    ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
    ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
    ba[j] <- (ba1 + ba2)/2
    models[[j]] <- outland$tree1st
    if (ba[j] > max01[j]) { 
      max01[j] <- ba[j]
#      print(max01)
      modmax[[j]] <- models[[j]]
    }
  }
##
# Joint assessment: 01
#  acc01 <- rowSums(equal)
#  ba[(anz+2)] <- sum(acc01) / (3*dim(data)[1])
#  acc01[acc01 < 3] <- 0
#  acc01[acc01 == 3] <- 1
#  ba[(anz+1)] <- sum(acc01)
#  ba[(anz+1)] <- ba[(anz+1)] / dim(data)[1]
}
}
#print(max01[1:anz])
#print(sum(max01[1:anz]) / 3)
#plot(unlist(modmax[[1]]))
#plot(unlist(modmax[[2]]))
#plot(unlist(modmax[[3]]))
baI <- c(max01[1:anz],mean(max01[1:anz]))
modelsI <- modmax
##########################
###########  mean maximization
message("\n","3rd stage: joint optimization")
set.seed(7654321)
dataMend <- data[,-inddep]  ## all exogenous variables used
bam <- 0
for (i in psize){
message("Number of predictors: ",i)
for (n in 1:M) {
  message("repetition ",n)
##  datas <- dataland[sample(rownames(dataland))[1:ssize],]  ## alternative: varying number of observations
  datas <- cbind(dataMend[,sample(colnames(dataMend))[1:psize]],data[,inddep])
  equal <- matrix(TRUE,ncol=anz,nrow=dim(data)[1])
  ba <- c(1:anz)*0
  for (j in 1:anz){
    perc <- percl[j]
    perc2 <- percs[j]
#    message(depnames[j])
#  outland <- PrInDT(dataland[,indall],"French",ctestv=NA,N,perc,perc2,conf.level=cf)
    outland <- PrInDT(datas,depnames[j],ctestv=ctestv,N=N,percl=perc,percs=perc2,conf.level=conf.level,seedl=FALSE,minsplit=minsplit,minbucket=minbucket)
#  print(outland)
#  plot(outland$tree1st)
    class.pred <- predict(outland$tree1st,newdata=data)
    class.pred <- relevel(class.pred,ref=levels(data[,inddep[j]])[1])
    conti <- table(class.pred,data[,inddep[j]])
    ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
    ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
    ba[j] <- (ba1 + ba2)/2
    models[[j]] <- outland$tree1st
  }
  if (mean(ba) > bam) { 
    bam <- mean(ba)
    max01[1:anz] <- ba[1:anz]
#    print(c(max01[1:anz],bam))
    for (j in 1:anz){
      modmax[[j]] <- models[[j]]
    }
  }
}
}
baJ <- c(max01[1:anz],bam)
#print(baJ)
baAll <- rbind(baP,baI,baJ)
result <- list(models1=modelsP,models2=modelsI,models3=modmax,classnames=depnames,baAll=baAll) 
class(result) <- "SimCPrInDT"
result
}
#' @export
print.SimCPrInDT <- function(x,...){
## output
cat("Accuracies of best models","\n")
colnames(x$baAll) <- c(x$classnames,"mean")
rownames(x$baAll) <- c("full sample","indiv. opt.","joint opt.")
print(x$baAll)
cat("\n\n")
for (i in 1:length(x$classnames)){
  cat("Best model 2nd stage: ",x$classnames[i], "\n")
  print(x$models2[[i]])
  cat("\n")
}
}
#' @export
plot.SimCPrInDT <- function(x,...){
for (i in 1:length(x$classnames)){
  plot(x$models2[[i]],main=x$classnames[i])
}
}