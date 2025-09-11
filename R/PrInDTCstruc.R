#' Structured subsampling for classification
#'
#' @description The function PrInDTCstruc applies structured subsampling for finding an optimal subsample to model
#' the relationship between the two-class factor variable 'classname' and all other factor and numerical variables
#' in the data frame 'datain' by means of 'N' repetitions of subsampling from a substructure and 'Ni' repetitions of subsampling from the predictors.  
#' The optimization citerion is the balanced accuracy on the validation sample 'valdat' (default is the full input sample 'datain'). 
#' Other criteria are possible (cf. parameter description of 'crit'). The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The substructure of the observations used for subsampling in modelling is specified by the list 'Struc' which consists of the factor variable 'name' representing the substructure,
#' the name 'check' of the factor variable with the information about the categories of the substructure, and the matrix 'labs' which specifies the values of 'check'
#' corresponding to two categories in its rows, i.e. in 'labs[1,]' and 'labs[2,]'. The names of the categories have to be specified by \code{rownames(labs)}.\cr
#' See parameter description of 'Struc' for its specification for 'vers="b"' and 'indrep > 0'.\cr
#' The number of predictors 'Pit' to be included in the model and the number of elements of the substructure 'Eit' have to be specified (lists allowed), and 
#' undersampling of the categories of 'classname' can be controlled by 'undersamp=TRUE/FALSE'.\cr
#' In structured subsampling, 'N' repetitions of subsampling of the variable 'name' with 'Eit' different elements of each category in 'check' are realized. If 'Eit' is a list, each entry is employed individually. If 'Eit' is larger than the maximum available number of elements with a certain value of 'check', the maximum possible number of elements is used.\cr
#' Four different versions of structured subsampling exist: \cr
#' a) just of the elements in the substructure (possibly with additional undersampling) with parameters 'N' and 'Eit',\cr
#' b) just of the predictors with parameters 'Ni' and 'Pit', \cr
#' c) of the predictors and for each subset of predictors subsampling of the elements of the substructure (possibly with additional undersampling) 
#' with parameters 'N', 'Ni', 'Eit', and 'Pit', and\cr
#' d) of the elements of the substructure (possibly with additional undersampling) and for each of these subsets subsampling of the predictors 
#' with the same parameters as version c).\cr
#' Sampling of the elements of the substructure can be influenced by using weights of the elements ('weights=TRUE') according to the number of appearances of the smaller 
#' class of 'classnames'. This way, elements with more realisations in the smaller class are preferred in modelling.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr\cr
#' The parameter 'indrep' indicates repeated measurement situations ('indrep > 0', only implemented for 'crit="ba"'); default = 0.\cr
#' Repeated measurements are multiple measurements of the same variable taken on the same subjects (or objects) either under different conditions or over two or more time periods.\cr
#' The variable with the repeatedly observed subjects (or objects) has to be specified by 'name' in 'Struc'. Only one value of 'classname' is allowed for each value of 'Struc$name'.
#' In case of 'indrep > 0' it is automatically assumed that the same number of 
#' subjects (or objects) of the two classes under study have to be used for model building. Possible such numbers can be specified by 'Eit'.\cr
#' For indrep=1, models are optimized for individual observations and the classes of the substructure elements are subsequently determined via the parameter 'thr'.\cr
#' For indrep=2, models are optimized so that the criterion "ba" is maximal for the substructure elements taking 'thr' into account.\cr
#' For indrep > 0 and 'valdat' unequal 'datain', the variable 'repvar' defines the substructure for 'valdat', length = dim(valdat)[1] necessary; default = NA'     
#'
#' @usage PrInDTCstruc(datain,classname,ctestv=NA,Struc=NA,vers="d",weight=FALSE,
#'                 Eit=NA,Pit=NA,N=99,Ni=99,undersamp=TRUE,crit="ba",ktest=0,
#'                 stest=integer(length=0),conf.level=0.95,indrep=0,
#'                 minsplit=NA,minbucket=NA,repvar=NA,valdat=datain,thr=0.5)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param classname Name of class variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param Struc = list(name,check,labs), cf. description for explanations; Struc not needed for vers="b"; for indrep=1, Struc = list(name);
#' @param vers Version of structured subsampling: "a", "b", "c", "d", cf. description; \cr
#'     default = "d"
#' @param weight Weights to be used for subsampling of elements of substructure (logical, TRUE or FALSE);
#'     default = FALSE
#' @param Eit List of number of elements of substructure (integers); \cr
#'     default = c((Cl-4):Cl), Cl = maximum number of elements in both categories 
#' @param Pit List of number of predictors (integers)\cr
#'     default = c(max(1,(D-3)):D), D = maximum number of predictors
#' @param N Number of repetitions of subsampling from substructure (integer)\cr
#'      default = 0 for vers="b", = 99 otherwise; if vers="b", any input is overwritten by the default
#' @param Ni Number of repetitions of subsampling from predictors\cr
#'      default = 0 for vers="a", = 99 for vers="b", = N otherwise; if vers="a", any input is overwritten by the default
#' @param undersamp Undersampling of categories of 'classname' to be used ((logical, TRUE or FALSE)\cr
#'     default = TRUE; if indrep=1 or vers="b", undersamp = FALSE
#' @param crit Optimisation criterion: "ba" for balanced accuracy, "bat" for balanced accuracy on test sets, "ta" for test accuracy, 
#'     "tal" for test accuracy of continuing parts of length 'ktest' in substructure elements 'stest';\cr
#'     default = "ba"
#' @param ktest Length of continuing parts to be tested (for crit="tal");\cr
#'     default = 0
#' @param stest Part of substructure to be tested (for crit="tal")(integer vector);\cr
#'     default = integer vector of length 0
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param indrep Indicator for repeated measurements, i.e. more than one observation with the same class for each element;\cr 
#'     for indrep=1, Struc=list(name) only; default = 0 
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#' @param valdat validation data; default = datain
#' @param repvar Values of variable defining the substructure corresponding to valdat in the case of repeated measurements, length = dim(valdat)[1] necessary; default = NA
#' @param thr threshold for element classfication: minimum percentage of correct class entries; default = 0.5
#'
#' @return
#' \describe{
#'   \item{modbest}{Best tree} 
#'   \item{interp}{Number of interpretable trees, overall number of trees} 
#'   \item{dmax}{Number of predictors in training set for best tree} 
#'   \item{ntmax}{Size of training set for best tree} 
#'   \item{acc1}{Accuracy of best tree on large class} 
#'   \item{acc2}{Accuracy of best tree on small class }  
#'   \item{bamax}{Balanced accuracy of best tree} 
#'   \item{tamax}{Test accuracy of best tree}
#'   \item{kumin}{Number of elements with misclassified parts longer than 'ktest' for best tree}
#'   \item{elems}{Elements with long misclassified parts for best tree}
#'   \item{mindlong}{Indices of long misclassified parts for best tree}
#'   \item{ind1max}{Elements of 1st category of substructure used by best tree} 
#'   \item{ind2max}{Elements of 2nd category of substructure used by best tree} 
#'   \item{indmax}{Predictors used by best tree} 
#'   \item{bestTrain}{Training set for best tree} 
#'   \item{bestTest}{Test set for best tree} 
#'   \item{labs}{labs from Struc} 
#'   \item{lablarge}{Label of large class} 
#'   \item{labsmall}{Label of small class} 
#'   \item{vers}{Version used for structured subsampling} 
#'   \item{acc1E}{Accuracy of large class on Elements of best tree} 
#'   \item{acc2E}{Accuracy of small class on Elements of best tree}  
#'   \item{bamaxE}{Balanced accuracy of best tree on Elements} 
#'   \item{nam1}{Names of misclassified Elements of large class}
#'   \item{nam2}{Names of misclassified Elements of small class}
#'   \item{thr}{Threshold for element classification}
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.1, for further information.
#'
#' Standard output can be produced by means of \code{print(name$besttree)} or just \code{name$besttree} as well as \code{plot(name$besttree)} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export PrInDTCstruc
#'
#' @examples
#' data <- PrInDT::data_zero
#' data <- na.omit(data) # cleaned full data: no NAs
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}','MLU == {1, 3}') # split exclusions
#' # substructure
#' name <- PrInDT::participant_zero
#' check <- "data$ETH"
#' labs <- matrix(1:6,nrow=2,ncol=3)
#' labs[1,] <- c("C1a","C1b","C1c")
#' labs[2,] <- c("C2a","C2b","C2c")
#' rownames(labs) <- c("Children 1","Children 2")
#' Struc <- list(name=name,check=check,labs=labs)
#' out <- PrInDTCstruc(data,"real",ctestv,Struc,vers="c",weight=TRUE,N=5,Pit=5,conf.level=0.99)
#' out
#' plot(out)
#' # indrep = 1
#' Struc <- list(name=name)
#' out <- PrInDTCstruc(data,"real",Struc=Struc,vers="c",Pit=5,Eit=5,N=5,crit="ba",indrep=1,
#'            conf.level=0.99)
#'
PrInDTCstruc <- function(datain,classname,ctestv=NA,Struc=NA,vers="d",weight=FALSE,Eit=NA,Pit=NA,N=99,Ni=99,undersamp=TRUE,
  crit="ba",ktest=0,stest=integer(length=0),conf.level=0.95,indrep=0,minsplit=NA,minbucket=NA,repvar=NA,valdat=datain,thr=0.5){
#
## input check
if ( typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character"))
    || !(typeof(Struc) %in% c("logical","list")) || typeof(vers) != "character" || typeof(weight) != "logical" 
    || !(typeof(Eit) %in% c("logical","integer","double")) || !(typeof(Pit) %in% c("logical","integer","double")) || !(typeof(N) %in% c("logical","integer","double")) 
    || !(typeof(Ni) %in% c("logical","integer","double")) ||  typeof(undersamp) != "logical" ||  typeof(crit) != "character" ||  !(typeof(ktest) %in% c("integer","double"))
    || !(typeof(stest) %in% c("integer","double")) || !(0 < conf.level & conf.level <= 1) || !(typeof(indrep) %in% c("integer","double"))
    || !(typeof(minsplit) %in% c("logical","double")) || !(typeof(minbucket) %in% c("logical", "double")) || typeof(valdat) != "list" || typeof(thr) != "double"){
  stop("irregular input")
}
  if (!(all(names(datain) %in% names(valdat)))){
    stop("validation data variables unequal input data variables")
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
if (indrep > 0 & crit != "ba"){
  stop("irregular input: repeated measurements only implemented for ('crit = ba')")
}
if (indrep == 0 & vers != "b"){
  check <- eval(parse(text = Struc$check))
  if (length(check) == 0){
    cat("\n irregular input: check variable not in data set\n")
    return()
  }
  if (is.factor(check) != TRUE){
    cat("\n Error: 'Struc$check' has to be a factor variable\n")
    return()
  }
}
#
## re-naming
dataname <- as.character(substitute(datain))
data <- datain  
name <- Struc$name
if (is.factor(name) != TRUE){
  cat("\n Error: 'Struc$name' has to be a factor variable\n")
  return()
}
names(data)[names(data)==classname] <- "class"
names(valdat)[names(valdat)==classname] <- "class"
if (indrep > 0 & all(is.na(repvar)) == TRUE & identical(data,valdat)) {
  repv <- name
} else {
  repv <- repvar
}
if (!(identical(levels(data$class),levels(valdat$class)))){
  stop("levels of input class variable unequal levels of validation class variable")
}
# indrep = 2 indicates direct optimization of ba for substructure elements
if (indrep > 0){
  Struc$check <- "data$class"
  check <- eval(parse(text = Struc$check))
  labs <- matrix(1:2,nrow=2,ncol=1)
  labs[,1] <- levels(data$class) # as.character(unique(data$class))
  rownames(labs) <- as.character(labs[,1])
  undersamp <- FALSE  
  weight <- FALSE
}
if ( crit == "ba" & indrep > 0 ) {
  if (all(is.na(repv)) == FALSE & length(repv) != dim(valdat)[1] ){
    stop("irregular input: no. of observations different in substructure and validation set")
  }
}
if ( length(name) != dim(data)[1] | length(name) != length(check)){
  stop("irregular input: no. of observations different in substructure and data set")
}
## indrep dependencies
## Repeated measurement analysis applicable?
if (indrep > 0){
  conti <- table(Struc$name,data$class)
  if ( (sum(conti[,1] > 0) + sum(conti[,2] > 0)) > dim(conti)[1]){
    message("\n","Substructure variable has more than one class per value",
             "\n","Repeated measurement analysis not applicable","\n")
    return()
  }
}
D <- dim(data)[2] - 1
if (all(is.na(Pit)) == TRUE){
  Pit <- c(max(2,(D-3)):D)  
}
if (max(Pit) > D){stop("Too many predictors requested: Maximum is ",D)}
#
## 
## vers dependencies
if (vers == "a"){
  Ni <- 0
  Pit <- D
}
if (vers == "b"){
  N <- 0
#  Ni <- 99
  Eit <- 0
  pS <- 0
  pE <- 0
  labs <- NA
  undersamp <- FALSE
}
if (crit == "tal"){
  undersamp <- FALSE
  weight <- FALSE
}
if (vers != "b"){
  if (Ni > 1) {Ni <- N} ## ?????????
  labs <- Struc$labs
  if (anyNA(Struc[[1]]) != TRUE & indrep == 0){
    if (Struc$check == paste0(dataname,"$",names(data)[names(data)==classname])) { 
      check <- eval(parse(text = Struc$check)) 
      Struc$check <- "data$class"
      indrep <- 1
      undersamp <- FALSE
      weight <- FALSE
    }
  }
## 
  n_class1 <- table(data$class)[1] # no. of elements of class 1
  n_class2 <- table(data$class)[2] # no. of elements of class 2
  if (min(n_class1,n_class2) == 0){
    stop("irregular input: only 1 class")
  }
#print(table(data$class))
  if (n_class1 < n_class2){
    # relevel of classes if smaller class first
    valdat$class <- stats::relevel(valdat$class, levels(data$class)[2]) # larger class now first
    data$class <- stats::relevel(data$class, levels(data$class)[2]) # larger class now first
  #  lablarge <- names(table(data$class))[1] # class 1 = large
  #  labsmall <- names(table(data$class))[2] # class 2 = small
    n_class1 <- table(data$class)[1] # no. of elements of larger class 1
    n_class2 <- table(data$class)[2] # no. of elements of smaller class 2
  }
  if (indrep > 0){
    pS <- table(Struc$name,data$class)[,1] #
    pS <- pS[pS != 0]
    pS <- pS / pS
#  pS[pS != 0] <- 1
#print(length(pS))
    pE <- table(Struc$name,data$class)[,2] 
    pE <- pE[pE != 0]
    pE <- pE / pE
#print(length(pE))
  } else {
#    check <- eval(parse(text = Struc$check))  ## check not equal classname!!!
    p <- table(Struc$name,data$class)[,2] # /sum(table(Struc$name,data$class)[,2]) * 100
    pS <- p[unique(Struc$name[check %in% Struc$labs[1,]])] 
    pE <- p[unique(Struc$name[check %in% Struc$labs[2,]])]
    pS[pS == 0] <- 5e-10
    pE[pE == 0] <- 5e-10
  }
  ## weights 
  if (weight == FALSE){
    pS[pS > 5e-10] <- 1 
    pE[pE > 5e-10] <- 1
  }
  if (all(is.na(Eit)) == TRUE){
    Cl <- min(length(pS),length(pE))
    Eit <- c((Cl-4):Cl)
  }
}
if (max(Eit) > min(length(pS),length(pE))){stop("Too many substructure elements requested: Maximum is ",min(length(pS),length(pE)))}
#
## PrInDT function calls 
if ((vers == "d") | (vers == "a")){
  out <- PrInDTstruc1(data,classname,ctestv=ctestv,name=name,check=check,labs=labs,N=N,Ni=Ni,Eit=Eit,Pit=Pit,crit=crit,ktest=ktest,stest=stest,p1=pS,p2=pE,undersamp=undersamp,co=conf.level,minsplit=minsplit,minbucket=minbucket,repvar=repv,valdat=valdat,indrep=indrep,thr=thr)
  nam1 <- out$nam1
  nam2 <- out$nam2
}
if ((vers == "c") | (vers == "b")){
  out <- PrInDTstruc2(data,classname,ctestv=ctestv,name=name,check=check,labs=labs,N=N,Ni=Ni,Pit=Pit,Eit=Eit,crit=crit,ktest=ktest,stest=stest,p1=pS,p2=pE,undersamp=undersamp,co=conf.level,minsplit=minsplit,minbucket=minbucket,repvar=repv,valdat=valdat,indrep=indrep,thr=thr) 
  nam1 <- out$nam1
  nam2 <- out$nam2
}
#
## indrep = 1: 'ba' calculation for elements for models optimized for individual observations 
acc1E <- 0
acc2E <- 0
bamaxE <- 0
if (crit == "ba" & indrep == 1){
  nam1 <-  list()
  nam2 <- list()
  pred <- predict(out$modbest,newdata=valdat)  ## data???
  ch <- table(Struc[[1]],pred)
  conti <- cbind(table(Struc[[1]],valdat$class),ch)  ## data???
  no1 <- sum(conti[,1] > 0)
  no2 <- sum(conti[,2] > 0)
  for (i in 1:dim(conti)[1]){
    if (conti[i,2] > 0 & conti[i,4] <= (conti[i,2]*thr)){ 
      nam2 <- c(nam2,rownames(ch)[i])
    }
    if (conti[i,1] > 0 & conti[i,3] < (conti[i,1]*thr)){ 
      nam1 <- c(nam1,rownames(ch)[i])
    }
  }
  acc2E <- 1- length(nam2) / no2
  acc1E <- 1 - length(nam1) / no1  
  bamaxE <- (acc1E + acc2E)/2
  if (length(nam1) == 0) {nam1 <- "-"}
  if (length(nam2) == 0) {nam2 <- "-"}
}
#
## results
lablarge <- names(table(data$class))[1] # class 1 = large
labsmall <- names(table(data$class))[2] # class 2 = small 
result <- list(modbest = out$modbest, interp=out$interp, dmax=out$dmax, ntmax = out$ntmax, acc1 = out$acc1, acc2 = out$acc2, 
     bamax=out$bamax,tamax=out$tamax,kumin=out$kumin,elems=out$elems,mindlong=out$mindlong,
     ind1max=out$ind1max, ind2max=out$ind2max, indmax=out$indmax, bestTrain=out$gmaxTrain, bestTest=out$gmaxTest,labs=labs,
     lablarge=lablarge,labsmall=labsmall,crit=crit,vers=vers,indrep=indrep,acc1E = acc1E, acc2E =  acc2E, bamaxE = bamaxE, 
     nam1 = nam1, nam2 = nam2, thr=thr)
class(result) <- "PrInDTCstruc"
result
}
#' @export
print.PrInDTCstruc <- function(x,...){
  cat("Number of interpretable trees: ",x$interp[1]," of ",x$interp[2]," trees","\n","\n")
  cat("Best interpretable tree")
  print(x$modbest)
  cat("\n","Number of predictors in training set for best tree: ",x$dmax,"\n")
  cat("\n","Size of training set for best tree: ",x$ntmax,"\n")
  cat("\n","Predictors in model construction: ",x$indmax,"\n")
  if (x$vers != "b"){
  cat("\n","Elements used in model construction")
#  cat("\n","Number of observations of",rownames(x$labs)[1],"in model construction","\n")  ## labs falsch???
  if (x$indrep == 0) {
    cat("\n","Number of observations of",rownames(x$labs)[1],"in model construction: ", 
          sum(x$bestTrain$check[x$bestTrain$SubStruc %in% x$ind1max] %in% x$labs[1,]),"\n")
  } else {
    cat("\n","Number of observations of",rownames(x$labs)[1],"in model construction : ", 
          sum(x$bestTrain$check[x$bestTrain$SubStruc %in% x$ind1max] > 0),"\n")
  }
  Elements <- cbind(as.character(x$bestTrain$check[x$bestTrain$SubStruc %in% x$ind1max]),
    as.character(x$bestTrain$SubStruc[x$bestTrain$SubStruc %in% x$ind1max]))
  print(table(Elements,exclude=x$labs))
  if (x$indrep == 1){
    print(table(Elements,exclude=x$ind1max))
  } else {
    print(table(Elements,exclude=c(x$ind1max,x$labs[2,])))
  }
#  cat("Number of observations of",rownames(x$labs)[2],"in model construction","\n")
  if (x$indrep == 0) {
    cat(" Number of observations of",rownames(x$labs)[2],"in model construction: ",
           sum(x$bestTrain$check[x$bestTrain$SubStruc %in% x$ind2max] %in% x$labs[2,]),"\n")
  } else {
    cat(" Number of observations of",rownames(x$labs)[2],"in model construction: ",
           sum(x$bestTrain$check[x$bestTrain$SubStruc %in% x$ind2max] > 0),"\n")
  }
  Elements <- cbind(as.character(x$bestTrain$check[x$bestTrain$SubStruc %in% x$ind2max]),
    as.character(x$bestTrain$SubStruc[x$bestTrain$SubStruc %in% x$ind2max]))
  print(table(Elements,exclude=x$labs))
  if (x$indrep > 0){
    print(table(Elements,exclude=x$ind2max))
  } else {
    print(table(Elements,exclude=c(x$ind2max,x$labs[1,])))
  }
#  cat("\n","Test accuracy of best tree","\n")
#  class.pred <- predict(x$modbest,newdata=x$bestTest)
#  cat(" ",sum(class.pred == x$bestTest$class) / length(class.pred))
}
#  cat("Accuracy of best tree")
#  class.pred <- predict(x$outmax$treeAll,newdata=data)
# print(sum(class.pred == data$class) / length(class.pred))
if (x$crit == "ba"){
  accs <- c(round(x$acc1,digits=6),round(x$acc2,digits=6),round(x$bamax,digits=6))
  names(accs) <- c(x$lablarge,x$labsmall,"balanced")
  if (x$indrep < 2){
    cat("\n","Validation sample accuracies of best tree (observations) ","\n")
    print(accs,row.names=FALSE)
  } else {
    cat("\n","Validation sample accuracies of best tree (Elements) ","\n")
    print(accs,row.names=FALSE)
  }    
#  cat("\n","Validation sample accuracies of best tree (observations)","\n",x$lablarge,"\t",x$labsmall,"\t","balanced","\n")
} else {
  accs <- c(round(x$acc1,digits=6),round(x$acc2,digits=6),round(x$bamax,digits=6))
  names(accs) <- c(substr(x$lablarge,1,12),substr(x$labsmall,1,12),"balanced")
  cat("\n","Test sample accuracies of best tree (observations) ","\n")
  print(accs,row.names=FALSE)
#  cat("\n","Test accuracies of best tree","\n"," ",substr(x$lablarge,1,12),"\t",substr(x$labsmall,1,12),"\t"," balanced","\n")
}
#  cat(" ",c(round(x$acc1,4),"\t",round(x$acc2,4),"\t",round(x$bamax,4)),"\n")
if (x$crit == "ta" | x$crit == "tal"){
  cat("\n","Overall test accuracy of best tree: ",x$tamax,"\n")
}
if (x$crit == "tal"){
  cat("\n","Optimal number of long misclassified parts: ",x$kumin)
  cat("\n","Elements with long misclassified parts:")
  cat("\n",x$elems,"\n")
}
if (x$indrep == 1 & x$crit == "ba"){
  cat("\n","Repeated measurements: Element classification: Threshold = ",x$thr)
  accs <- c(round(x$acc1E,digits=6),round(x$acc2E,digits=6),round(x$bamaxE,digits=6))
  names(accs) <- c(x$lablarge,x$labsmall,"balanced")
  cat("\n","Input sample accuracies of best tree (Elements)","\n")
  print(accs,row.names=FALSE)
#  cat("\n","Input sample accuracies of best tree (Elements)","\n",x$lablarge,"\t",x$labsmall,"\t","balanced","\n") 
#  cat(" ",c(round(x$acc1E,4),"\t","  ",round(x$acc2E,4),"\t","  ",round(x$bamaxE,4)),"\n\n")
  cat("Wrongly predicted Elements:","\n")
  cat(x$lablarge,"\n")
  cat(unlist(x$nam1),"\n")
  cat(x$labsmall,"\n")
  cat(unlist(x$nam2),"\n")
}
if (x$indrep == 2 & x$crit == "ba"){
    cat("Wrongly predicted Elements:","\n")
  cat(x$lablarge,"\n")
  cat(unlist(x$nam1),"\n")
  cat(x$labsmall,"\n")
  cat(unlist(x$nam2),"\n")
}
}
#' @export
plot.PrInDTCstruc <- function(x,...){
  plot(x$modbest,main="Best interpretable tree")
}