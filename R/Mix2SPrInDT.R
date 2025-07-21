#' Two-stage estimation for classification-regression mixtures
#'
#' @description The function Mix2SPrInDT applies 'N' repetitions of subsampling for finding an optimal subsample to model
#' the relationship between the dependent variables specified in the sublist 'targets' of the list 'datalist' and all other factor and numerical variables
#' in the corresponding data frame specified in the sublist 'datanames' of the list 'datalist' in the same order as 'targets'. \cr
#' The function is prepared to handle classification tasks with 2 or more classes and regression tasks. At first stage, the targets are estimated only on
#' the basis of the exogenous predictors. At second stage, summaries of the predictions of the endogenous features from the first stage models are added as new predictors.\cr
#' Subsampling of observations and predictors uses the percentages specified in the matrix 'percent', one row per estimation taks. For classification tasks, 
#' percentages 'percl' and 'percs' for the larger and the smaller class have to be specified. For regression tasks, 'pobs' and 'ppre' have to specified 
#' for observations and predictors, respectively.\cr
#' For generating summaries, the variables representing the substructure have to be specified in the sublist 'datastruc' of 'datalist'.\cr
#' The sublist 'summ' of 'datalist' includes for each discrete target the classes for which you want to calculate the summary percentages and for each 
#' continuous target just NA (in the same order as in the sublist 'targets'). For a discrete target in this list, you can provide a sublist of classes to be combined 
#' (for an example see the below example).\cr
#' The optimization citerion is the balanced accuracy for classification and R2 for regression, both evaluated on the full sample.  \cr
#' The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage Mix2SPrInDT(datalist,ctestv=NA,N=99,percent=NA,conf.level=0.95,minsplit=NA,minbucket=NA)
#'
#' @param datalist list(datanames,targets,datastruc,summ) Input data: For specification see the above description
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param N Number of repetitions of subsampling for predictors (integer); default = 99
#' @param percent matrix of percent spefications: For specification see the above description; default: 'percent = NA' meaning default values for percentages.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1); default = 0.95
#' @param minsplit Minimum number of elements in a node to be splitted; default = 20
#' @param minbucket Minimum number of elements in a node; default = 7
#'
#' @return
#' \describe{
#'   \item{models1}{Best trees at stage 1} 
#'   \item{models2}{Best trees at stage 2} 
#'   \item{depnames}{names of dependent variables}
#'   \item{nmod}{number of models of tasks}
#'   \item{nlev}{levels of tasks}
#'   \item{accAll}{accuracies of best trees at both stages}  
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.6.1, for further information.
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{name} as well as \code{plot(name} where 'name' is the output data 
#' frame of the function.\cr
#'
#' @export Mix2SPrInDT
#' @importFrom stats aggregate
#' @importFrom stringr str_split
#' @importFrom gdata mv
#'
#' @examples
#' # zero data
#' datazero <- PrInDT::data_zero
#' datazero <- na.omit(datazero) # cleaned full data: no NAs
#' names(datazero)[names(datazero)=="real"] <- "zero"
#' CHILDzero <- PrInDT::participant_zero
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}','MLU == {1, 3}') # split exclusions
#' ##
#' # multi-level data
#' datastrat <- PrInDT::data_zero
#' datamult <- na.omit(datastrat)
#' # ctestv <- NA
#' datamult$mult[datamult$ETH %in% c("C1a","C1b","C1c") & datamult$real == "zero"] <- "zero1"
#' datamult$mult[datamult$ETH %in% c("C2a","C2b","C2c") & datamult$real == "zero"] <- "zero2"
#' datamult$mult[datamult$real == "realized"] <- "real"
#' datamult$mult <- as.factor(datamult$mult) # mult is new class variable
#' datamult$real <- NULL # remove old class variable
#' CHILDmult <- CHILDzero
#' ##
#' # vowel data
#' data <- PrInDT::data_vowel
#' data <- na.omit(data)
#' CHILDvowel <- data$Nickname
#' data$Nickname <- NULL
#' syllable <- 3 - data$syllables
#' data$syllabels <- NULL
#' data$syllables <- syllable
#' data$speed <- data$word_duration / data$syllables
#' names(data)[names(data) == "target"] <- "vowel"
#' datavowel <- data
#' ##
#' # function preparation and call
#' datanames <- list("datazero","datamult","datavowel")
#' targets <- c("zero","mult","vowel")
#' datastruc <- list(CHILDzero,CHILDmult,CHILDvowel)
#' summult <- paste("zero1","zero2",sep=",")  
#' summ <- c("zero",summult,NA)
#' datalist <- list(datanames=datanames,targets=targets,datastruc=datastruc,summ=summ)
#' percent <- matrix(NA,nrow=3,ncol=2)
#' percent[1,] <- c("percl=0.075","percs=0.9") # percentages for datazero
#' # no percentages needed for datapast
#' percent[3,] <- c("pobs=0.9","ppre=c(0.9,0.8)") # percentages for datavowel
#' out2SMix <- Mix2SPrInDT(datalist,ctestv=ctestv,N=19,percent=percent,conf.level=0.99)
#' out2SMix
#' plot(out2SMix)
#'
Mix2SPrInDT <- function(datalist,ctestv=NA,N=99,percent=NA,conf.level=0.95,minsplit=NA,minbucket=NA){
## input check
if ( typeof(datalist) != "list" || !(typeof(ctestv) %in% c("logical", "character")) 
    || !(typeof(percent) %in% c("logical","character"))
    || !(typeof(N) %in% c("integer","double")) || !(0 < conf.level & conf.level <= 1) || !(typeof(minsplit) %in% c("logical","double")) 
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
######################
#### interdependent system of classification and regression models
## estimate with observables only
#
  anz <- length(datalist$datanames)
  if (all(is.na(percent) == TRUE)){
    percent <- matrix(NA,nrow=anz,ncol=2)
  }
  nmod <- c(1:anz)
  nlev <- list()
  acc <- c(1:anz)*0
  datat <- datalist[[1]]
  summ <- datalist[[4]]
  inds <- list()
  extdata <- list()
  modelsP <- list()
  for (K in 1:anz){
    data <- as.data.frame(eval(parse(text = datat[K])))    
    dataname <- datat[K]
    target <- datalist$targets[K]
    ind <- c(1:dim(data)[1])
#print(levels(as.factor(data[,target])))
message(dataname,": 1st stage")
# 2 classes: call of PrInDT
    if (length(levels(data[,target])) == 2){
      nmod[K] <- 1
      nlev <- c(nlev,"class")
# assigning percl, percs
      percl <- NA
      percs <- 1
      x <- str_split(percent[K,1],"=")
      assign(x[[1]][1],eval(parse(text=percent[K,1])))
      x <- str_split(percent[K,2],"=")
      assign(x[[1]][1],eval(parse(text=percent[K,2])))
      if (!(all(0 < percl & percl <= 1) || is.na(percl)==TRUE) || !(all(0 < percs & percs <= 1))){
        stop("irregular input for percl or percs")
      } 
      out2S <- PrInDT(data,target,ctestv,N,percl=percl,percs=percs,conf.level,minsplit=minsplit,minbucket=minbucket)
      sub <- unlist(datalist$datastruc[K])
      sub <- droplevels(sub)
      s <- unlist(str_split(summ[K],","))
#      si <- which(levels(data[,target])==s)
#      snotreal <- table(sub,data[,target])[,si]
#      sc <- as.character(sub)
#      z1 <- aggregate(data[,target] ~ sc,data=cbind(sc,data[,target]),FUN="length")
      ctpreds <- predict(out2S$tree1st,newdata=data)  ## !!
      si <- which(levels(ctpreds)==s)
      snotreal <- table(sub,ctpreds)[,si]
      sc <- as.character(sub)
      z1 <- aggregate(ctpreds ~ sc,data=cbind(sc,ctpreds),FUN="length")
      ind <- snotreal / z1[,2]
      inds[[K]] <- ind
#      print(length(inds[[K]]))
      extdata[[K]] <- paste0(dataname,"i")
      acc[K] <- out2S$ba1st[3]
      modelsP <- c(modelsP,out2S$tree1st)
    } 
# more than 2 classes: call of PrInDTMulev
    if (length(levels(data[,target])) > 2){
      nmod[K] <- length(levels(data[,target]))
      nlev <- c(nlev,levels(data[,target]))
      out2S <- PrInDTMulev(data,target,ctestv,N,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket) # percentages chosen automatically
      sub <- unlist(datalist$datastruc[K])
      sub <- droplevels(sub)
      s <- unlist(str_split(summ[K],","))
#      si <- which(levels(data[,target]) %in% s)
#      snotmarked <- table(sub,data[,target])[,si[1]] + table(sub,data[,target])[,si[2]]
#      sc <- as.character(sub)
#      p1 <- aggregate(data[,target] ~ sc,data=cbind(sc,data[,target]),FUN="length")
      ctpreds <- predict(out2S,target,newdata=data)
#print(table(ctpreds,data[,target]))
      si <- which(levels(ctpreds) %in% s)
      snotmarked <- table(sub,ctpreds)[,si[1]]
      if (length(s) > 1){
        for (i in 2:length(s)){
          snotmarked <- snotmarked + table(sub,ctpreds)[,si[i]]
        }
      }
      sc <- as.character(sub)
      ctpreds <- as.character(ctpreds)
      p1 <- aggregate(ctpreds ~ sc,data=cbind(sc,ctpreds),FUN="length")
      ind <- snotmarked / p1[,2]
      inds[[K]] <- ind
#     print(length(inds[[K]]))
      extdata[[K]] <- paste0(dataname,"i")
      acc[K] <- out2S$ba
      modelsP <- c(modelsP,out2S$trees)
     }
# regression: call of PrInDTreg
     if (length(levels(data[,target])) == 0){
       nmod[K] <- 1
       nlev <- c(nlev,"real")
# assigning pobs, ppre
       pobs <- c(0.9,0.7)
       ppre <- c(0.9,0.7) 
       x <- str_split(percent[K,1],"=")
       assign(x[[1]][1],eval(parse(text=percent[K,1])))
       x <- str_split(percent[K,2],"=")
      assign(x[[1]][1],eval(parse(text=percent[K,2])))
      if (!(all(0 < pobs & pobs <= 1)) || !(all(0 < ppre & ppre <= 1))){
          stop("irregular input for pobs or ppre")
      } 
      out2S <- PrInDTreg(data,target,ctestv,N,pobs=pobs,ppre=ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
      sub <- unlist(datalist$datastruc[K])
      sc <- as.character(sub)
#      v1 <- aggregate(data[,target] ~ sc,data=cbind(sc,data[,target]),FUN="mean")
      ctpreds <- predict(out2S$ctmax,newdata=data)
      v1 <- aggregate(ctpreds ~ sc,data=cbind(sc,ctpreds),FUN="mean")
      ind <- v1[,2]
      names(ind) <- v1[,1]
      inds[[K]] <- ind
#      print(length(inds[[K]]))
      extdata[[K]] <- paste0(dataname,"i")
      acc[K] <- out2S$maxR2
      modelsP <- c(modelsP,out2S$ctmax)
    }
#   cat("Best model for",target,"\n")
#   print(out2S)
#   plot(out2S)
 }
accP <- acc
#print(accP)
#print(extdata)
#################
## 2nd stage
#################
## extended datasets
for (K in 1:anz){
  datai <- as.data.frame(eval(parse(text = datat[K])))    
  mv("datai",extdata[[K]])
}
## name definitions
for (K in 1:anz){
  datai <- eval(parse(text = extdata[K]))
  for (i in 1:anz){
    if (i != K){
      ind <- inds[[i]]
      D <- dim(datai)[2]
      datai$indi <- mean(ind) ## NA
      for (j in 1:length(ind)) {
        datai$indi[unlist(datalist$datastruc[K]) == names(ind)[j]] <- ind[j]
      }
   names(datai)[(D+1)] <- paste0("p",datalist$targets[i])
   }
 }
mv("datai",extdata[[K]])
#print(str(eval(parse(text = extdata[K]))))
} 
#
#### estimation
#
acc <- c(1:anz)*0
models <- list()
message("\n","Second stage")
for (K in 1:anz){
   datai <- eval(parse(text = extdata[K]))
   target <- datalist$targets[K]
message(datat[K],": 2nd stage")
## 2 classes
   if (length(levels(datai[,target])) == 2){
     x <- str_split(percent[K,1],"=")
     assign(x[[1]][1],eval(parse(text=percent[K,1])))
     x <- str_split(percent[K,2],"=")
     assign(x[[1]][1],eval(parse(text=percent[K,2])))
     out2S <- PrInDT(datai,target,ctestv,N,percl=percl,percs=percs,conf.level,minsplit=minsplit,minbucket=minbucket) # best percentages found
     acc[K] <- out2S$ba1st[3]
     models <- c(models,out2S$tree1st)
  }
## more than 2 classes
  if (length(levels(datai[,target])) > 2){
      datai[,target] - as.factor(datai[,target])
      out2S <- PrInDTMulev(datai,target,ctestv,N,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket) # percentages chosen automatically
      message("balanced accuracies of class models: ",paste(round(out2S$acc,3),collapse=", "))  ## individual balanced accuracies
      acc[K] <- out2S$ba
      models <- c(models,out2S$trees)
  }
## regression
  if (length(levels(datai[,target])) == 0){
    x <- str_split(percent[K,1],"=")
    assign(x[[1]][1],eval(parse(text=percent[K,1])))
    x <- str_split(percent[K,2],"=")
    assign(x[[1]][1],eval(parse(text=percent[K,2])))
    out2S <- PrInDTreg(datai,target,ctestv,N,pobs=pobs,ppre=ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
    acc[K] <- out2S$maxR2
    models <- c(models,out2S$ctmax)
  }
#   cat("Best model for",target,"\n")
#   print(out2S)
#   plot(out2S)
}
accP <- cbind(accP,acc)
#rownames(accP) <- datalist$targets
#colnames(accP) <- c("acc","acc2S")
result <- list(models1=modelsP,models2=models,depnames=unlist(datalist$target),nmod=nmod,nlev=nlev,accAll=accP)
class(result) <- "Mix2SPrInDT"
result
}
#' @export
print.Mix2SPrInDT <- function(x,...){
## output
cat("Accuracies of best models","\n")
rownames(x$accAll) <- c(x$depnames)
colnames(x$accAll) <- c("1st stage","2nd stage")
print(x$accAll)
cat("\n\n")
J <- 0
for (i in 1:length(x$depnames)){
  if (x$nmod[i] == 1){
    J <- J + 1
    cat("Best model 2nd stage: ",x$depnames[i], "\n") 
    print(x$models2[[i]])
    cat("\n\n")
  } else {
    K <- x$nmod[i]
    for (k in 1:K){
       J <- J + 1
       cat("Best model 2nd stage: ",x$depnames[i],":",x$nlev[[J]],"\n") 
       print(x$models2[[i]])
       cat("\n\n")
    }
  }
}}
#' @export
plot.Mix2SPrInDT <- function(x,...){
J <- 0
for (i in 1:length(x$depnames)){
  if (x$nmod[i] == 1){
    J <- J + 1
    plot(x$models2[[i]],main=paste0("2nd stage: ",x$depnames[i]))
  } else {
    K <- x$nmod[i]
    for (k in 1:K){
       J <- J + 1
       plot(x$models2[[i]],main=paste0("2nd stage: ",x$depnames[i],":",x$nlev[[J]]))
    }
  }
}}