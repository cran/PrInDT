#' Two-stage estimation for regression
#'
#' @description The function R2SPrInDT applies 'N' repetitions of subsampling for finding an optimal subsample to model
#' the relationship between the continuous variables with indices 'inddep' and all other factor and numerical variables
#' in the data frame 'datain'. \cr
#' Subsampling of observations and predictors uses the percentages in 'pobs1' and 'ppre1', respectively, at stage 1, and the percentages 'pobs2' and 'ppre2' 
#' at stage 2, accordingly.
#' The optimization criterion is the goodness of fit R2 on the full sample.  \cr
#' The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage R2SPrInDT(data,ctestv=NA,inddep,N=99,pobs1=c(0.90,0.70),ppre1=c(0.90,0.70),
#'                 pobs2=pobs1,ppre2=ppre1,conf.level=0.95,minsplit=NA,minbucket=NA)
#'
#' @param data Input data frame with continuous target variable 'regname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param inddep Column indices of target variables in datain
#' @param N Number of repetitions of subsampling from predictors (integer) in versions "b" and "c";\cr
#'      default = 99
#' @param pobs1 Percentage(s) of observations for  subsampling at stage 1;\cr
#'      default=c(0.9,0.7)
#' @param ppre1 Percentage(s) of predictors for subsampling at stage 1;\cr
#'      default=c(0.9,0.7)
#' @param pobs2 Percentage(s) of observations for  subsampling at stage 2";\cr
#'      default=pobs1
#' @param ppre2 Percentage(s) of predictors for subsampling at stage 2;\cr
#'      default=ppre1
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
#'   \item{depnames}{names of dependent variables}
#'   \item{R2both}{R2s of best trees at both stages}  
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.6.1, for further information.
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{name} as well as \code{plot(name} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export R2SPrInDT
#'
#' @examples
#' data <- PrInDT::data_vowel
#' data <- na.omit(data)
#' CHILDvowel <- data$Nickname
#' data$Nickname <- NULL
#' syllable <- 3 - data$syllables
#' data$syllabels <- NULL
#' data$syllables <- syllable
#' data$speed <- data$word_duration / data$syllables
#' names(data)[names(data) == "target"] <- "vowel_length"
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}','MLU == {1, 3}') # split exclusions
#' inddep <- c(13,9) 
#' out2SR <- R2SPrInDT(data,ctestv=ctestv,inddep=inddep,N=9,conf.level=0.99)
#' out2SR
#' plot(out2SR)
#'
R2SPrInDT <- function(data,ctestv=NA,inddep,N=99,pobs1=c(0.90,0.70),ppre1=c(0.90,0.70),pobs2=pobs1,ppre2=ppre1,conf.level=0.95,minsplit=NA,minbucket=NA){
# input check
if ( typeof(data) != "list" || !(typeof(ctestv) %in% c("logical", "character")) || !(typeof(inddep) %in% c("integer","double"))
    || !all(0 < pobs1 & pobs1 <= 1) || !all(0 < ppre1 & ppre1 <= 1)
    || !all(0 < pobs2 & pobs2 <= 1) || !all(0 < ppre2 & ppre2 <= 1) 
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
depnames <- colnames(data)[inddep]
anz <- length(inddep)
D <- dim(data)[2]
## only exogenous predictors
message("1st stage")
set.seed(7654321)
models <- list()
ctpreds <- matrix(0,ncol=anz,nrow=dim(data)[1])
R2 <- c(1:anz)*0
for (i in 1:anz){
  indpre <- c(1:(dim(data)[2]))[-inddep]
  lind <- length(indpre)
  indmid <- c(min(inddep):max(inddep))[-(inddep-min(inddep)+1)]
  indpre <- c(indpre[1:(min(inddep)-1)],indmid,inddep[i],indpre[(min(inddep)+length(indmid)):lind])
#  indpre <- c(indpre[1:9],inddep[i],indpre[10:17])
# print(indpre)
  message(depnames[i],": 1st stage")
#outreg <- PrInDTreg(datavowel,"vowel_length",ctestv,N,pobs,ppre,conf.level=cf) ## all vowel related variables
  outreg <- PrInDTreg(data[,indpre],depnames[i],ctestv,N=N,pobs=pobs1,ppre=ppre1,conf.level=conf.level,seedl=FALSE,minsplit=minsplit,minbucket=minbucket)
#  print(outreg) 
#  plot(outreg)
  models <- c(models,outreg$ctmax)
  ctpredsVL <- predict(outreg$ctmax,newdata=data)
  ctpreds[,i] <- ctpredsVL
  R2[i] <- outreg$maxR2
# print(R2[i])
}
#
## 2nd stage estimations
#
message("\n","2nd stage")
set.seed(7654321)
dataP <- data
dataP <- cbind(dataP,ctpreds)
R2P <- R2
modelsP <- models
for (i in 1:anz){
  names(dataP)[(D+i)] <- paste0(depnames[i],"_2S")
}
models <- list()
ctpreds <- matrix(0,ncol=anz,nrow=dim(dataP)[1])
R2 <- c(1:anz)*0
for (i in 1:anz){
  indpre <- c(1:(dim(dataP)[2]))[-(D+i)]
  indpre <- indpre[-inddep]
  indpre <- c(inddep[i],indpre)
#print(indpre)
  message(depnames[i],": 2nd stage")
  outreg <- PrInDTreg(dataP[,indpre],depnames[i],ctestv,N=N,pobs=pobs2,ppre=ppre2,conf.level=conf.level,seedl=FALSE,minsplit=minsplit,minbucket=minbucket)
  models <- c(models,outreg$ctmax)
  ctpredsVL <- predict(outreg$ctmax,newdata=dataP)
  ctpreds[,i] <- ctpredsVL
  R2[i] <- outreg$maxR2
}
R2P <- rbind(R2P,R2)
colnames(R2P) <- depnames
result <- list(models1=modelsP,models2=models,depnames=depnames,R2both=R2P) 
class(result) <- "R2SPrInDT"
result
}
#' @export
print.R2SPrInDT <- function(x,...){
## output
cat("Accuracies of best models","\n")
colnames(x$R2both) <- c(x$depnames)
rownames(x$R2both) <- c("1st stage","2nd stage")
print(x$R2both)
cat("\n\n")
for (i in 1:length(x$depnames)){
  cat("Best model 2nd stage: ",x$depnames[i], "\n")
  print(x$models2[[i]])
  cat("\n")
}
}
#' @export
plot.R2SPrInDT <- function(x,...){
for (i in 1:length(x$depnames)){
  plot(x$models2[[i]],main=x$depnames[i])
}
}