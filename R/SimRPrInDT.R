#' Interdependent estimation for regression
#'
#' @description The function \code{\link{SimRPrInDT}} applies structured subsampling for finding an optimal subsample to model
#' the relationship between the continuous variables with indices 'inddep' and all other factor and numerical variables
#' in the data frame 'datain'. \cr
#' The substructure of the observations used for subsampling is specified by the list 'Struc' which consists of the variable 'name' representing the substructure,
#' the name 'check' of the variable with the information about the categories of the substructure, and the matrix 'labs' which specifies the values of 'check'
#' corresponding to two categories in its rows, i.e. in 'labs[1,]' and 'labs[2,]'. The names of the categories have to be specified by \code{rownames(labs)}.\cr
#' In structured subsampling first 'M' repetitions of subsampling of the variable 'name' with 'nsub' different elements of the substructure are realized. If 'nsub' is a list, each entry is employed individually. Then,
#' for each of the subsamples 'N' repetitions of subsampling of 'ppre' percentages of the predictors are carried out.\cr 
#' Subsampling of observations can additionally be restricted to 'pobs' percentages.\cr
#' The optimization citerion is the goodness of fit R2 on the full sample. At stage 2, the models are optimized individually. 
#' At stage 3, the mean of accuracies is optimized over all models.\cr
#' Struc=NA causes random subsampling of observations instead of structured subsampling.\cr
#' The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage SimRPrInDT(data,ctestv=NA,Struc=NA,inddep,N=99,pobs=0.9,ppre=c(0.9,0.7),
#'                 M=1,nsub=1,conf.level=0.95,minsplit=NA,minbucket=NA)
#'
#' @param data Input data frame with continuous target variables with column indices 'inddep' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param Struc = list(name,check,labs), cf. description for explanations
#' @param inddep Column indices of target variables in datain
#' @param N Number of repetitions of subsampling (integer) of predictors; default = 99
#' @param pobs Percentage(s) of observations for subsampling; default=0.9
#' @param ppre Percentage(s) of predictors for subsampling; default=c(0.9,0.7)
#' @param M Number of repetitions of subsampling of elements of substructure
#' @param nsub (List of) numbers of different elements of substructure per subsample
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#'
#' @return
#' \describe{
#'   \item{modelsF}{Best trees at stage 1} 
#'   \item{modelsI}{Best trees for the different values of 'nsub' at stage 2}
#'   \item{modelsJ}{Best trees for the different values of 'nsub' after mean optimization}
#'   \item{Struc}{Used structure}
#'   \item{msub}{Best numbers of elements in substructure: 2nd stage, 3rd stage}
#'   \item{depnames}{names of dependent variables}
#'   \item{R2All}{R2s of best trees at stages 1, 2, mean max.}  
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.6.2, for further information.
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{name} as well as \code{plot(name} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export SimRPrInDT
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
#' # structure definition
#' name <- CHILDvowel
#' check <- "data$ETH"
#' labs <- matrix(1:6,nrow=2,ncol=3)
#' labs[1,] <- c("C1a","C1b","C1c")
#' labs[2,] <- c("C2a","C2b","C2c")
#' rownames(labs) <- c("children 1","children 2")
#' Struc <- list(name=name,check=check,labs=labs)
#' # column indices of dependent variables
#' inddep <- c(13,9) 
#' outSimR <- SimRPrInDT(data,ctestv=ctestv,Struc=Struc,inddep=inddep,N=3,M=2,
#'                    nsub=c(19,20),conf.level=0.99)
#' outSimR
#' plot(outSimR)
#'
SimRPrInDT <- function(data,ctestv=NA,Struc=NA,inddep,N=99,pobs=0.9,ppre=c(0.9,0.7),M=1,nsub=1,conf.level=0.95,minsplit=NA,minbucket=NA){
## input check
if ( typeof(data) != "list" || !(typeof(ctestv) %in% c("logical", "character"))
    || !(typeof(Struc) %in% c("logical","list")) || !all(0 < pobs & pobs <= 1) || !all(0 < ppre & ppre <= 1) 
    || !(typeof(N) %in% c("integer","double")) || !(typeof(M) %in% c("integer","double")) 
    || !(typeof(inddep) %in% c("integer","double")) || !(typeof(nsub) %in% c("integer","double"))
    || !(0 < conf.level & conf.level <= 1) || !(typeof(minsplit) %in% c("logical","double")) || 
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
depnames <- colnames(data)[inddep]
anz <- length(inddep)
D <- dim(data)[2]
if (any(is.na(Struc) != TRUE)){
  name <- Struc$name
#check <- Struc$check
  check <- eval(parse(text = Struc$check))
  labs <- Struc$labs
} else {
  M <- 1
  nsub <- 1
}
#message ("Simultaneous modeling")
models <- list()
modmax <- as.list(1:(anz*(1+length(nsub)))) 
dim(modmax) <- c((1+length(nsub)),anz) 
ctpreds <- matrix(0,ncol=anz,nrow=dim(data)[1])
R2 <- c(1:anz)*0
max01 <- matrix(0,nrow=(length(nsub)+1),ncol=(anz+1)) # nrow = full, all sizes
##
message("1st stage: full sample")
for (i in 1:anz){
#  message(depnames[i])
  #outreg <- PrInDTreg(datavowel,"vowel_length",ctestv,N,pobs,ppre,conf.level=cf) ## all vowel related variables
  outreg <- PrInDTregAll(data,depnames[i],ctestv,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
  models[[i]] <- outreg$treeAll
  ctpredsVL <- predict(outreg$treeAll,newdata=data)
  ctpreds[,i] <- ctpredsVL
  R2[i] <- outreg$R2All
  modmax[[1,i]] <- models[[i]]
}
max01[1,] <- c(R2,mean(R2))
colnames(max01) <- c(depnames,"mean")
#print(max01[1,])
R2P <- max01[1,]
modelsP <- models
##
message("\n","2nd stage: individual optimization")
set.seed(7654321)
k <- 1
for (i in nsub){
  if (any(is.na(Struc) != TRUE)){
    message("Number of elements: ",i)
  }
  k <- k + 1
  max01[k,] <- max01[1,]
  for (n in 1:M){
    if (any(is.na(Struc) != TRUE)){
      message("repetition ",n)
      ind1 <- sample(unique(name[check %in% labs[1,]]))[1:i] 
      ind2 <- sample(unique(name[check %in% labs[2,]]))[1:i] 
      datas <- data[(as.factor(name) %in% ind1) | (as.factor(name) %in% ind2) ,]
##  datas <- cbind(dataMend[,sample(colnames(dataMend))[1:psize]],datavowel[,c(8,9,10,13)])
    } else {
      datas <- data
    }
    for (j in 1:anz){
      outland <- suppressMessages(PrInDTreg(datas,depnames[j],ctestv=ctestv,N,pobs,ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket))
      ctpredsVL <- predict(outland$ctmax,newdata=data)
      R2[j] <- 1 - sum((data[,inddep[j]] - ctpredsVL)^2) / sum((data[,inddep[j]]-mean(data[,inddep[j]]))^2) ## R2 !!!!!
      models[[j]] <- outland$ctmax
      if (R2[j] > max01[k,j]) { 
        max01[k,j] <- R2[j]
        max01[k,(anz+1)] <- mean(max01[k,1:anz])
#        print(max01[k,])
        modmax[[k,j]] <- models[[j]]
      }
   }  
 } 
}
colnames(max01) <- c(depnames,"mean")
#R2I <- max01
#modelsI <- modmax
k2 <- which.max(max01[2:k,(anz+1)])
R2I <- max01[(k2+1),]
modelsI <- modelsP
for (i in 1:anz){
  modelsI[[i]] <- modmax[[(k2+1),i]]
}
#print(k)
#for (i in 1:anz){
#  plot(unlist(modmax[[k,i]]))
#}
###  optimization of mean
message("\n","3rd stage: joint optimization")
k  <- 1
for (i in nsub){
  if (any(is.na(Struc) != TRUE)){
    message("Number of elements: ",i)
  }
  k <- k + 1
  max01[k,] <- max01[1,]
  R21m <- sum(max01[1,1:anz]) / anz
  for (n in 1:M) {
    set.seed(4321 * n + n)
    if (any(is.na(Struc) != TRUE)){
      message("repetition ",n)
      ind1 <- sample(unique(name[check %in% labs[1,]]))[1:i] 
      ind2 <- sample(unique(name[check %in% labs[2,]]))[1:i] 
      datas <- data[(as.factor(name) %in% ind1) | (as.factor(name) %in% ind2) ,]
##  datas <- cbind(dataMend[,sample(colnames(dataMend))[1:psize]],datavowel[,c(8,9,10,13)])
    } else {
      datas <- data
    }
    for (j in 1:anz){
      outland <- suppressMessages(PrInDTreg(datas,depnames[j],ctestv=ctestv,N,pobs,ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket))
      ctpredsVL <- predict(outland$ctmax,newdata=data)
      R2[j] <- 1 - sum((data[,inddep[j]] - ctpredsVL)^2) / sum((data[,inddep[j]]-mean(data[,inddep[j]]))^2) ## R2 !!!!!
      models[[j]] <- outland$ctmax
    } 
    R2m <- sum(R2[1:anz]) / anz
    if (R2m > R21m) { 
      max01[k,1:anz] <- R2
      max01[k,(anz+1)] <- mean(max01[k,1:anz])
#      print(max01[k,])
      for (j in 1:anz){
        modmax[[k,j]] <- models[[j]]
      }
      R21m <- R2m
    }
  }
}
colnames(max01) <- c(depnames,"mean")
k3 <- which.max(max01[2:k,(anz+1)])
R2J <- max01[(k3+1),]
modelsJ <- modelsP
for (i in 1:anz){
  modelsJ[[i]] <- modmax[[(k3+1),i]]
}
#for (i in 1:anz){
#  plot(unlist(modmax[[k,i]]))
#}
R2All <- rbind(R2P,R2I,R2J)
result <- list(modelsF=modelsP,modelsI=modelsI,modelsJ=modelsJ,Struc=Struc,msub=c(nsub[k2],nsub[k3]),depnames=depnames,R2All=R2All)
class(result) <- "SimRPrInDT"
result
}
#' @export
print.SimRPrInDT <- function(x,...){
## output
#anz <- length(x$depnames)
if (any(is.na(x$Struc) != TRUE)){
  cat("Best numbers of elements in substructure","\n")
  cat("2nd stage: ",x$msub[1],"\n")
  cat("3rd stage: ",x$msub[2],"\n","\n")
}
cat("Accuracies of best models","\n")
colnames(x$R2All) <- c(x$depnames,"mean")
rownames(x$R2All) <- c("full sample","indiv. opt.","joint opt.")
print(x$R2All)
cat("\n\n")
#k <- which.max(modelsI[,(anz+1)])
#print(k)
for (i in 1:length(x$depnames)){
  cat("Best model 2nd stage: ",x$depnames[i], "\n")
  print(x$modelsI[[i]])
  cat("\n")
}
cat("\n")
#k <- which.max(modelsJ[,(anz+1)])
#print(k)
for (i in 1:length(x$depnames)){
  cat("Best model after joint opt.: ",x$depnames[i], "\n")
  print(x$modelsJ[[i]])
  cat("\n")
}
}
#' @export
plot.SimRPrInDT <- function(x,...){
#anz <- length(x$depnames)
#k <- which.max(R2All[,(anz+1)])
#print(k)
for (i in 1:length(x$depnames)){
  plot(x$modelsI[[i]],main=paste0("2nd stage: ",x$depnames[i]))
}
for (i in 1:length(x$depnames)){
  plot(x$modelsJ[[i]],main=paste0("mean max.: ",x$depnames[i]))
}
}