#' Interdependent estimation for classification-regression mixtures
#'
#' @description The function SimMixPrInDT applies structured subsampling for finding an optimal subsample to model
#' the relationship between the dependent variables specified in the sublist 'targets' of the list 'datalist' and all other factor and numerical variables
#' in the corresponding individual data frame specified in the sublist 'datanames' of the list 'datalist' in the same order as 'targets'. The data frames have to be different of each other! \cr
#' The function is prepared to handle classification tasks with 2 or more classes and regression tasks. At first stage, the targets are estimated based on
#' the full sample of all exogenous variables and on summaries of the observed endogenous variables.\cr
#' For generating summaries, the variables representing the substructure have to be specified in the sublist 'datastruc' of 'datalist'.\cr
#' The sublist 'summ' of 'datalist' includes for each discrete target the classes for which you want to calculate the summary percentages and for each 
#' continuous target just NA (in the same order as in the sublist 'targets'). For a discrete target in this list, you can provide a sublist of classes to be combined 
#' (for an example see the below example).\cr
#' At second stage, structured subsampling is used for improving the models from first stage.
#' The substructure of the observations used for structured subsampling is specified by the list 'Struc' which here only consists of 
#' the name 'check' of the variable with the information about the categories of the substructure (without specification of the dataset names already specified 
#' in 'datanames', see example below), and the matrix 'labs' which specifies the values of 'check' corresponding to two categories in its rows, i.e. in 'labs[1,]' and 'labs[2,]'. 
#' The names of the categories have to be specified by  \code{rownames(labs)}.\cr
#' In structured subsampling, first 'M' repetitions of subsampling of the variable 'name' with 'nsub' different elements of each category in 'check' are realized. If 'nsub' is a list, each entry is employed individually. If 'nsub' is larger than the maximum available number of elements with a certain value of 'check', the maximum possible number of elements is used. Then,
#' for each of the subsamples 'N' repetitions of subsampling in classification or regression with the specified percentages of classes, observations, and predictors are carried out.\cr 
#' These percentages are specified in the matrix 'percent', one row per estimation task. For binary classification tasks, 
#' percentages 'percl' and 'percs' for the larger and the smaller class have to be specified. For multilevel classification tasks, NA is specified (see the below example) 
#' since the percentages are generated automatically. For regression tasks, 'pobs' and 'ppre' have to be specified for observations and predictors, respectively.\cr
#' The optimization citerion is balanced accuracy for classification and goodness of fit R2 for regression on the full sample, respectively.  
#' At stage 2, the models are optimized individually. At stage 3, the models are optimized on the maximum of joint elements in their substructures.\cr
#' The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr
#'
#' @usage SimMixPrInDT(datalist,ctestv=NA,Struc,M=12,N=99,nsub,percent=NA,conf.level=0.95,
#'                                       minsplit=NA,minbucket=NA)
#'
#' @param datalist list(datanames,targets,datastruc,summ) Input data: For specification see the above description
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param Struc list(name,check,labs) Paprametes for structured subsampling, as explained in the desciption above.
#' @param M Number of repetitions of subsampling of elements of substructure; default = 12
#' @param N Number of repetitions of subsampling for predictors (integer); default = 99
#' @param nsub (List of) numbers of different elements of substructure per subsample
#' @param percent matrix of percent spefications: For specification see the above description; default: 'percent = NA' meaning default values for percentages.
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1); default = 0.95
#' @param minsplit Minimum number of elements in a node to be splitted; default = 20
#' @param minbucket Minimum number of elements in a node; default = 7
#'
#' @return
#' \describe{
#'   \item{modelsF}{Best trees at stage 1 (Full sample)} 
#'   \item{modelsI}{Best trees at stage 2 (Individual optimization)} 
#'   \item{modelsJ}{Best trees at stage 3 (Joint optimization)} 
#'   \item{depnames}{names of dependent variables}
#'   \item{nmod}{number of models of tasks}
#'   \item{nlev}{levels of tasks}
#'   \item{accAll}{accuracies of best trees at both stages}  
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.6.2, for further information.
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{name} as well as \code{plot(name} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export SimMixPrInDT
#' @importFrom stats aggregate relevel
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
#' CHILDvowel <- as.factor(gsub("Nick","P",CHILDvowel))
#' data$Nickname <- NULL
#' syllable <- 3 - data$syllables
#' data$syllabels <- NULL
#' data$syllables <- syllable
#' data$speed <- data$word_duration / data$syllables
#' names(data)[names(data) == "target"] <- "vowel"
#' datavowel <- data
#' ##
#' # function preparation and call
#' # datalist and percent
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
#' # substructures
#' labs <- matrix(1:6,nrow=2,ncol=3)
#' labs[1,] <- c("C1a","C1b","C1c")
#' labs[2,] <- c("C2a","C2b","C2c")
#' rownames(labs) <- c("children 1","children 2")
#' Struc <- list(check="ETH",labs=labs)
#' outSimMix <- SimMixPrInDT(datalist,ctestv=ctestv,Struc=Struc,M=2,N=9,nsub=c(19,20),
#'                        percent=percent,conf.level=0.99)
#' outSimMix
#' plot(outSimMix)+
#'
SimMixPrInDT <- function(datalist,ctestv=NA,Struc,M=12,N=99,nsub,percent=NA,conf.level=0.95,minsplit=NA,minbucket=NA){
## input check
if ( typeof(datalist) != "list" || !(typeof(ctestv) %in% c("logical", "character")) 
    || !(typeof(percent) %in% c("logical","character")) || typeof(Struc) != "list" || !(typeof(M) %in% c("integer","double"))
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
#### interdependent system
#
  anz <- length(datalist$datanames)
  if (all(is.na(percent) == TRUE)){
    percent <- matrix(NA,nrow=anz,ncol=2)
  }
  nmod <- c(1:anz)
  nlev <- list()
  datat <- datalist[[1]]
  summ <- datalist[[4]]
  inds <- list()
  extdata <- list()
#####
## elements of substructure in all datasets
#####
nl <- list()
for (K in 1:anz){
  nc <- unique(datalist$datastruc[K])
  nl <- c(nl,nc)
}
ninter <- Reduce(intersect,nl)
ninter <- as.factor(ninter)
#####
## estimate with observables only
models <- list()
modmax <- list()
  for (K in 1:anz){
    data <- as.data.frame(eval(parse(text = datat[K])))   
    dataname <- datat[K]
    target <- datalist$targets[K]
    if (length(levels(data[,target])) == 2){
      nmod[K] <- 1
      nlev <- c(nlev,"class")
      sub <- unlist(datalist$datastruc[K])
      if (is.factor(sub) != TRUE){
        cat("\n Error: 'datastruc' has to contain factor variables only\n")
        return()
      }
      sub <- droplevels(sub)
      s <- unlist(str_split(summ[K],","))
      si <- which(levels(data[,target])==s)
      snotreal <- table(sub,data[,target])[,si]
      sc <- as.character(sub)
      z1 <- aggregate(data[,target] ~ sc,data=cbind(sc,data[,target]),FUN="length")
      ind <- snotreal / z1[,2]
      inds[[K]] <- ind
#      print(length(inds[[K]]))
      extdata[[K]] <- paste0(dataname,"i")
    } 
    if (length(levels(data[,target])) > 2){
      nmod[K] <- length(levels(data[,target]))
      nlev <- c(nlev,levels(data[,target]))
      sub <- unlist(datalist$datastruc[K])
      if (is.factor(sub) != TRUE){
        cat("\n Error: 'datastruc' has to contain factor variables only\n")
        return()
      }
      sub <- droplevels(sub)
      s <- unlist(str_split(summ[K],","))
      si <- which(levels(data[,target]) %in% s)
      snotmarked <- table(sub,data[,target])[,si[1]]
      if (length(s) > 1){
        for (i in 2:length(s)){
          snotmarked <- snotmarked + table(sub,data[,target])[,si[i]]
        }
      }
      sc <- as.character(sub)
      p1 <- aggregate(data[,target] ~ sc,data=cbind(sc,data[,target]),FUN="length")
      ind <- snotmarked / p1[,2]
      inds[[K]] <- ind
#      print(length(inds[[K]]))
      extdata[[K]] <- paste0(dataname,"i")
     }
     if (length(levels(data[,target])) == 0){
      nmod[K] <- 1
      nlev <- c(nlev,"real")
      sub <- unlist(datalist$datastruc[K])
      if (is.factor(sub) != TRUE){
        cat("\n Error: 'datastruc' has to contain factor variables only\n")
        return()
      }
      sc <- as.character(sub)
      v1 <- aggregate(data[,target] ~ sc,data=cbind(sc,data[,target]),FUN="mean")
      ind <- v1[,2]
      names(ind) <- v1[,1]
      inds[[K]] <- ind
#      print(length(inds[[K]]))
      extdata[[K]] <- paste0(dataname,"i")
    }
  }
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
  modelsP <- list()
  acc <- c(1:anz)*0
## first stage
  for (K in 1:anz){
    data <- eval(parse(text = extdata[K]))
    l <- colnames(data)
    if (length(unique(l)) != length(l)){
      cat("\n Error: Duplicated variable name generated: Make sure that the data frames in 'datanames' are different!\n")
      return()
    }
    dataname <- datat[K]
    target <- datalist$targets[K]
#print(levels(as.factor(data[,target])))
message(dataname,": 1st stage")
# 2 classes: call of PrInDT
    if (length(levels(data[,target])) == 2){
      out2S <- PrInDTAll(data,target,ctestv,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
      acc[K] <- out2S$baAll
      modelsP <- c(modelsP,out2S$treeAll)
    } 
# more than 2 classes: call of PrInDTMulev
    if (length(levels(data[,target])) > 2){
      out2S <- PrInDTMulevAll(data,target,ctestv,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket) # percentages chosen automatically
      acc[K] <- out2S$baAll
      modelsP <- c(modelsP,out2S$treeAll)
     }
# regression: call of PrInDTreg
     if (length(levels(data[,target])) == 0){
#       x <- str_split(percent[K,1],"=")
#       assign(x[[1]][1],eval(parse(text=percent[K,1])))
#       x <- str_split(percent[K,2],"=")
#      assign(x[[1]][1],eval(parse(text=percent[K,2])))
      out2S <- PrInDTregAll(data,target,ctestv,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)
      acc[K] <- out2S$R2All
      modelsP <- c(modelsP,out2S$treeAll)
    }
 }
accP <- acc
#print(accP)
accF <- acc
modelsF <- modelsP
#print(table(datazeroi$ETH[(!(datazeroi$MLU == 1)) & datazeroi$SEX == "male" & 
#  datazeroi$pvowel <= 168.441 & datazeroi$PRN %in% c("he","I","she","they")]))
#print(table(datazeroi$ETH[(!(datazeroi$MLU == 1)) & datazeroi$SEX == "male" & 
#  datazeroi$pvowel <= 168.441 & datazeroi$PRN %in% c("he","I","she","they") & datazeroi$zero == "zero"]))
#print(table(datazeroi$ETH[(!(datazeroi$MLU == 1)) & datazeroi$SEX == "male" & 
#  datazeroi$pvowel <= 168.441 & datazeroi$PRN %in% c("he","I","she","they") & datazeroi$zero == "real"]))
#print(table(datazeroi$ETH[(!(datazeroi$MLU == 1)) & datazeroi$SEX == "male" & 
#  datazeroi$PRN %in% c("he","I","she","they")]))
#print(table(datazeroi$ETH[(!(datazeroi$MLU == 1)) & datazeroi$SEX == "male" & 
#  datazeroi$PRN %in% c("he","I","she","they") & datazeroi$zero == "zero"]))
####################
## 2nd stage
ssize <- nsub
modmax <- modelsP
#set.seed(7654321)
 acc <- c(1:anz)*0
 NM <- N  ## NEW version
# NM <- N / M * 2 # 333
 labs <- Struc$labs
 check <- Struc$check
message("\n","Second stage")
# call of PrInDT functions
for (i in ssize){
  message("Number of elements: ",i)
  for (n in 1:M){
    message("replication ",n)
    for (K in 1:anz){
      dataname <- datat[K]
message(dataname,": 2nd stage")
      set.seed(54321*n + K*n)                                        
      datai <- eval(parse(text = extdata[K]))
      if (!(check %in% colnames(datai))){
        cat("\n Error: 'check' is not in data set\n")
        return()
      }
      z <- datai[,check]
      if (is.factor(z) != TRUE){
        cat("\n Error: 'check' has to be factor variable\n")
        return()
      }
      sub <- unlist(datalist$datastruc[K])
      sub <- droplevels(sub)
      target <- datalist$targets[K]
## 2 classes
      if (length(levels(datai[,target])) == 2){
        ind1 <- sample(unique(sub[datai[,check ]%in% labs[1,]]))[1:i] 
        ind2 <- sample(unique(sub[datai[,check] %in% labs[2,]]))[1:i] 
        datas <- datai[(as.factor(sub) %in% ind1 | as.factor(sub) %in% ind2) ,]
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
# outzero <- suppressMessages(PrInDT(datazeroi,"real",ctestv, N,percl=0.075,percs=0.9,conf.level,seedl=FALSE)) # best percentages found
        outSim <- suppressMessages(PrInDT(datas,target,ctestv, NM,percl=percl,percs=percs,conf.level,seedl=FALSE,minsplit=minsplit,minbucket=minbucket))
        class.pred <- predict(outSim$tree1st,newdata=datai)
        class.pred <- relevel(class.pred,ref=levels(datai[,target])[1])
        conti <- table(class.pred,datai[,target])
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        acc[K] <- (ba1 + ba2)/2
        models[[K]] <- outSim$tree1st
        if (acc[K] > accP[K]){ 
         accP[K] <- acc[K]
#         print(c(n,accP))
         modmax[[K]] <- models[[K]]
#         accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3  ## NOT GENERAL
       }}
## more than 2 classes
       if (length(levels(datai[,target])) > 2){
#sample(unique(name[check %in% labs[1,]]),prob=p1)[1:M] 
         ind1 <- sample(unique(sub[datai[,check ]%in% labs[1,]]))[1:i] 
         ind2 <- sample(unique(sub[datai[,check] %in% labs[2,]]))[1:i] 
         datas <- datai[(as.factor(sub) %in% ind1 | as.factor(sub) %in% ind2) ,]
         outSim <- suppressMessages(PrInDTMulev(datas,target,ctestv,NM,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)) # percentages chosen automatically
         models[[K]] <- outSim$trees
         predclass <- predict(outSim,target,newdata=datai)
         conf <- table(predclass,datai[,target]) # confusion matrix
         if (dim(conf)[1] == 1){                                                ## NEWNEW
           conf <- matrix(c(conf,0,0),nrow=2,byrow=TRUE)
         }
#           accx <- 0
#           for (i in 1:length(levels(datas[,target]))){
#             accx <- accx + conf[i,i] / colSums(conf)[i]
#           }
#           acc[K] <- accx / length(levels(datas[,target])
           accx <- rep(0,length(levels(datas[,target])))
           for (j in 1:length(rownames(conf))){
            for (l in 1:length(colnames(conf))){
              if (rownames(conf)[j] == colnames(conf)[l]){
                accx[l] <- conf[j,l] / colSums(conf)[l]
              }
           }
         }
          acc[K] <- sum(accx) / length(levels(datas[,target]))
          if (acc[K] > accP[K]) { 
            accP[K] <- acc[K]
#            print(accP)
            modmax[[K]] <- models[[K]]
#            accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3  ## not general
          }
       }
## regression
       if (length(levels(datai[,target])) == 0){
#sample(unique(name[check %in% labs[1,]]),prob=p1)[1:M] 
         ind1 <- sample(unique(sub[datai[,check ]%in% labs[1,]]))[1:i] 
         ind2 <- sample(unique(sub[datai[,check] %in% labs[2,]]))[1:i] 
         datas <- datai[(as.factor(sub) %in% ind1 | as.factor(sub) %in% ind2) ,]
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
         outSim <- suppressMessages(PrInDTreg(datas,target,ctestv,NM,pobs,ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket))
         models[[K]] <- outSim$ctmax
         ctpreds <- predict(outSim$ctmax,newdata=datai)
         acc[K] <- 1 - sum((datai[,target] - ctpreds)^2) / sum((datai[,target]-mean(datai[,target]))^2) ## R2 !!!!!
         if (acc[K] > accP[K]) { 
           accP[K] <- acc[K]
#           print(accP)
           modmax[[K]] <- models[[K]]
#           accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3       ## not general
         }
       }          
}}}
#print(accP)
accI <- accP
accF <- cbind(accF,accI)
modelsI <- modmax
#plot(unlist(modmax[[1]]))
#plot(unlist(modmax[[2]][[1]]))
#plot(unlist(modmax[[2]][[2]]))
#plot(unlist(modmax[[2]][[3]]))
#plot(unlist(modmax[[3]]))
##############################################
# call on subsamples: joined version
##
#ssize <- length(ninter)/2
set.seed(7654321)
acc <- c(1:anz)*0
#accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3  ## not general
n <- 1
message("\n","Third stage: Joint version")
message(length(ninter)," elements in substructure")
#for (i in ssize){
#  print(i)
#  for (n in 1:M){
#    print(n)
    for (K in 1:anz){  
      dataname <- datat[K]
      set.seed(54321*n + K*n)                                        
      datai <- eval(parse(text = extdata[K]))
      sub <- unlist(datalist$datastruc[K])
      sub <- droplevels(sub)
      if (K == 1){
        ind1 <- sample(unique(sub[sub %in% ninter & datai[,check ] %in% labs[1,]]))  ## [1:i] 
        ind2 <- sample(unique(sub[sub %in% ninter & datai[,check ] %in% labs[2,]]))  ## [1:i] 
        ind1 <- droplevels(ind1)
        ind2 <- droplevels(ind2)
        lind1 <- length(ind1)
        lind2 <- length(ind2)
message("number of elements in categories: ",lind1," , ",lind2)
     }
      target <- datalist$targets[K]
message(dataname,": joint version")
## 2 classes
      if (length(levels(datai[,target])) == 2){
#        ind1 <- sample(unique(sub[sub %in% ninter & datai[,check ] %in% labs[1,]]))  ## [1:i] 
#        ind2 <- sample(unique(sub[sub %in% ninter & datai[,check] %in% labs[2,]]))   ## [1:i] 
        datas <- datai[(as.factor(sub) %in% ind1) | (as.factor(sub) %in% ind2) ,]
#if (K == 1){
#print(table(datas$ETH[(!(datas$MLU == 1)) & datas$SEX == "male" & 
#  datas$pvowel <= 168.441 & datas$PRN %in% c("he","I","she","they") & datas$zero == "zero"]))
#}
# assigning percl, percs
        percl <- NA
        percs <- 1
        x <- str_split(percent[K,1],"=")
        assign(x[[1]][1],eval(parse(text=percent[K,1])))
        x <- str_split(percent[K,2],"=")
        assign(x[[1]][1],eval(parse(text=percent[K,2])))
        outSim <- suppressMessages(PrInDT(datas,target,ctestv, N,percl=percl,percs=percs,conf.level,seedl=FALSE,minsplit=minsplit,minbucket=minbucket)) 
        class.pred <- predict(outSim$tree1st,newdata=datai)
        class.pred <- relevel(class.pred,ref=levels(datai[,target])[1])
        conti <- table(class.pred,datai[,target])
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        acc[K] <- (ba1 + ba2)/2
        models[[K]] <- outSim$tree1st
#        if (acc[K] > accP[K]){ 
#         accP[K] <- acc[K]
#         print(c(n,accP))
#         modmax[[K]] <- models[[K]]
#         accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3  ## NOT GENERAL
#       }
      }
## more than 2 classes
       if (length(levels(datai[,target])) > 2){
#         ind1 <- sample(unique(sub[sub %in% ninter & datai[,check ]%in% labs[1,]]))  ## [1:i] 
#         ind2 <- sample(unique(sub[sub %in% ninter & datai[,check] %in% labs[2,]]))  ## [1:i] 
         datas <- datai[(as.factor(sub) %in% ind1) | (as.factor(sub) %in% ind2) ,]
         outSim <- suppressMessages(PrInDTMulev(datas,target,ctestv,N,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket)) # percentages chosen automatically
         models[[K]] <- outSim$trees
         predclass <- predict(outSim,target,newdata=datai)
         conf <- table(predclass,datai[,target]) # confusion matrix
         if (dim(conf)[1] == 1){                                                ## NEWNEW
           conf <- matrix(c(conf,0,0),nrow=2,byrow=TRUE)
         }
#print(conf)
#           accx <- 0
#           for (i in 1:length(levels(datas[,target]))){
#             accx <- accx + conf[i,i] / colSums(conf)[i]
#           }
#           acc[K] <- accx / length(levels(datas[,target])
           accx <- rep(0,length(levels(datas[,target])))
           for (j in 1:length(rownames(conf))){
            for (l in 1:length(colnames(conf))){
              if (rownames(conf)[j] == colnames(conf)[l]){
                accx[l] <- conf[j,l] / colSums(conf)[l]
              }
           }
         }
          acc[K] <- sum(accx) / length(levels(datas[,target]))
#          if (acc[K] > accP[K]) { 
#            accP[K] <- acc[K]
#            print(accP)
#            modmax[[K]] <- models[[K]]
#            accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3  ## not general
#          }
       }
## regression
       if (length(levels(datai[,target])) == 0){
#         ind1 <- sample(unique(sub[sub %in% ninter & datai[,check ] %in% labs[1,]]))  ## [1:i] 
#         ind2 <- sample(unique(sub[sub %in% ninter & datai[,check] %in% labs[2,]]))  ## [1:i] 
         datas <- datai[(as.factor(sub) %in% ind1) | (as.factor(sub) %in% ind2) ,]
# assigning pobs, ppre
         pobs <- c(0.9,0.7)
         ppre <- c(0.9,0.7)
         x <- str_split(percent[K,1],"=")
         assign(x[[1]][1],eval(parse(text=percent[K,1])))
         x <- str_split(percent[K,2],"=")
         assign(x[[1]][1],eval(parse(text=percent[K,2]))) 
         NM <- N*M/2  ## NEW version
         outSim <- suppressMessages(PrInDTreg(datas,target,ctestv,NM,pobs,ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket))
         models[[K]] <- outSim$ctmax
         ctpreds <- predict(outSim$ctmax,newdata=datai)
         acc[K] <- 1 - sum((datai[,target] - ctpreds)^2) / sum((datai[,target]-mean(datai[,target]))^2) ## R2 !!!!!
#         if (acc[K] > accP[K]) { 
#           accP[K] <- acc[K]
#           print(accP)
#           modmax[[K]] <- models[[K]]
#           accPmean <- (accP[1] + accP[2] + sqrt(accP[3]))/3       ## not general
#         }
       }          
}#}
#accmean <- (acc[1] + acc[2] + sqrt(acc[3]))/3  ## not general
#  if (accmean > accPmean) { 
#    accP[1:anz] <- acc[1:anz]
#    print(accP)
#    for (i in 1:anz){
#      modmax[[i]] <- models[[i]]
#    }
#   accPmean <- accmean
#  }
#print(acc)
#print(accP)
#print(accPmean)
accJ <- acc
accF <- cbind(accF,acc)
modelsJ <- models
#plot(modelsJ[[1]])
#plot(modelsJ[[2]][[1]])
#plot(modelsJ[[2]][[2]])
#plot(modelsJ[[2]][[3]])
#plot(modelsJ[[3]])
result <- list(modelsF=modelsF,modelsI=modelsI,modelsJ=modelsJ,depnames=unlist(datalist$target),nmod=nmod,nlev=nlev,accAll=accF)
class(result) <- "SimMixPrInDT"
result
}
#' @export
print.SimMixPrInDT <- function(x,...){
## output
cat("Accuracies of best models","\n")
rownames(x$accAll) <- c(x$depnames)
colnames(x$accAll) <- c("1st stage","2nd stage","3rd stage")
print(x$accAll)
cat("\n\n")
J <- 0
for (i in 1:length(x$depnames)){
  if (x$nmod[i] == 1){
    J <- J + 1
    cat("Best model 2nd stage: ",x$depnames[i], "\n") 
    print(x$modelsI[[i]])
    cat("\n\n")
  } else {
    K <- x$nmod[i]
    for (k in 1:K){
       J <- J + 1
       cat("Best model 2nd stage: ",x$depnames[i],":",x$nlev[[J]],"\n") 
       print(x$modelsI[[i]][[k]])
       cat("\n\n")
    }
  }
}}
#' @export
plot.SimMixPrInDT <- function(x,...){
J <- 0
for (i in 1:length(x$depnames)){
  if (x$nmod[i] == 1){
    J <- J + 1
    plot(x$modelsI[[i]],main=paste0("2nd stage: ",x$depnames[i]))
  } else {
    K <- x$nmod[i]
    for (k in 1:K){
       J <- J + 1
       plot(x$modelsI[[i]][[k]],main=paste0("2nd stage: ",x$depnames[i],":",x$nlev[[J]]))
    }
  }
}}