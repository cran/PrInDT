#' Structured subsampling for regression
#'
#' @description The function PrInDTRstruc applies structured subsampling for finding an optimal subsample to model
#' the relationship between the continuous variable 'regname' and all other factor and numerical variables
#' in the data frame 'datain' by means of 'M' repetitions of subsampling from a substructure and 'N' repetitions of subsampling from the predictors. 
#' The optimization citerion is the goodness of fit R2 on the validation sample 'valdat' (default = 'datain').  
#' The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The substructure of the observations used for subsampling is specified by the list 'Struc' which consists of the 'name' of the variable representing the substructure,
#' the name 'check' of the variable with the information about the categories of the substructure, and the matrix 'labs' which specifies the values of 'check'
#' corresponding to two categories in its rows, i.e. in 'labs[1,]' and 'labs[2,]'. The names of the categories have to be specified by \code{rownames(labs)}.\cr
#' The number of predictors 'Pit' to be included in the model and the number of elements of the substructure 'Mit' have to be specified (lists allowed).\cr
#' The percentages of involved observations and predictors can be controlled by the parameters 'pobs' and 'ppre', respectively.\cr
#' The parameter 'Struc' is needed for all versions of subsampling except "b". Four different versions of structured subsampling exist:\cr 
#' a) just of the elements in the substructure with parameters 'M' and 'Mit',\cr
#' b) just of the predictors with parameters 'N' and 'Pit',\cr
#' c) of the predictors and for each subset of predictors subsampling of the elements of the substructure with parameters 'M', 'N', 'Mit', 'Pit', 'pobs', and 'ppre', and\cr
#' d) of the elements of the substructure and for each of these subsets subsampling of the predictors with the same parameters as version c).\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr\cr
#' Repeated measurements can also be handled by this function (indrep=1). They are multiple measurements of the same variable taken on the same subjects (or objects) either under different conditions 
#' or over two or more time periods.\cr
#' The name of the variable with the repeatedly observed subjects (or objects) has to be specified by 'name' in 'Struc'.\cr 
#' Additionally, in 'check' and 'labs' in 'Struc' a binary variable has to be specified for which the same number of 
#' subjects (or objects) of the two classes of this variable should be used for model building. Possible such numbers can be specified by 'Mit'.  
#'
#' @usage PrInDTRstruc(datain,regname,ctestv=NA,Struc=NA,vers="d",M=NA,Mit=NA,N=99,Pit=NA,
#'                pobs=c(0.9,0.7),ppre=c(0.9,0.7),conf.level=0.95,minsplit=NA,minbucket=NA,
#'                valdat=datain,indrep=0)
#'
#' @param datain Input data frame with continuous target variable 'regname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param regname Name of target variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in variable 1 in a tree.\cr
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     If no restrictions exist, the default = NA is used.
#' @param Struc = list(name,check,labs), cf. description for explanations; Struc not needed for vers="b"
#' @param vers Version of structured subsampling: "a", "b", "c", "d", cf. description; \cr
#'     default = "d"
#' @param Mit List of number of elements of substructure (integers); \cr
#'     default = c((Cl-4):Cl), Cl = maximum number elements in both categories 
#' @param Pit List of number of predictors (integers)\cr
#'     default = c(max(1,(D-3)):D), D = maximum number of predictors
#' @param M Number of repetitions of subsampling from substructure (integer) in versions "a" and "d";\cr
#'      default = 99
#' @param N Number of repetitions of subsampling from predictors (integer) in versions "b" and "c";\cr
#'      default = 99
#' @param pobs Percentage(s) of observations for subsampling in versions "c" and "d";\cr
#'      default=c(0.9,0.7)
#' @param ppre Percentage(s) of predictors for subsampling in versions "c" and "d";\cr
#'      default=c(0.9,0.7)
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#' @param valdat Validation data; default = datain
#' @param indrep Indicator for repeated measurements, i.e. more than one observation with the same class for each element;\cr 
#'     indrep=1: Struc=list(name) only; default = 0 
#'
#' @return
#' \describe{
#'   \item{outmax}{Best tree} 
#'   \item{interp}{Number of interpretable trees, overall number of trees} 
#'   \item{ntmax}{Size of training set for best tree} 
#'   \item{R2max}{R squared of best tree}
#'   \item{R2sub}{Mean R squared of objects in substructure}
#'   \item{MAEmax}{MAE (Mean Absolute Error) of best tree}
#'   \item{MAEsub}{Mean MAE of objects in substructure}
#'   \item{ind1max}{Elements of 1st category of substructure used by best tree} 
#'   \item{ind2max}{Elements of 2nd category of substructure used by best tree} 
#'   \item{indmax}{Predictors used by best tree} 
#'   \item{gmaxTrain}{Training set for best tree} 
#'   \item{labs}{labs from Struc} 
#'   \item{vers}{Version used for structured subsampling} 
#'   \item{lMit}{Number of different numbers of substructure elements}
#'   \item{lPit}{Number of different numbers of predictors}
#'   \item{M}{Number of repetitions of selection of substructure elements} 
#'   \item{N}{Number of repetitions of selection of predictors}
#'   \item{indrep}{Indicator of repeated measurements: indrep=1} 
#' }
#'
#' @details
#' See Buschfeld & Weihs (2025), Optimizing decision trees for the analysis of World Englishes and sociolinguistic data. Cambridge University Press, section 4.5.4, for further information.
#'
#' Standard output can be produced by means of \code{print(name$besttree)} or just \code{name$besttree} as well as \code{plot(name$besttree)} where 'name' is the output data 
#' frame of the function.\cr
#
#' @export PrInDTRstruc
#'
#' @examples
#' data <- PrInDT::data_vowel
#' data <- na.omit(data)
#' CHILDvowel <- data$Nickname
#' data$Nickname <- NULL
#' data$syllables <- 3 - data$syllables
#' data$speed <- data$word_duration / data$syllables  ## NEW NEW
#' names(data)[names(data) == "target"] <- "vowel_length"
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}','MLU == {1, 3}') # split exclusions
#' name <- CHILDvowel
#' check <- "data$ETH"
#' labs <- matrix(1:6,nrow=2,ncol=3)
#' labs[1,] <- c("C1a","C1b","C1c")
#' labs[2,] <- c("C2a","C2b","C2c")
#' rownames(labs) <- c("children 1","children 2")
#' Struc <- list(name=name,check=check,labs=labs)
#' outstruc <- PrInDTRstruc(data,"vowel_length",ctestv=ctestv,Struc=Struc,vers="d",
#'                   M=3,Mit=21,N=9,pobs=c(0.95,0.7),ppre=c(1,0.7),conf.level=0.99)
#' outstruc
#' plot(outstruc)
#'
PrInDTRstruc <- function(datain,regname,ctestv=NA,Struc=NA,vers="d",M=NA,Mit=NA,N=99,Pit=NA,pobs=c(0.9,0.7),ppre=c(0.9,0.7),conf.level=0.95,
                                         minsplit=NA,minbucket=NA,valdat=datain,indrep=0){
## input check
if ( typeof(datain) != "list" || typeof(regname) != "character" || !(typeof(ctestv) %in% c("logical", "character"))
    || !(typeof(Struc) %in% c("logical","list")) || typeof(vers) != "character" || !all(0 < pobs & pobs <= 1) || !all(0 < ppre & ppre <= 1) 
    || !(typeof(Mit) %in% c("logical","integer","double")) || !(typeof(Pit) %in% c("logical","integer","double")) || !(typeof(N) %in% c("logical","integer","double")) 
    || !(typeof(M) %in% c("logical","integer","double")) || !(0 < conf.level & conf.level <= 1) || !(typeof(minsplit) %in% c("logical","double")) 
    || !(typeof(minbucket) %in% c("logical", "double")) || typeof(valdat) != "list" || !(typeof(indrep) %in% c("integer","double")) ) {
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
data <- datain
if ((vers == "a") | (vers == "c") | (vers == "d")){ 
  name <- Struc$name
#  check <- Struc$check
# print(names(data))
  check <- eval(parse(text = Struc$check))
  if (length(check) == 0){
    stop("irregular input: check variable not in data set")
  }
  labs <- Struc$labs
  if (all(is.na(Mit)) == TRUE){stop("No number of elements in substructure specified")}
  lMit <- length(Mit)
} else {
  lMit <- NA
  gmaxTrain <- NA
}
##
if ((vers == "a") | (vers == "d")){
  R2max <- 0
  interp <- 0
  indmax <- NA
  set.seed(7654321)
  for (Mi in Mit){
#    set.seed(7654321)
    message("Number of substructure elements per category: ",Mi)
    for (i in 1:M){
      ind1 <- sample(unique(name[check %in% labs[1,]]))[1:Mi]
#print(ind1)
      ind2 <- sample(unique(name[check %in% labs[2,]]))[1:Mi]
#print(ind2)
      gTrain <- data[(name %in% ind1 | name %in% ind2),]
      gTest <- data[!(name %in% ind1 | name %in% ind2),]
      gcheck <- check[(name %in% ind1 | name %in% ind2)]
      gname <- name[(name %in% ind1 | name %in% ind2)]
      nt <- nrow(gTrain)
      if (vers == "d"){
#print(dim(gTrain))
        suppressMessages(outAll <- PrInDTreg(gTrain,regname,ctestv=ctestv,N,pobs=pobs,ppre=ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket) )
        if (outAll$interpmax == FALSE){
          interp <- interp + 1
        }
        ctpreds <- predict(outAll$ctmax,newdata=valdat)
        R2 <- 1 - sum((valdat[,regname] - ctpreds)^2) / sum((valdat[,regname]-mean(valdat[,regname]))^2)
        MAE <- sum(abs(valdat[,regname] - ctpreds))/ dim(valdat)[1]
#print(R2)
        if (R2 > R2max & outAll$interpmax == FALSE){
          R2max <- R2
          MAEmax <- MAE
#print(c(i,R2max))
          ntmax <- nt 
          ind1max <- ind1
          ind2max <- ind2
          outmax <- outAll
          gmaxTrain <- cbind(gTrain,as.vector(gname),as.vector(gcheck))
          colnames(gmaxTrain) <- c(colnames(gTrain),"SubStruc","check")
          gmaxTest <- gTest
        }
      }
      if (vers == "a"){
        suppressMessages(outAll <- PrInDTregAll(gTrain,regname,ctestv=ctestv,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket) )
        if (outAll$interpAll == FALSE){
          interp <- interp + 1
        }
        ctpreds <- predict(outAll$treeAll,newdata=data)
        R2 <- 1 - sum((data[,regname] - ctpreds)^2) / sum((data[,regname]-mean(data[,regname]))^2)
        MAE <- sum(abs(data[,regname] - ctpreds))/ dim(data)[1]
        if (R2 > R2max & outAll$interpAll == FALSE){
          R2max <- R2
          MAEmax <- MAE
#print(c(i,R2max))
          ntmax <- nt 
          ind1max <- ind1
          ind2max <- ind2
          outmax <- outAll
          gmaxTrain <- cbind(gTrain,as.vector(gname),as.vector(gcheck))
          colnames(gmaxTrain) <- c(colnames(gTrain),"SubStruc","check")
          gmaxTest <- gTest
        }
      }
    }
message("current optimal R2 value: ",round(R2max,digits=4))
  }
}
##
D <- dim(data)[2] - 1 
if (all(is.na(Pit)) == TRUE){
  Pit <- c(D-3,D-2,D-1)
}
lPit <- length(Pit)
if ((vers == "b") | (vers == "c")){
  interp <- 0
  R2max <- 0
  set.seed(7654321)
  for (J in Pit) {
    message("number of predictors: ",J)
    if (vers == "b"){
      for (n in 1:N){
        datat <- data[,names(data)!=regname] 
        ind <- sample(1:D)[1:J] 
#print(ind)
        datat <- datat[,ind]
        datat <- cbind(datat,data[,regname])
        colnames(datat)[J+1] <- regname
        suppressMessages(outAll <- PrInDTregAll(datat,regname,ctestv,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket))
        if (outAll$interpAll == FALSE){
          interp <- interp + 1
        }
        ctpreds <- predict(outAll$treeAll,newdata=data)
        R2 <- 1 - sum((data[,regname] - ctpreds)^2) / sum((data[,regname]-mean(data[,regname]))^2)
        MAE <- sum(abs(data[,regname] - ctpreds))/ dim(data)[1]
        if (R2 > R2max & outAll$interpAll == FALSE){
          outmax <- outAll
          indmax <- ind
          R2max <- R2
          MAEmax <- MAE
#print(c(J,n,R2max))
        }
      }
      ntmax <- nrow(datat)
      ind1max <- NA
      ind2max <- NA
    }
    if (vers == "c"){
       for (i in Mit){
#         for (m in 1:N){
           datat <- data[,names(data)!=regname] 
           ind <- sample(1:D)[1:J] 
           datat <- datat[,ind]
           datat <- cbind(datat,data[,regname])
           colnames(datat)[J+1] <- regname
           ind1 <- sample(unique(name[check %in% labs[1,]]))[1:i] # select observations
           ind2 <- sample(unique(name[check %in% labs[2,]]))[1:i]
           gTrain <- datat[(name %in% ind1 | name %in% ind2),]
           gname <- name[(name %in% ind1 | name %in% ind2)]
           gcheck <- check[(name %in% ind1 | name %in% ind2)]
           nt <- nrow(gTrain)
          suppressMessages(outAll <- PrInDTreg(gTrain,regname,ctestv=ctestv,N,pobs=pobs,ppre=ppre,conf.level=conf.level,minsplit=minsplit,minbucket=minbucket) )
          if (outAll$interpmax == FALSE){
            interp <- interp + 1
          }
          ctpreds <- predict(outAll$ctmax,newdata=valdat)
          R2 <- 1 - sum((valdat[,regname] - ctpreds)^2) / sum((valdat[,regname]-mean(valdat[,regname]))^2)
          MAE <- sum(abs(valdat[,regname] - ctpreds))/ dim(valdat)[1]
          if (R2 > R2max & outAll$interpmax == FALSE){
            outmax <- outAll
            indmax <- ind
            R2max <- R2
            MAEmax <- MAE
            ntmax <- nt 
            ind1max <- ind1
            ind2max <- ind2
            gmaxTrain <- cbind(gTrain,as.vector(gname),as.vector(gcheck))
            colnames(gmaxTrain) <- c(colnames(gTrain),"SubStruc","check")
#print(c(J,i,R2max))
          }
#        }
      }
    }
message("current optimal R2 value: ",round(R2max,digits=4))
  }
}
R2sub <- 0
MAEsub <- 0
if ((vers == "a") | (vers == "c") | (vers == "d")){ 
  pred <- predict(outmax,newdata=data)
  diff <- abs(pred - data[,regname])
  agg <- aggregate(diff,list(name=name),mean)
  MAEsub <- mean(agg[,2])
##
  diff <- (pred - data[,regname])^2
  agg <- aggregate(diff,list(name=name),sum)
  diff <- (data[,regname]-mean(data[,regname]))^2
  agg <- 1 - agg[,2] / aggregate(diff,list(name=name),sum)[,2]
  R2sub <- mean(agg)
#  meansub <- vector()
#  aggsub <- vector()
#  diffsub <- matrix(0,nrow=dim(valdat)[1],ncol=length(unique(name)))
#  for (l in 1:length(unique(name))){
#    meansub[l] <- mean(valdat[as.integer(name)==l,regname])
#    diffsub[1:length(valdat[as.integer(name)==l,regname]),l] <- (valdat[as.integer(name)==l,regname] - meansub[l])^2
#    aggsub[l] <- sum(diffsub[,l])
#  }
#print(agg[,2])
#print(aggsub)  ## much smaller: WHY??
#  agg <- agg[,2] / aggsub
#  R2sub <- 1 - mean(agg)
}
## 
## results
result <- list(outmax=outmax,interp=interp,ind1max=ind1max,ind2max=ind2max,indmax=indmax,vers=vers,R2max=R2max,MAEmax=MAEmax,
     ntmax=ntmax,lMit=lMit,lPit=lPit,M=M,N=N,gmaxTrain=gmaxTrain,labs=labs,MAEsub=MAEsub,R2sub=R2sub,indrep=indrep)
class(result) <- "PrInDTRstruc"
result
}
#' @export
print.PrInDTRstruc <- function(x,...){
  if ((x$vers == "d") | (x$vers == "a")){
    cat("Number of interpretable best trees: ",x$interp," of ",(x$lMit*x$M)," trees","\n","\n")
#    cat(x$interp," of ",(x$lMit*x$M)," trees","\n","\n")
    cat("Best interpretable tree from subsampling","\n")
    if (x$vers == "d") {print(x$outmax$ctmax)}
    if (x$vers == "a") {print(x$outmax$treeAll)}
##plot(outmax)
    cat("\n","Best R2 on validation data: ",x$R2max)
    cat("\n","Best MAE on validation data: ",x$MAEmax,"\n")
    cat("\n","Size of full training set for best tree: ",x$ntmax,"\n")
    cat("\n","Number of observations of ",rownames(x$labs)[1]," in model construction: ",sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)] %in% x$labs[1,]),"\n")
#    cat(sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)] %in% x$labs[1,]),"\n")
    Elements <- cbind(as.character(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)]),
      as.character(x$gmaxTrain$SubStruc[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)]))
    print(table(Elements,exclude=x$labs[1,]))
    print(table(Elements,exclude=x$ind1max))
    cat("\n","Number of observations of ",rownames(x$labs)[2]," in model construction: ",sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)] %in% x$labs[2,]),"\n")
#    cat(sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)] %in% x$labs[2,]),"\n")
    Elements <- cbind(as.character(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)]),
      as.character(x$gmaxTrain$SubStruc[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)]))
    print(table(Elements,exclude=x$labs[2,]))
    print(table(Elements,exclude=x$ind2max))
#cat("\n","Test accuracy of best tree","\n")
#ctpreds <- predict(outmax$treeAll,newdata=gmaxTest)
#print(1 - sum((data[,regname] - ctpreds)^2) / sum((data[,regname]-mean(data[,regname]))^2))
#cat("Accuracy of best tree")
#ctpreds <- predict(outmax$treeAll,newdata=data)
#print(1 - sum((data[,regname] - ctpreds)^2) / sum((data[,regname]-mean(data[,regname]))^2))
#cat("Accuracies of best tree","\n"," ",lablarge," ",labsmall," balanced","\n")
#cat(" ",c(acc1,acc2,bamax))
  }
  if ((x$vers == "b")){
    cat("Number of interpretable best trees: ",x$interp," of ",(x$N*x$lPit)," trees","\n","\n")
#    cat(x$interp," of ",(x$N*x$lPit)," trees","\n","\n") 
    cat("Best interpretable tree","\n")
    print(x$outmax$treeAll)
##    plot(outmax)
    cat("\n","Best R2 on validation data: ",x$R2max)
    cat("\n","Best MAE on validation data: ",x$MAEmax,"\n")
  }
  if ((x$vers == "c")){
    cat("Number of interpretable best trees: ",x$interp," of ",(x$N*x$lPit*x$lMit)," trees","\n","\n")
#    cat(x$interp," of ",(x$N*x$Pit*x$lMit)," trees","\n","\n")
    cat("Best interpretable tree","\n")
    print(x$outmax$ctmax)
##    plot(outmax)
    cat("\n","Best R2 on validation data: ",x$R2max)
    cat("\n","Best MAE on validation data: ",x$MAEmax,"\n")
    cat("\n","Size of training set for best tree: ",x$ntmax,"\n")
#    cat(x$ntmax,"\n")
    cat("\n","Tokens of ",rownames(x$labs)[1]," in model construction: ",sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)] %in% x$labs[1,]),"\n")
#    cat(sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)] %in% x$labs[1,]),"\n")
    Elements <- cbind(as.character(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)]),
      as.character(x$gmaxTrain$SubStruc[x$gmaxTrain$SubStruc %in% droplevels(x$ind1max)]))
    print(table(Elements,exclude=x$labs[1,]))
    print(table(Elements,exclude=x$ind1max))
    cat("\n","Tokens of ",rownames(x$labs)[2]," in model construction: ",sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)] %in% x$labs[2,]),"\n")
#    cat(sum(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)] %in% x$labs[2,]),"\n")
    Elements <- cbind(as.character(x$gmaxTrain$check[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)]),
      as.character(x$gmaxTrain$SubStruc[x$gmaxTrain$SubStruc %in% droplevels(x$ind2max)]))
    print(table(Elements,exclude=x$labs[2,]))
    print(table(Elements,exclude=x$ind2max))
  }
  if (x$indrep == 1){
    if ((x$vers == "a") | (x$vers == "c") | (x$vers == "d")){
      cat("\n","Repeated measurements","\n") 
      cat("Measures for objects in substructure","\n")
#    cat("Mean R2 of objects in substructure: ",x$R2sub,"\n")
      cat("Mean MAE on full original data: ",x$MAEsub,"\n\n")
    }
  }
}
#' @export
plot.PrInDTRstruc <- function(x,...){
  plot(x$outmax,main="Best interpretable tree")
}