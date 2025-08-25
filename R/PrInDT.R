#' The basic undersampling loop for classification
#'
#' @description The function PrInDT uses ctrees (conditional inference trees from the package "party") for optimal modeling of
#' the relationship between the two-class factor variable 'classname' and all other factor and numerical variables
#' in the data frame 'datain' by means of 'N' repetitions of undersampling. The optimization citerion is the balanced accuracy 
#' on the validation sample 'valdat' (default = full input sample 'datain'). The trees generated from undersampling can be restricted by not accepting trees 
#' including split results specified in the character strings of the vector 'ctestv'.\cr
#' The undersampling percentages are 'percl' for the larger class and 'percs' for the smaller class (default = 1).\cr
#' A probability threshold 'thres' for the prediction of the smaller class may be specified (default = 0.5).\cr
#' Undersampling may be stratified in two ways by the option 'stratvers'.\cr
#' The parameters 'conf.level', 'minsplit', and 'minbucket' can be used to control the size of the trees.\cr\cr
#' In the case of repeated measurements ('indrep=1'), the values of the substructure variable have to be given in 'repvar'. 
#' Only one value of 'classname' is allowed for each value of 'repvar'.
#' If for a value of 'repvar' the percentage 'thr' of the observed occurence of a value of 'classname' is not reached by the number of predictions of the value of 'classname', a misclassification is detected.
#'
#' @usage PrInDT(datain,classname,ctestv=NA,N,percl,percs=1,conf.level=0.95,thres=0.5,
#'          stratvers=0,strat=NA,seedl=TRUE,minsplit=NA,minbucket=NA,repvar=NA,indrep=0,
#'          valdat=datain,thr=0.5)
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
#' @param N Number (> 2) of repetitions (integer)
#' @param percl Undersampling percentage of larger class (numerical, > 0 and <= 1);\cr
#'     if percs = default, default of percl = percentage of large class resulting in the same number of observations as in the small class
#' @param percs Undersampling percentage of smaller class (numerical, > 0 and <= 1);\cr
#'     default = 1
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#' @param thres Probability threshold for prediction of smaller class (numerical, >= 0 and < 1); default = 0.5
#' @param stratvers Version of stratification;\cr
#'     = 0: none (default),\cr
#'     = 1: stratification according to the percentages of the values of the factor variable 'strat',\cr
#'     > 1: stratification with minimum number "stratvers" of observations per value of "strat"
#' @param strat Name of one (!) stratification variable for undersampling (character);\cr
#'     default = NA (no stratification)
#' @param seedl Should the seed for random numbers be set (TRUE / FALSE)?\cr
#'     default = TRUE
#' @param minsplit Minimum number of elements in a node to be splitted;\cr
#'     default = 20
#' @param minbucket Minimum number of elements in a node;\cr
#'     default = 7
#' @param repvar Values of variable defining the substructure in the case of repeated measurements, length = dim(valdat)[1] necessary; default = NA
#' @param indrep Indicator of repeated measurements ('indrep=1'); default = 0
#' @param valdat Validation data; default = datain
#' @param thr threshold for element classification: minimum percentage of correct class entries; default = 0.5
#'
#' @return
#' \describe{
#' \item{tree1st}{best tree on validation sample}
#' \item{tree2nd}{2nd-best tree on validation sample}
#' \item{tree3rd}{3rd-best tree on validation sample}
#' \item{treet1st}{best tree on test sample}
#' \item{treet2nd}{2nd-best tree on test sample}
#' \item{treet3rd}{3rd-best tree on test sample}
#' \item{ba1st}{accuracies: largeClass, smallClass, balanced of 'tree1st', both for validation and test sample}
#' \item{ba2nd}{accuracies: largeClass, smallClass, balanced of 'tree2nd', both for validation and test sample}
#' \item{ba3rd}{accuracies: largeClass, smallClass, balanced of 'tree3rd', both for validation and test sample}
#' \item{baen}{accuracies: largeClass, smallClass, balanced of ensemble of all interpretable, 3 best acceptable, and all acceptable trees on validation sample}
#' \item{bafull}{vector of balanced accuracies of all trees from undersampling}
#' \item{batest}{vector of test accuracies of all trees from undersampling}
#' \item{treeAll}{tree based on all observations}
#' \item{baAll}{balanced accuracy of 'treeAll'}
#' \item{interpAll}{criterion of interpretability of 'treeall' (TRUE / FALSE)}
#' \item{confAll}{confusion matrix of 'treeAll'}
#' \item{acc1AE}{Accuracy of full sample tree on Elements of large class} 
#' \item{acc2AE}{Accuracy of full sample tree on Elements of small class }  
#' \item{bamaxAE}{balanced accuracy of full sample tree on Elements} 
#' \item{namA1}{Names of misclassified Elements by full sample tree of large class}
#' \item{namA2}{Names of misclassified Elements by full sample tree of small class}
#' \item{acc1E}{Accuracy of best tree on Elements of large class} 
#' \item{acc2E}{Accuracy of best tree on Elements of small class }  
#' \item{bamaxE}{balanced accuracy of best tree on Elements} 
#' \item{nam1}{Names of misclassified Elements by best tree of large class}
#' \item{nam2}{Names of misclassified Elements by best tree of small class}
#' \item{lablarge}{Label of large class} 
#' \item{labsmall}{Label of small class} 
#' \item{valdat}{Validation set after edit}
#' \item{thr}{Threshold for repeated measurements}
#' }
#'
#' @details
#' For the optimzation of the trees, we employ a method we call Sumping (Subsampling umbrella of 
#' model parameters), a variant of Bumping (Bootstrap umbrella of model parameters) (Tibshirani 
#' & Knight, 1999) which uses subsampling instead of bootstrapping. The aim of the 
#' optimization is to identify conditional inference trees with maximum predictive power
#' on the full sample under interpretability restrictions.
#'
#' \strong{References} \cr
#' -- Tibshirani, R., Knight, K. 1999. Model Search and Inference By Bootstrap "bumping".
#' Journal of Computational and Graphical Statistics, Vol. 8, No. 4 (Dec., 1999), pp. 671-686 \cr
#' -- Weihs, C., Buschfeld, S. 2021a. Combining Prediction and Interpretation in  Decision Trees (PrInDT) - 
#' a Linguistic Example. arXiv:2103.02336
#'
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.\cr
#' The plot function will produce a series of more than one plot. If you use R, you might want to specify \code{windows(record=TRUE)} before 
#' \code{plot(name)} to save the whole series of plots. In R-Studio this functionality is provided automatically.
#
#' @export PrInDT
#'
#' @examples
#' datastrat <- PrInDT::data_zero
#' data <- na.omit(datastrat) # cleaned full data: no NAs
#' # interpretation restrictions (split exclusions)
#' ctestv <- rbind('ETH == {C2a, C1a}','MLU == {1, 3}') # split exclusions
#' N <- 41  # no. of repetitions
#' conf.level <- 0.99 # 1 - significance level (mincriterion) in ctree
#' percl <- 0.08  # undersampling percentage of the larger class
#' percs <- 0.95 # undersampling percentage of the smaller class
#' # calls of PrInDT
#' out <- PrInDT(data,"real",ctestv,N,percl,percs,conf.level) # unstratified
#' out # print best model and ensembles as well as all observations
#' plot(out)
#' out <- PrInDT(data,"real",ctestv,N,percl,percs,conf.level,stratvers=1,
#'               strat="SEX") # percentage stratification
#' out <- PrInDT(data,"real",ctestv,N,percl,percs,conf.level,stratvers=50,
#'               strat="SEX") # stratification with minimum no. of tokens
#' out <- PrInDT(data,"real",ctestv,N,percl,percs,conf.level,thres=0.4) # threshold = 0.4
#'
#' @importFrom stats relevel median predict
#' @importFrom party ctree ctree_control
#' @importFrom splitstackshape stratified
#' @importFrom graphics hist abline par
#'
PrInDT <- function(datain,classname,ctestv=NA,N,percl=NA,percs=1,conf.level=0.95,thres=0.5,stratvers=0,strat=NA,seedl=TRUE,
                               minsplit=NA,minbucket=NA,repvar=NA,indrep=0,valdat=datain,thr=0.5){
  ## input check
  if (typeof(datain) != "list" || typeof(classname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || N <= 0 ||
      !( (0 < percl & percl <= 1) | typeof(percl)=="logical") || !(0 < conf.level & conf.level <= 1) || !(0 <= thres & thres < 1) ||
      !(0 <= stratvers) || !(typeof(strat) %in% c("logical", "character")) || typeof(seedl) != "logical" || !(typeof(minsplit) %in% c("logical","double")) || 
      !(typeof(minbucket) %in% c("logical", "double")) || !(typeof(repvar) %in% c("logical","integer","character")) || !(typeof(indrep) %in% c("integer","double")) ||
      typeof(valdat) != "list" || typeof(thr) != "double") {
    stop("irregular input")
  }
  if (!(all(names(datain) %in% names(valdat)))){
    stop("validation data variables unequal input data variables")
  }
  if (indrep == 1 & all(is.na(repvar)) == TRUE){
    message("\n","Note: repeated measurement analysis disabled: no repeated measurement variable specified","\n")
    indrep <- 0
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
  if (N <= 2){
      stop("Number of repetitions too low (must be > 2)")
    }
  if (seedl){
    set.seed(7654321)  # set seed of random numbers
  }
  ###
  #print(.Random.seed[1:5])  # possible check of random seed
  data <- datain
  repv <- repvar
  names(data)[names(data)==classname] <- "class"
  n_class1 <- table(data$class)[1] # no. of elements of larger class 1
  n_class2 <- table(data$class)[2] # no. of elements of smaller class 2
  if (min(n_class1,n_class2) == 0){
   stop("irregular input: only 1 class")
  }
  if (n_class1 < n_class2){
    # relevel of classes if smaller class first
    data$class <- stats::relevel(data$class, levels(data$class)[2]) # larger class now first
    n_class1 <- table(data$class)[1] # no. of elements of larger class 1
    n_class2 <- table(data$class)[2] # no. of elements of smaller class 2
  }
  if (is.na(percl) == TRUE){
    percl <- percs*n_class2/n_class1
  }
  n <- n_class1 + n_class2 # overall no. of observations in data
  ## reordering of classes: smaller class first in observations
#print(c(n_class1,n_class2))
 # if (n_class1 > n_class2){
   order_class <- order(as.numeric(data$class),decreasing=TRUE) # re-ordered line numbers
 # } else {
 #   order_class <- order(as.numeric(data$class))
 # }
  if (percl*n_class1 < 3 || percs*n_class2 < 3){
      stop("Number of observations in one of the class levels too low (< 3)")
  }
  data <- data[order_class,] # data now reordered: smaller class first
  ## validation set re-ordered
  if (identical(valdat, datain)){
    valdat <- data
    if (indrep == 1) {repv <- repv[order_class]}
  }
  names(valdat)[names(valdat)==classname] <- "class"
  valdat$class <- stats::relevel(valdat$class, levels(data$class)[1]) # larger class now first  ## NEWNEW
  if (!(identical(levels(data$class),levels(valdat$class)))){
    stop("levels of input class variable unequal levels of validation class variable")
  }
  n_class1v <- table(valdat$class)[1] # no. of elements of larger class 1
  n_class2v <- table(valdat$class)[2] # no. of elements of smaller class 2
  if (min(n_class1v,n_class2v) == 0){
   stop("irregular input: only 1 class in valdation set")
  }
  nv <- n_class1v + n_class2v
  if (indrep == 1 & length(repv) != dim(valdat)[1]){
    stop("irregular input: no. of observations different in 'repvar' and validation set")
  }
######
  ## model with all original observations
  ######
  ct <- party::ctree(class ~ ., data = data, control = party::ctree_control(mincriterion=conf.level,minsplit=minsplit,minbucket=minbucket))
  if (thres != 0.5){
    predsprob <- stats::predict(ct,type="prob",newdata=valdat)
    ctpreds <- as.factor(sapply( 1:dim(valdat)[1],
                                 function(x)
                                   ifelse(predsprob[[x]][2] > thres,
                                          levels(valdat$class)[2], levels(valdat$class)[1] )))
  } else {
    ctpreds <- stats::predict(ct,newdata=valdat) # predictions for all validation observations
  }
  if (dim(table(ctpreds))[1] > 1){
    confAll <- table(ctpreds, valdat$class)  # confusion matrix calculation
  # balanced accuracy
    crit2 <- (confAll[1,1] / n_class1v + confAll[2,2] / n_class2v)/2
  } else {
    confAll <- matrix(c(table(valdat$class),0,0),nrow=2,byrow=TRUE)
    colnames(confAll) <- names(table(valdat$class))
    rownames(confAll) <- colnames(confAll)
    crit2 <- 0.5
  }
  crit1All <- "FALSE"
  if (all(is.na(ctestv)) == FALSE) {
    crit1All <- FindSubstr(ct,ctestv) # call of FindSubstr for overall tree
  }
  treeAll <- ct
  baAll <- crit2
acc1AE <- 0
acc2AE <- 0
bamaxAE <- 0
namA1 <-  list()
namA2 <- list()
if (indrep == 1){
  pred <- predict(ct,newdata=valdat)  ## NEWNEW
#print(c(length(repv),length(pred)))
  ch <- table(repv,pred)
  conti <- cbind(table(repv,valdat$class),ch)  ## NEWNEW
#print(conti)
  no1 <- sum(conti[,1] > 0)
  no2 <- sum(conti[,2] > 0)
  if ((no1 + no2) > dim(conti)[1]){
#    cat("Note: repeated measurement variable has more than one class per value")
    acc2AE <- -1
    acc1AE <- -1
  } else {
    for (i in 1:dim(conti)[1]){
      if (conti[i,2] > 0 & conti[i,4] <= (conti[i,2]*thr)){ 
        namA2 <- c(namA2,rownames(ch)[i])
      }
      if (conti[i,1] > 0 & conti[i,3] < (conti[i,1]*thr)){ 
        namA1 <- c(namA1,rownames(ch)[i])
      }
    }
    acc2AE <- 1 - length(namA2) / no2
    acc1AE <- 1 - length(namA1) / no1  
    bamaxAE <- (acc1AE + acc2AE)/2
    if (length(namA1) == 0) {namA1 <- "-"}
    if (length(namA2) == 0) {namA2 <- "-"} 
  }
}
################################################
  ## resampling loop
  ######
  ## initialization
  chs <- matrix(0,nrow=N,ncol=8) # matrix of characteristics
  s_ind <- matrix(0,nrow=N,ncol=n) # matrix for resampled indices
  preds <- matrix(0,nrow=N,ncol=n) # matrix for predictions
  preds_ordered <- preds   # matrix for ordered predictions
  preds_orderedv <- matrix(0,nrow=N,ncol=nv)  # matrix for validation predictions
  crit1 <- rep("FALSE",N) # vector for crit1 (interpretability)
  crit2 <- rep(0,N)      # vector for crit2 (predictive power on validation sample)
  crit3 <- rep(0,N)      # vector for crit3 (predictive power on test sample)
  trees <- ct              # initialization of vector of trees (to be concatenated)
  if (is.na(strat) == FALSE){
    df1 <- data.frame( cbind( 1:n_class1,data[data$class == levels(data$class)[1],names(data)==strat] ) )
    colnames(df1) <- c("x","strat")
    df2 <- data.frame( cbind( 1:n_class2,data[data$class == levels(data$class)[2],names(data)==strat] ) )
    colnames(df2) <- c("x","strat")
    stratv <- as.factor(as.numeric(data[,strat]))
    nmin <- min(sum(df1[,2] == levels(stratv)[1]) + sum(df2[,2] == levels(stratv)[1]),
                sum(df1[,2] == levels(stratv)[2]) + sum(df2[,2] == levels(stratv)[2]))
    if (nmin < stratvers){
      stop("Minimum number of tokens of predictor value for stratification too big")
    }
  }
  ## column names of confusion matrix elements
  colnames(chs) <- c("train11","train21","train12","train22","test11","test21","test12","test22")
  if (is.na(strat) == FALSE & stratvers == 1){
#    message("\n")
    message("Percentage stratification")
    message("\n")
  }
  if (is.na(strat) == FALSE & stratvers > 1){
#    message("\n")
    message("Stratification with minimum number of tokens per predictor value = ",stratvers)
    message("\n")
  }
  ## start of resampling loop
  for (i in 1:N) {  # N repetitions
    # generation of sample indices for undersampling
    if (is.na(strat) == TRUE){
      sample_ind <- sample.int(n_class1,round(percl*n_class1)) + n_class2
      if (percs < 1){
        sample_ind <- c(sample.int(n_class2,round(percs*n_class2)),sample_ind)
      } else {
        sample_ind <- c(seq(1:n_class2),sample_ind)
      }
    }
    if (is.na(strat) == FALSE){
      if (stratvers == 1){
        sample_ind <- t(splitstackshape::stratified( df1,"strat",percl )[,1] + n_class2)
        if (percs < 1){
          suppressMessages(sample_ind <- c(t( splitstackshape::stratified( df2,"strat",percs )[,1] ), sample_ind))
        } else {
          suppressMessages(sample_ind <- c( seq(1:n_class2), sample_ind ))
        }
      } else {
        iter <- 0
        nmin <- 0
        while (nmin < stratvers & iter < 10){
          iter <- iter + 1
          sample_ind <- sample.int(n_class1,round(percl*n_class1)) + n_class2
          if (percs < 1){
            sample_ind <- c(sample.int(n_class2,round(percs*n_class2)),sample_ind)
          } else {
            sample_ind <- c(seq(1:n_class2),sample_ind)
          }
          nmin <- min(sum(data[sample_ind,strat] == levels(data[,strat])[1]),
                      sum(data[sample_ind,strat] == levels(data[,strat])[2]))
        }
        if (iter == 10){
          stop("Minimum number of tokens of predictor value specified for stratification not realizable")
        }
      }
    }
    ls <- length(sample_ind)
    s_ind[i,1:ls] <- sample_ind # chosen indices of larger class
    s_ind[i,(ls + 1):n] <- (seq(1:n))[-sample_ind]
    # generation of train and test samples
    data_train <- data[sample_ind,] # undersamples of both classes
    data_test <- data[-sample_ind,] # rest of sample
    ## modeling
    ct <- party::ctree(class ~ ., data = data_train,control = party::ctree_control(mincriterion=conf.level,minsplit=minsplit,minbucket=minbucket)) 
    trees <- c(trees,ct) # storing of new tree
    if (all(is.na(ctestv)) == FALSE) {
      crit1[i] <- FindSubstr(ct,ctestv) # deciding on interpretability
    }
    ## classification accuracy
    ## training data
    if (thres != 0.5){
      predsprob <- stats::predict(ct,type="prob")
      ctpreds <- as.factor(sapply( 1:dim(data_train)[1], function(x) ifelse( predsprob[[x]][2] > thres,
                                                                             levels(data$class)[2],levels(data$class)[1] )))
    } else {
      ctpreds <- stats::predict(ct)
    }
    if (dim(table(ctpreds))[1] > 1){
      conf_train <- table(ctpreds, data_train$class)
    }
    else {
      conf_train <- matrix(c(table(data_train$class),0,0),nrow=2,byrow=TRUE)
      colnames(conf_train) <- names(table(data$class))
      rownames(conf_train) <- colnames(conf_train)
    }
    preds[i,1:ls] <- ctpreds # training predictions
    chs[i,1:4] <- as.vector(conf_train) # confusion matrix from training
    ## test data 
    if (thres != 0.5){
      predsprob <- stats::predict(ct,newdata=data_test,type="prob")
      ctpreds <- factor(sapply( 1:dim(data_test)[1], function(x) ifelse( predsprob[[x]][2] > thres,
                                                                         levels(data$class)[2],levels(data$class)[1] )))
    } else {
      ctpreds <- stats::predict(ct,newdata=data_test)
    }
    if (dim(table(ctpreds))[1] > 1){
      conf_test <- table(ctpreds, data_test$class)
    }
    else {
      conf_test <- matrix(c(table(data_test$class),0,0),nrow=2,byrow=TRUE)
      colnames(conf_test) <- names(table(data$class))
      rownames(conf_test) <- colnames(conf_test)
    }
    preds[i,(ls+1):n] <- ctpreds # test predictions
    preds_ordered[i,] <- preds[i,order(s_ind[i,])] # correct order of predictions
    chs[i,5:8] <- as.vector(conf_test) # confusion matrix from testing
#   crit2[i] <- 0.5 * ( (chs[i,1] + chs[i,5]) / n_class1 + (chs[i,4] + chs[i,8]) / n_class2 ) # overall balanced accuracy
#if (i > 1 & crit2[i] > max(crit2[1:(i-1)])){
#print(c(i,crit2[i]))
#}
    if ( percs < 1 ){
      crit3[i] <- 0.5 * ( chs[i,5] / (chs[i,5] + chs[i,6]) + chs[i,8] / (chs[i,7] + chs[i,8]) ) # test balanced accuracy
    } else {
      crit3[i] <- chs[i,5] / (chs[i,5] + chs[i,6])  # test balanced accuracy
    }
    ###########################
    ## validation data 
    if (thres != 0.5){
      predsprob <- stats::predict(ct,newdata=valdat,type="prob")
      ctpreds <- factor(sapply( 1:dim(valdat)[1], function(x) ifelse( predsprob[[x]][2] > thres,
                                                                         levels(valdat$class)[2],levels(valdat$class)[1] )))
    } else {
      ctpreds <- stats::predict(ct,newdata=valdat)
    }
    if (dim(table(ctpreds))[1] > 1){
#print(length(ctpreds))
#print(str(valdat))
      conf_test <- table(ctpreds, valdat$class)
    }
    else {
      conf_test <- matrix(c(table(valdat$class),0,0),nrow=2,byrow=TRUE)
      colnames(conf_test) <- names(table(valdat$class))
      rownames(conf_test) <- colnames(conf_test)
    }
    ba1 <- conf_test[1,1] / (conf_test[1,1] + conf_test[2,1])
    ba2 <- conf_test[2,2] / (conf_test[1,2] + conf_test[2,2])
    crit2[i] <- (ba1 + ba2)/2
    preds_orderedv[i,] <- ctpreds ## for validation set
#    preds[i,(ls+1):n] <- ctpreds # test predictions
#    preds_ordered[i,] <- preds[i,order(s_ind[i,])] # correct order of predictions
#    chs[i,5:8] <- as.vector(conf_test) # confusion matrix from testing
#    crit2[i] <- 0.5 * ( (chs[i,1] + chs[i,5]) / n_class1 + (chs[i,4] + chs[i,8]) / n_class2 ) # overall balanced accuracy
  }
  #### end of resampling loop
  ####
  ## analysis of resampled trees
  ####
  trees <- trees[-1]  # analyze without the tree from all data
  ## 3 best balanced accurary trees on validation sample
  crit2o <- sort(crit2,decreasing=TRUE,index.return=TRUE) # sorting the balanced accuracies
  three <- crit2o$ix[crit1[crit2o$ix] == "FALSE"][1:3]    # indices of the 3 best interpretable trees: overall
  crit3o <- sort(crit3,decreasing=TRUE,index.return=TRUE) # sorting the balanced accuracies
  threet <- crit3o$ix[crit1[crit3o$ix] == "FALSE"][1:3]    # indices of the 3 best intrepretable trees: test
  ba1st <- c(crit2[three[1]],crit3[threet[1]])
  ba2nd <- c(crit2o$x[crit1[crit2o$ix] == "FALSE"][2],crit3[threet[2]])
  ba3rd <- c(crit2o$x[crit1[crit2o$ix] == "FALSE"][3],crit3[threet[3]])
  ## accuracies on validation sample
  ## best tree
  classE <- as.integer(valdat$class)
  preds_orderedE <- preds_orderedv[three[1],]
  preds1 <- preds_orderedE[classE == 1]
  preds2 <- preds_orderedE[classE == 2]
  accE1 <- sum (preds1 == 1) / n_class1
  accE2 <- sum (preds2 == 2) / n_class2
  ## accuracies on test sample
  classE <- as.integer(data$class[s_ind[threet[1],(ls + 1):n]])
  preds_orderedE <- preds[threet[1],(ls+1):n]
  preds1 <- preds_orderedE[classE == 1]
  preds2 <- preds_orderedE[classE == 2]
  accT1 <- sum (preds1 == 1) / (chs[threet[1],5] + chs[threet[1],6])
  if ( percs < 1 ){
    accT2 <- sum (preds2 == 2) / (chs[threet[1],7] + chs[threet[1],8])
  } else {
    accT2 <- NA
  }
  ba1st <- c(accE1,accE2,ba1st[1],accT1,accT2,ba1st[2])
  names(ba1st) <- c(names(ba1st)[1:2],"balanced",paste0("test:",names(ba1st)[1]),paste0("test:",names(ba1st)[2]),"balanced")
  ## 2nd best tree
  classE <- as.integer(valdat$class)
  preds_orderedE <- preds_orderedv[three[2],]
  preds1 <- preds_orderedE[classE == 1]
  preds2 <- preds_orderedE[classE == 2]
  accE1 <- sum (preds1 == 1) / n_class1
  accE2 <- sum (preds2 == 2) / n_class2
  ##
  classE <- as.integer(data$class[s_ind[threet[2],(ls + 1):n]])
  preds_orderedE <- preds[threet[2],(ls+1):n]
  preds1 <- preds_orderedE[classE == 1]
  preds2 <- preds_orderedE[classE == 2]
  accT1 <- sum (preds1 == 1) / (chs[threet[2],5] + chs[threet[2],6])
  if ( percs < 1 ){
    accT2 <- sum (preds2 == 2) / (chs[threet[2],7] + chs[threet[2],8])
  } else {
    accT2 <- NA
  }
  ba2nd <- c(accE1,accE2,ba2nd[1],accT1,accT2,ba2nd[2])
  names(ba2nd) <- names(ba1st)
  ## 3rd best tree
  classE <- as.integer(valdat$class)
  preds_orderedE <- preds_orderedv[three[3],]
  preds1 <- preds_orderedE[classE == 1]
  preds2 <- preds_orderedE[classE == 2]
  accE1 <- sum (preds1 == 1) / n_class1
  accE2 <- sum (preds2 == 2) / n_class2
  ##
  classE <- as.integer(data$class[s_ind[threet[3],(ls + 1):n]])
  preds_orderedE <- preds[threet[3],(ls+1):n]
  preds1 <- preds_orderedE[classE == 1]
  preds2 <- preds_orderedE[classE == 2]
  accT1 <- sum (preds1 == 1) / (chs[threet[3],5] + chs[threet[3],6])
  if ( percs < 1 ){
    accT2 <- sum (preds2 == 2) / (chs[threet[3],7] + chs[threet[3],8])
  } else {
    accT2 <- NA
  }
  ba3rd <- c(accE1,accE2,ba3rd[1],accT1,accT2,ba3rd[2])
  names(ba3rd) <- names(ba1st)
  ####
  ## ensembles: accuracies on validation and test sample
  ####
  baen <- matrix(0,nrow=3,ncol=5)
  ## all acceptable trees
  ## Mode function definition
  modus <- function(x){which.max(tabulate(x))} # definition of Mode function
  ##
  classE <- as.integer(valdat$class)
  no1 <- sum(classE == 1)
  no2 <- sum(classE == 2)
  ##
  treesE <- trees[crit1 == 'FALSE'] # trees fulfilling criterion 1
  preds_orderedEv <- preds_orderedv[crit1 == 'FALSE',]
  preds1 <- preds_orderedEv[,classE == 1]
  preds2 <- preds_orderedEv[,classE == 2]
  accE1 <- sum (apply(preds1,2,modus) == 1) / no1
  accE2 <- sum (apply(preds2,2,modus) == 2) / no2
  accE <- (accE1 + accE2) / 2
  accET <- mean(crit3[crit1 == 'FALSE'])
  baen[1,] <- c(length(treesE),accE1,accE2,accE,accET)  # save for print
  ## 3 best trees
 # classE <- as.integer(data$class)
  treesE <- trees[three]
  preds_orderedEv <- preds_orderedv[three,]
  preds1 <- preds_orderedEv[,classE == 1]
  preds2 <- preds_orderedEv[,classE == 2]
  accE1 <- sum (apply(preds1,2,modus) == 1) / no1
  accE2 <- sum (apply(preds2,2,modus) == 2) / no2
  accE <- (accE1 + accE2) / 2
  accET <- mean(crit3[threet])
  baen[2,] <- c(length(treesE),accE1,accE2,accE,accET)  # save for print
  ## trees with overall accuracies better than median
  cthres <- stats::median(crit2) # trees with overall balanced accuracy greater or equal median
  cthrest <- stats::median(crit3) # trees with test balanced accuracy greater of equal median
 # classE <- as.integer(data$class)
  treesE <- trees[crit1 == 'FALSE' & crit2 >= cthres]
  treesET <- trees[crit1 == 'FALSE' & crit3 >= cthrest]
  ##treesEU <- unique(treesE) # (restiction to) unique trees
  ##length(treesEU) # no. of unique trees
  preds_orderedEv <- preds_orderedv[crit1 == 'FALSE' & crit2 >= cthres,]
  preds1 <- preds_orderedEv[,classE == 1]
  preds2 <- preds_orderedEv[,classE == 2]
  accE1 <- sum (apply(preds1,2,modus) == 1) / no1
  accE2 <- sum (apply(preds2,2,modus) == 2) / no2
  accE <- (accE1 + accE2) / 2
  accET <- mean(crit3[crit1 == 'FALSE' & crit3 >= cthrest])
  baen[3,] <- c(length(treesE),accE1,accE2,accE,accET)  # save for print
lablarge <- names(table(data$class))[1] # class 1 = large
labsmall <- names(table(data$class))[2] # class 2 = small 
  ####
acc1E <- 0
acc2E <- 0
bamaxE <- 0
nam1 <-  list()
nam2 <- list()
###
if (indrep == 1){ 
  pred <- predict(trees[[three[1]]],newdata=valdat)  ## NEWNEW
#print(table(pred))
  ch <- table(repv,pred)
#print(dim(ch))
  conti <- cbind(table(repv,valdat$class),ch)  ## NEWNEW
#print(conti)
  no1 <- sum(conti[,1] > 0)
  no2 <- sum(conti[,2] > 0)
  if ((no1 + no2) > dim(conti)[1]){
#    cat("Note: repeated measurement variable has more than one class per value")
    acc2E <- -1
    acc1E <- -1 
  } else {
    for (i in 1:dim(conti)[1]){
      if (conti[i,2] > 0 & conti[i,4] <= (conti[i,2]*thr)){ 
        nam2 <- c(nam2,rownames(ch)[i])
      }
      if (conti[i,1] > 0 & conti[i,3] < (conti[i,1]*thr)){ 
        nam1 <- c(nam1,rownames(ch)[i])
      }
    }
    acc2E <- 1 - length(nam2) / no2
    acc1E <- 1 - length(nam1) / no1  
    bamaxE <- (acc1E + acc2E)/2
    if (length(nam1) == 0) {nam1 <- "-"}
    if (length(nam2) == 0) {nam2 <- "-"}
  }
}
  result <- list(tree1st = trees[[three[1]]], treet1st = trees[[threet[1]]], ba1st = ba1st, tree2nd = trees[[three[2]]], treet2nd = trees[[threet[2]]],
                 ba2nd = ba2nd, tree3rd = trees[[three[3]]], treet3rd = trees[[threet[3]]], ba3rd = ba3rd, baen = baen, bafull= crit2, batest=crit3, 
                 treeAll=treeAll, baAll=baAll, interpAll=crit1All, confAll = confAll,lablarge=lablarge,labsmall=labsmall,
                 acc1AE=acc1AE,acc2AE=acc2AE,bamaxAE=bamaxAE,namA1=namA1,namA2=namA2,
                 acc1E=acc1E,acc2E=acc2E,bamaxE=bamaxE,nam1=nam1,nam2=nam2,indrep=indrep,valdat=valdat,thr=thr)
  class(result) <- c("PrInDT", "PrInDTAll")
#  gc(full=TRUE)
  result
}
#' @export
print.PrInDTAll <- function(x, ...){
  # output for model on all observations
  cat("******************************","\n")
  cat("Tree built on all observations","\n")
  cat("******************************","\n")
  print(x$treeAll) # print of tree
  cat("\n")
  cat("****************","\n")
  cat("Confusion matrix","\n")
  cat("****************","\n")
  print(x$confAll)
  cat("***********************************************","\n")
  cat("Criteria: interpretable    balanced accuracy","\n")
  cat("***********************************************","\n")
  cat("             ",unname(!as.logical(x$interpAll)),"              ",unname(x$baAll),"\n")
if (x$indrep == 1){
  cat("\n","Repeated measurements: Threshold = ",x$thr)
  if (x$acc1AE == -1){
    cat("\n","no analysis: repeated measurement variable has more than one class per value","\n")
  } else {  
    accs <- c(round(x$acc1AE,digits=6),round(x$acc2AE,digits=6),round(x$bamaxAE,digits=6))
    names(accs) <- c(colnames(x$confAll)[1:2],"balanced")
    cat("\n","Full sample accuracies of best tree (Elements)","\n")
    print(accs,row.names=FALSE)
 # cat("\n","Full sample accuracies of best tree (Elements)","\n"," ",colnames(x$confAll)[1]," \t",colnames(x$confAll)[2],"\t","balanced","\n")
 # cat(" ",c(round(x$acc1AE,digits=6),"\t",round(x$acc2AE,digits=6),"\t",round(x$bamaxAE,digits=6)),"\n\n")
    cat("\n","Wrongly predicted Elements:","\n")
    cat(x$lablarge,"\n")
    cat(unlist(x$namA1),"\n")
    cat(x$labsmall,"\n")
    cat(unlist(x$namA2),"\n")
  }
}
  invisible(x)
}
#' @export
print.PrInDT <- function(x, ...){
  NextMethod()
  cat("\n")
  cat("********************************************","\n")
  cat("Best acceptable tree from undersampling according to full sample accuracies","\n")
  cat("********************************************","\n")
  ## best tree: accuracies on full sample
  print(x$tree1st) # structure of tree with best overall accuracy
  accs <- c(round(x$ba1st[1],digits=6),round(x$ba1st[2],digits=6),round(x$ba1st[3],digits=6))
  names(accs) <- c(names(x$ba1st)[1:2],"balanced")
  cat("\n","Validation sample accuracies","\n")
  print(accs,row.names=FALSE)
#  cat("Validation sample accuracies: ",names(x$ba1st)[1],"\t","  ",names(x$ba1st)[2],"\t","  balanced","\n")
#  cat("                              ",unname(x$ba1st)[1],"\t","  ",unname(x$ba1st)[2],"\t","  ",unname(x$ba1st)[3],"\n") 
  cat("\n")
if (x$indrep == 1){
  cat("\n","Repeated measurements: Threshold = ",x$thr)
  if (x$acc1AE == -1){
    cat("\n","no analysis: repeated measurement variable has more than one class per value","\n")
  } else {  
    accs <- c(round(x$acc1E,digits=6),round(x$acc2E,digits=6),round(x$bamaxE,digits=6))
    names(accs) <- c(colnames(x$confAll)[1:2],"balanced")
    cat("\n","Input sample accuracies of best tree (Elements)","\n")
    print(accs,row.names=FALSE)
    cat("\n","Wrongly predicted Elements:","\n")
    cat(x$lablarge,"\n")
    cat(unlist(x$nam1),"\n")
    cat(x$labsmall,"\n")
    cat(unlist(x$nam2),"\n")
  }
}
  cat("\n","*************************************","\n")
  cat("Best acceptable tree from undersampling according to test accuracies","\n")
  cat("*************************************","\n")
  print(x$treet1st) # structure of tree with best test accuracy
  accs <- c(round(x$ba1st[4],digits=6),round(x$ba1st[5],digits=6),round(x$ba1st[6],digits=6))
  names(accs) <- c(names(x$ba1st)[4:5],"balanced")
  cat("\n","Test sample accuracies","\n")
  print(accs,row.names=FALSE)
#  cat("Test sample accuracies: ",names(x$ba1st)[4],"  ",names(x$ba1st)[5],"   balanced","\n")
#  cat("                        ",unname(x$ba1st)[4],"      ",unname(x$ba1st)[5],"         ",unname(x$ba1st)[6],"\n") # balanced accuracy of best tree printed
  cat("\n")
  cat("\n")
  cat("************************************************************************************************","\n")
  cat("Ensemble of all acceptable trees","\n")
  accs <- c(x$baen[1,1:5])
  names(accs) <- c("no. of trees",names(x$ba1st)[1:2],"balanced","  mean test balanced")
  print(accs,row.names=FALSE)
  #,names(x$ba1st)[1],"  ",names(x$ba1st)[2],"    balanced (mean test balanced)","\n"
  #cat("     ",unname(x$baen)[1,1],"                              ",unname(x$baen)[1,2],unname(x$baen)[1,3],unname(x$baen)[1,4],"   (",unname(x$baen)[1,5],")\n") # print no. of trees and accuracies
  cat("\n")
  cat("************************************************************************************************","\n")
  cat("Ensemble of 3 best acceptable trees","\n")
  accs <- c(x$baen[2,1:5])
  names(accs) <- c("no. of trees",names(x$ba1st)[1:2],"balanced","  mean test balanced")
  print(accs,row.names=FALSE)
 # cat("************************************************************************************************","\n")
 # cat("Ensemble of acceptable 3 best trees on full or test sample","\n"," no. of trees; full sample accuracies:",names(x$ba1st)[1],"  ",names(x$ba1st)[2],"   balanced (mean test balanced)","\n")
 # cat("******************************************************************************************","\n")
 # cat("     ",unname(x$baen)[2,1],"                               ",unname(x$baen)[2,2],unname(x$baen)[2,3],unname(x$baen)[2,4],"  (",unname(x$baen)[2,5],")\n") # print no. of trees and accuracies
  cat("\n")
  cat("************************************************************************************************","\n")
  cat("Ensemble of acceptable trees with full or test sample accuracies >= median ","\n")
  accs <- c(x$baen[3,1:5])
  names(accs) <- c("no. of trees",names(x$ba1st)[1:2],"balanced","  mean test balanced")
  print(accs,row.names=FALSE)
#  cat("*******************************************************************************************************","\n")
#  cat("Ensemble of acceptable trees with full or test sample accuracies >= median","\n"," no. of trees: full (test); full sample accuracies:",names(x$ba1st)[1],"  ",names(x$ba1st)[2],"   balanced (mean test #balanced)","\n")
#  cat("*******************************************************************************************************","\n")
#  cat("     ",unname(x$baen)[3,1],"                                          ",unname(x$baen)[3,2],unname(x$baen)[3,3],unname(x$baen)[3,4],"  (",unname(x$baen)[3,5],")\n") # print no. of trees and accuracies
  invisible(x)
}
#' @export
plot.PrInDTAll <- function(x, ...){
  plot(x$treeAll,main="Tree based on all observations")
  invisible(x)
}
#' @export
plot.PrInDT <- function(x, ...){
#  NextMethod()
  ## plot histograms of balanced accuracies
#  par(mfrow=c(1,2))
  hist(x$bafull,main=" ",xlab="Full sample balanced accuracies",cex.axis=1.5,cex.lab=1.5)
  abline(v = stats::median(x$bafull),col='red',lwd = 3)
  hist(x$batest,main=" ",xlab="Test sample balanced accuracies",cex.axis=1.5,cex.lab=1.5)
  abline(v = stats::median(x$batest),col='red',lwd = 3)
## Interpretable tree with highest full sample accuracy
  plot(x$tree1st,main="Interpretable tree from undersampling with highest full sample accuracy")
## Interpretable tree with highest test sample accuracy
  plot(x$treet1st,main="Interpretable tree from undersampling with highest test sample accuracy") 
  invisible(x)
}
#' @export
predict.PrInDT <- function(object,...){
  predict(object$tree1st,...)
}
#' @export
predict.PrInDTAll <- function(object,...){
  predict(object$treeAll,...)
}

