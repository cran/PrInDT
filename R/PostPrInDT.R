#' Posterior analysis of conditional inference trees: distribution of a specified variable in the terminal nodes.
#'
#' @description The conditional inference tree 'ct' is analyzed according to the distribution of a variable 'var' in its terminal nodes.\cr
#' In the case of a discrete variable 'var', the appearance of the different levels is considered for each terminal node. \cr
#' In the case of a continuous variable 'var', means and standard deviations of 'var' or the target variable are considered for each terminal node.\cr
#' In particular, this function can be used for the posterior analysis of a tree regarding the distribution of a variable not present in the tree.
#' 
#' @usage PostPrInDT(datain, ct, target, var, vardata, vt)
#'
#' @param datain input data frame with the observatios of all variables used in the model 
#' @param ct conditional inference tree to be analyzed
#' @param target name of target variable of 'ct' (character)
#' @param var name of variable of interest (character)
#' @param vardata observations of 'var'
#' @param vt type of variables: 'dd' for discrete target (classification) and discrete variable 'var', 'dc' for discrete target (classification) and continuous 'var',
#' 'cd' for continuous target (regression) and discrete 'var', and 'cc' for continuous target (regression) and continuous 'var'.
#'
#' @return None: Relevant output is produced by the function.
#'
#' @export PostPrInDT 
#'
#' @examples
#' data <- PrInDT::data_zero
#' data <- na.omit(data)
#' outAll <- PrInDTAll(data,"real") 
#' PostPrInDT(data,outAll$treeAll,"real","ETH",data$ETH,vt="dd")
#' PostPrInDT(data,outAll$treeAll,"real","AGE",data$AGE,vt="dc")
#' datareg <- PrInDT::data_vowel
#' outregAll <- PrInDTregAll(datareg,"target") 
#' PostPrInDT(datareg,outregAll$treeAll,"target","Nickname",datareg$Nickname,vt="cd")
#' PostPrInDT(datareg,outregAll$treeAll,"target","AGE",datareg$AGE,vt="cc")
#'
#' @import party
#'
PostPrInDT <- function(datain,ct,target,var,vardata,vt){
  ## input check
  if (typeof(datain) != "list" || typeof(ct) != "S4" || typeof(var) != "character" || typeof(target) != "character" ||
      !(typeof(vardata) %in% c("integer","double")) || typeof(vt) != "character"){
    stop("irregular input")
  }
####
  data <- datain
  colnames(data)[colnames(data)==target] <- "target"
  if (vt == "dd" | vt == "cd"){
    check <- cbind( as.character(predict(ct,newdata=data)) , where(ct,newdata = data), as.character(vardata))
    tvar <- table(check[,2],check[,3])
    cat("\n","Occurence of levels of variable *** ",var," *** in terminal nodes of the following tree")
    print(ct)
    cat("\n","Nodes in rows, levels in columns")
    print(tvar)
  }
  if (vt == "dc" | vt == "cc"){
    cat("\n","Means and standard deviations of values of variable *** ",var," *** in terminal nodes of the following tree","\n")
    print(ct)
    tt <- table(where(ct,newdata=data))
    means <- rep( 0, dim(tt)[1] )
    stdevs <- means
    for (i in 1:dim(tt)[1]) {
        mm <- mean( vardata[ where(ct,newdata=data) == names(tt)[i] ])
        if (is.na(mm) == FALSE) {means[i] <- mm}
        st <- sd( vardata[ where(ct,newdata=data) == names(tt)[i] ] )
        if (is.na(st) == FALSE) {stdevs[i] <- st}
    }
    tvar <- rbind(means,stdevs)
    rownames(tvar) <- c("means","stdevs")
    colnames(tvar) <- names(tt)
    cat("\n","Means and standard deviations in nodes","\n")
    print(tvar)
  }
  if (vt == "cd" |vt == "cc"){
    cat("\n","Means and standard deviations of values of target variable *** ",target," *** in terminal nodes of the following tree","\n")
    print(ct)
    tt <- table(where(ct,newdata=data))
    means <- rep( 0, dim(tt)[1] )
    stdevs <- means
    for (i in 1:dim(tt)[1]) {
        mm <- mean( data$target[ where(ct,newdata=data) == names(tt)[i] ])
        if (is.na(mm) == FALSE) {means[i] <- mm}
        st <- sd( data$target[ where(ct,newdata=data) == names(tt)[i] ] )
        if (is.na(st) == FALSE) {stdevs[i] <- st}
    }
    tvar <- rbind(means,stdevs)
    rownames(tvar) <- c("means","stdevs")
    colnames(tvar) <- names(tt)
    cat("\n","Means and standard deviations in nodes","\n")
    print(tvar)
  }
}