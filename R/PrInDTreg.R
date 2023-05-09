#' Regression tree resampling by the PrInDT method
#'
#' @description Regression tree optimzation to identify the best interpretable tree; interpretability is checked (see 'ctestv').\cr
#' The relationship between the target variable 'regname' and all other factor and numerical variables
#' in the data frame 'datain' is optimally modeled by means of 'N' repetitions of subsampling.\cr 
#' The optimization criterion is the R2 of the model on the full sample.\cr
#' Multiple subsampling percentages of observations and predictors can be specified (in 'pobs' and 'ppre', correspondingly).\cr
#' The trees generated from undersampling can be restricted by
#' rejecting unacceptable trees which include split results specified in the character strings of the vector 'ctestv'.\cr
#'
#' @usage PrInDTreg(datain, regname, ctestv=NA, N, pobs, ppre, conf.level=0.95)
#'
#' @param datain Input data frame with class factor variable 'classname' and the\cr
#'    influential variables, which need to be factors or numericals (transform logicals and character variables to factors) 
#' @param regname name of regressand variable (character)
#' @param ctestv Vector of character strings of forbidden split results;\cr
#'     {see function \code{\link{PrInDT}} for details.}\cr
#'     If no restrictions exist, the default = NA is used.
#' @param N Number of repetitions (integer > 0)
#' @param pobs Vector of resampling percentages of observations (numerical, > 0 and <= 1)
#' @param ppre Vector of resampling percentages of predictor variables (numerical, > 0 and <= 1)
#' @param conf.level (1 - significance level) in function \code{ctree} (numerical, > 0 and <= 1);\cr
#'     default = 0.95
#'
#' @return 
#' \describe{
#' \item{meanint}{Mean number of interpretable trees over the combinations of individual percentages in 'pobs' and 'ppre'}
#' \item{R2mean}{Mean R2 on test sets}
#' \item{ctmax}{best resampled regression tree according to R2 on the full data set}
#' \item{percmax}{Maximum R2 achieved for \%observations}
#' \item{perfeamax}{Maximum R2 achieved for \%predictors}
#' \item{maxR2}{best R2 on the full data set for resampled regression trees (for 'ctmax')} 
#' \item{interpmax}{interpretability of best tree 'ctmax'}
#' \item{ctmax2}{second best resampled regression tree according to R2 on the full data set}
#' \item{percmax2}{second best R2 achieved for \%observations}
#' \item{perfeamax2}{second best R2 achieved for \%features}
#' \item{max2R2}{second best R2 on the full data set for resampled regression trees (for 'ctmax2')}
#' \item{interp2max}{interpretability of second-best tree 'ctmax2'}
#' }
#'
#' @details
#' For the optimzation of the trees, we employ a method we call Sumping (Subsampling umbrella of 
#' model parameters), a variant of Bumping (Bootstrap umbrella of model parameters) (Tibshirani 
#' & Knight, 1999) which use subsampling instead of bootstrapping. The aim of the 
#' optimization is to identify conditional inference trees with maximum predictive power
#' on the full sample under interpretability restrictions. 
#'
#' \strong{Reference}\cr Tibshirani, R., Knight, K. 1999. Model Search and Inference By Bootstrap "bumping".
#'            Journal of Computational and Graphical Statistics, Vol. 8, No. 4 (Dec., 1999), pp. 671-686 
#'  
#' Standard output can be produced by means of \code{print(name)} or just \code{ name } as well as \code{plot(name)} where 'name' is the output data 
#' frame of the function.\cr
#' The plot function will produce a series of more than one plot. If you use R, you might want to specify \code{windows(record=TRUE)} before 
#' \code{plot(name)} to save the whole series of plots. In R-Studio this functionality is provided automatically.
#'
#' @exportS3Method print PrInDTreg
#' @exportS3Method plot PrInDTreg
#' @export PrInDTreg
#'
#' @examples
#' data <- PrInDT::data_vowel
#' data <- na.omit(data)
#' ctestv <- 'vowel_maximum_pitch <= 320'
#' N <- 30 # no. of repetitions
#' pobs <- c(0.70,0.60) # percentages of observations
#' ppre <- c(0.90,0.70) # percentages of predictors
#' outreg <- PrInDTreg(data,"target",ctestv,N,pobs,ppre)
#' outreg
#' plot(outreg)
#'
#' @importFrom party ctree ctree_control
#' @importFrom stats formula
#'
PrInDTreg <- function(datain,regname,ctestv=NA,N,pobs,ppre,conf.level=0.95){
  ## input check
  if (typeof(datain) != "list" || typeof(regname) != "character" || !(typeof(ctestv) %in% c("logical", "character")) || N <= 0 ||
      !all(0 < pobs & pobs <= 1) || !all(0 < ppre & ppre <= 1) || !(0 < conf.level & conf.level <= 1)){
    stop("irregular input")
  }
####
  data <- datain
  maxR2 <- 0
  max2R2 <- 0
  R2mean <- matrix(0,nrow = length(pobs), ncol = length(ppre))
  interpmax <- FALSE
  interp2max <- FALSE
  meanint <- 0
  set.seed(1234567)
  if (any(pobs*dim(data)[1] < 10)){
      stop("For one of the percentages: Number of tokens too small (< 10)")
  }
  if (any(ppre*dim(data)[2] < 1)){
      stop("For one of the percentages: Number of predictors too small (< 1)")
  }
  message("combinations of percentages of","\n","observations & predictors")
  for (i in 1:length(pobs)) {
    perc <- pobs[i]
    for (j in 1:length(ppre)) {
      perfea <- ppre[j]
      message("    ",perc,"         ",perfea)
      R2 <- 0
      crit1 <- rep(FALSE,N) # vector for crit1 (interpretability)
      for (k in 1:N) {
        samplerows <- sample(1:nrow(data), perc*dim(data)[1],replace=FALSE)
        datasub <- data[samplerows,]
        datatest <- data[-samplerows,]
        datasub2 <- datasub[,sample(1:ncol(data), perfea*dim(data)[2],replace=FALSE)]
        datasub2[,regname] <- datasub[,regname]
        ct <- party::ctree(formula(paste(regname, "~ .")), data = datasub2,control = party::ctree_control(mincriterion=conf.level))
        if (all(is.na(ctestv)) == FALSE) {
          crit1[k] <- FindSubstr(ct,ctestv) # deciding on interpretability
        }
        if (crit1[k] == FALSE){
          meanint <- meanint + 1
        }
        if (maxR2 == 0 & crit1[k] == FALSE) {
          ctmax <- ct
          percmax <- perc
          perfeamax <- perfea
          interpmax <- crit1[k]
        }
        ctpreds <- predict(ct,newdata=data)
        R2 <- 1 - sum((data[,regname] - ctpreds)^2) / sum((data[,regname]-mean(data[,regname]))^2)
        if (R2 > maxR2 & crit1[k] == FALSE) {
          max2R2 <- maxR2
          maxR2 <- R2
          ctmax2 <- ctmax
          ctmax <- ct
          percmax2 <- percmax
          perfeamax2 <- perfeamax
          percmax <- perc
          perfeamax <- perfea
          interp2max <- interpmax
          interpmax <- crit1[k]
        }
        ctpreds <- predict(ct,newdata=datatest)
        R2 <- 1 - sum((datatest[,regname] - ctpreds)^2) / sum((datatest[,regname]-mean(datatest[,regname]))^2)
        R2mean[i,j] <- R2mean[i,j] + R2
      }
      R2mean[i,j] <- R2mean[i,j] / N
    }
  }
  meanint <- meanint / (length(pobs)*length(ppre)) # mean number of interpretable trees over the combinations
  colnames(R2mean) <- as.character(ppre)
  rownames(R2mean) <- as.character(pobs)
  result <- list(meanint = meanint, R2mean = R2mean, ctmax = ctmax, percmax = percmax, perfeamax = perfeamax, maxR2 = maxR2, 
    interpmax = interpmax, ctmax2 = ctmax2, percmax2 = percmax2, perfeamax2 = perfeamax2, max2R2 = max2R2, interp2max = interp2max)
  class(result) <- "PrInDTreg"
  result
}
## print function
print.PrInDTreg <- function(x,...){
  cat("\n","Mean number of interpretable trees over the combinations: ","\n")
  print(x$meanint)
  cat("\n","Mean R2 on test sets over the combinations: ","\n")
  print(x$R2mean)
#
  cat("\n\n","****Best interpretable regression tree****","\n")
  cat("\n","Characteristics","\n")
  cat("  %obs","   %pred","   : ","R2 on full sample","\n") #      interpretable","\n")
  cat("  ",x$percmax,"   ",x$perfeamax,"    :    ",x$maxR2,"\n") #"    ",!as.logical(x$interpmax),"\n")
  print(x$ctmax)
#
  cat("\n\n","****Second best interpretable regression tree****","\n")
  cat("\n","Characteristics","\n")
  cat("  %obs","   %pred","   : ","R2 on full sample","\n") #      interpretable","\n")
  cat("  ",x$percmax2,"   ",x$perfeamax2,"    :    ",x$max2R2,"\n") #"    ",!as.logical(x$interp2max),"\n")
  print(x$ctmax2)
}
## plot function
plot.PrInDTreg <- function(x,...){
  plot(x$ctmax,main="Best interpretable regression tree")
  plot(x$ctmax2,main="Second best interpretable regression tree")
}
