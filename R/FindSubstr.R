#' Check for forbidden split results in trees
#'
#' Check whether one of the character strings in the vector 'ctestv' appears as a split result in the conditional inference tree 'ct';
#' ctestv is a vector of character strings of forbidden split results.\cr
#'     Example: ctestv <- rbind('variable1 == \{value1, value2\}','variable2 <= value3'), where
#'     character strings specified in 'value1', 'value2' are not allowed as results of a splitting operation in 'variable 1' in a tree.
#'     For restrictions of the type 'variable <= xxx', all split results in a tree are excluded with 'variable <= yyy' and yyy <= xxx.\cr
#'     Trees with split results specified in 'ctestv' are not accepted during optimization.\cr
#'     A concrete example is: 'ctestv <- rbind('ETH == \{C2a, C1a\}','AGE <= 20')' for variables 'ETH' and 'AGE' and values 'C2a','C1a', and '20';\cr
#'     For an application, please refer to, e.g., the functions \code{\link{PrInDT}} and \code{\link{PrInDTreg}}.\cr
#'     If no restrictions exist, the default = NA is used.
#'
#' @usage FindSubstr(ct, ctestv)
#'
#' @param ct Tree to be checked
#' @param ctestv Vector with character strings of excluded split results
#'
#' @return
#' \describe{ 
#' \item{testt}{TRUE if any of the split results in 'ctestv' appears in 'ct'; FALSE otherwise}
#' }
#'
#' @importFrom utils capture.output
#'
#' @noRd
#'
FindSubstr <- function(ct,ctestv){
  ## input check
  if (typeof(ct) != "S4" || typeof(ctestv) != "character"){
    stop("irregular input")
  }
  ##
  J <- length(ctestv)
  ctvec <- utils::capture.output(print(ct))
  K <- length(ctvec)
  testc <- array(rep("FALSE",J*K),c(J,K))
  for(j in 1:J){             # loop over strings to be tested
     nj <- nchar(ctestv[j])
     for(k in 1:K){     # loop over parts of tree description
         nc <- nchar(ctvec[k]) - nj     # nc = length of part of tree - length of tested string
         if (nc >= 0) {
           if (grepl( "<=",ctestv[j],fixed=TRUE ) == FALSE){
             testc[j,k] <- grepl( ctestv[j],ctvec[k],fixed=TRUE )
             }
           else {
            if (grepl( ") ",ctvec[k],fixed=TRUE) == TRUE){
             if ( strsplit(strsplit(ctvec[k],"<=")[[1]][1],") ")[[1]][2] == strsplit(ctestv[j],"<=")[[1]][1] &
               as.numeric(gsub(".*?([0-9]+).*", "\\1",ctvec[k])) <= as.numeric(gsub(".*?([0-9]+).*", "\\1",ctestv[j])) ){
               testc[j,k] <- "TRUE"
             }
            }
           }
         }
     }
  }
  testt <- any(as.logical(testc))
  return(testt)
}
