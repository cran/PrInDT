#' @noRd
#'
PrInDTstruc1 <- function(data,class,ctestv=NA,name,check,labs,Ni,N=NA,Eit,Pit,undersamp=TRUE,co,crit="ba",ktest=0,stest=integer(length=0),p1,p2,minsplit=20,minbucket=7,valdat){
#names(data)[names(data)==class] <- "class"
#names(valdat)[names(valdat)==class] <- "class"
#name <- Struc$name
##check <- Struc$check
#check <- eval(parse(text = Struc$check)) 
#labs <- Struc$labs
lMit <- length(Eit)
lDit <- length(Pit)
if (Ni > 1){
  N <- Ni
}
if (is.na(N) == TRUE){
  N <- 99
}
lablarge <- names(table(data$class))[1] # class 1 = large
labsmall <- names(table(data$class))[2] # class 2 = small 
bamax <- 0
tamax <- 0
D <- dim(data)[2] - 1
kumin <-  length(stest)
summin <- length(stest)
mindlong <- list()
interp <- 0
# print(c(Eit,N,Ni))
for (M in Eit){
  message("Number of substructure elements per category: ",M)
  set.seed(7654321)
  for (i in 1:N){
#          ind1 <- sample(unique(name[check %in% labs[1,]]),prob=p1)[1:M] 
#          ind2 <- sample(unique(name[check %in% labs[2,]]),prob=p2)[1:M]
    ind1 <- sample(names(p1),prob=p1)[1:M] 
    ind2 <- sample(names(p2),prob=p2)[1:M]
    gTrain <- data[(name %in% ind1 | name %in% ind2),]
    gTest <- data[!(name %in% ind1 | name %in% ind2),]
    gname <- name[(name %in% ind1 | name %in% ind2)]
    gnameT <- name[!(name %in% ind1 | name %in% ind2)]
    gcheck <- check[(name %in% ind1 | name %in% ind2)]
    if (undersamp == TRUE & crit != "tal"){  ## for "tal" no effect of undersamp! NEWNEW
      gTrainr <- gTrain[gTrain$class==lablarge,]
      gnamer <- gname[gTrain$class==lablarge]
      gcheckr <- gcheck[gTrain$class==lablarge]
      gTrainz <- gTrain[gTrain$class==labsmall,]
      gnamez <- gname[gTrain$class==labsmall]
      gcheckz <- gcheck[gTrain$class==labsmall]
      samp <- sample(1:nrow(gTrainr))
      zTrain <- gTrainr[samp,]
      zname <- gnamer[samp]
      zcheck <- gcheckr[samp]
      gTest <- data[!(name %in% ind1 | name %in% ind2),]
      zTTrain <- rbind(zTrain[1:dim(gTrainz)[1],],gTrainz) # check if dim(gTrainz)[1] > 0 !!!!
      zTname <- c(as.character(zname[1:dim(gTrainz)[1]]),as.character(gnamez))
      nt <- nrow(zTTrain)
      zTcheck <- c(as.character(zcheck[1:dim(gTrainz)[1]]),as.character(gcheckz))
      gTest <- rbind(gTest,zTrain[(dim(gTrainz)[1]+1):(dim(gTrainr)[1]),])
    } else {
        zTTrain <- gTrain
        nt <- nrow(zTTrain)
        zTname <- gname
        zTcheck <- gcheck
    }
    if (Ni == 0){  
      datat <- zTTrain
      outAll <- PrInDTAll(datat,class,ctestv=ctestv,conf.level=co,minsplit=minsplit,minbucket=minbucket)
      if (outAll$interpAll == FALSE){
        interp <- interp + 1
      }
      indmax <- colnames(datat)
      if (crit == "ba"){
        class.pred <- predict(outAll$treeAll,newdata=valdat)
        class.pred <- relevel(class.pred,ref=levels(valdat$class)[1])
      # sum(class.pred == data$class) / length(class.pred)
        conti <- table(class.pred,valdat$class)
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        ba <- (ba1 + ba2)/2
        if (ba > bamax & outAll$interpAll == FALSE){
          bamax <- ba
          acc1 <- ba1
          acc2 <- ba2
          ntmax <- nt
          dmax <- D
          indmax <- colnames(datat)[1:D]
          ind1max <- ind1
          ind2max <- ind2
          outmax <- outAll
          gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
          colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
          gmaxTest <- gTest
        }
      }
      if (crit == "bat"){
        class.pred <- predict(outAll$treeAll,newdata=gTest)
        class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
      # sum(class.pred == gTest$class) / length(class.pred)
        conti <- table(class.pred,gTest$class)
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        ba <- (ba1 + ba2)/2
        if (ba > bamax & outAll$interpAll == FALSE){
          bamax <- ba
          acc1 <- ba1
          acc2 <- ba2
          ntmax <- nt 
          dmax <- D
          indmax <- colnames(datat)[1:D]
          ind1max <- ind1
          ind2max <- ind2
          outmax <- outAll
          gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
          colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
          gmaxTest <- gTest
        }
      }
      if (crit == "ta"){
        class.pred <- predict(outAll$treeAll,newdata=gTest)
        class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
        conti <- table(class.pred,gTest$class)
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ta1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ta2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        ta <- sum(class.pred == gTest$class) / length(class.pred)
        bat <- (ta1 + ta2)/2
        if (ta > tamax & outAll$interpAll == FALSE){
          bamax <- bat
          tamax <- ta
          acc1 <- ta1
          acc2 <- ta2
          ntmax <- nt 
          dmax <- D
          indmax <- colnames(datat)[1:D]
          ind1max <- ind1
          ind2max <- ind2
          outmax <- outAll
          gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
          colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
          gmaxTest <- gTest
        }
      }
      if (crit == "tal"){
        sTest <- stest[!(name %in% ind1 | name %in% ind2)]  ### NEWNEW
        splitTest <- table(gnameT)  
        for (i in 2:length(splitTest)){
          splitTest[i] <- splitTest[i] + splitTest[i-1]  #### splits according to gTest
        }
        class.pred <- predict(outAll$treeAll,newdata=gTest)
        class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
        conti <- table(class.pred,gTest$class)
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ta1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ta2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        bat <- (ta1 + ta2)/2
        ta <- sum(class.pred == gTest$class) / length(class.pred)
        predind <- class.pred == gTest$class
        indfal <- which(predind == FALSE)
        indfal <- indfal[!(indfal %in% splitTest)]
        indsplit <- split(indfal, cumsum(c(1, diff(indfal) != 1)))
        j <- 0
        indlong <- list()
        indstart <- vector()
        indlen <- vector()
        sumlen <- 0
        for (l in 1:length(indsplit)){
          if (length(unlist(indsplit[as.character(l)])) > ktest){ 
            j <- j + 1
            indlong[j] <- indsplit[as.character(l)]
            indstart[j] <- unlist(indlong[j])[1]
            indlen[j] <- length(unlist(indlong[j]))
            sumlen <- sumlen + indlen[j]
          }
        }
        kstart <- length(indstart)
        kustart <- length(unique(sTest[indstart]))
        if ((kustart < kumin & outAll$interpAll == FALSE) |
            (kustart == kumin & ta > tamax & outAll$interpAll == FALSE) |
            (kustart == kumin & ta == tamax & sumlen < summin & outAll$interpAll == FALSE) ){
          kumin <- kustart
          mindstart <- indstart  ## NEWNEW
          mindlong <- indlong
          elems <- unique(sTest[mindstart])
          summin <- sumlen
#print(c(i,kumin))
          tamax <- ta
          bamax <- bat
          acc1 <- ta1
          acc2 <- ta2
          ntmax <- nt 
          dmax <- D
          indmax <- colnames(datat)[1:D]
          ind1max <- ind1
          ind2max <- ind2
          outmax <- outAll
          gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
          colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
          gmaxTest <- gTest
        } 
      } 
    } else {
      for (d in Pit){
        for (n in 1:Ni){
#           datat <- zTTrain[,-1]
          datat <- zTTrain
          datat$class <- NULL
          ind <- sample(1:D)[1:d]
          datat <- datat[,ind]
#            datat <- cbind(datat,zTTrain[,1])
          datat <- cbind(datat,zTTrain$class)
          colnames(datat)[d+1] <- "class"
#print(colnames(datat))
          outAll <- PrInDTAll(datat,class,ctestv=ctestv,conf.level=co,minsplit=minsplit,minbucket=minbucket)
          if (outAll$interpAll == FALSE){
            interp <- interp + 1
          }
          if (crit == "ba"){
            class.pred <- predict(outAll$treeAll,newdata=valdat)
            class.pred <- relevel(class.pred,ref=levels(valdat$class)[1])
      # sum(class.pred == data$class) / length(class.pred)
            conti <- table(class.pred,valdat$class)
            if (dim(conti)[1] == 1){                                                ## NEWNEW
              conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
            }
            ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
            ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
            ba <- (ba1 + ba2)/2
            if (ba > bamax & outAll$interpAll == FALSE){
              bamax <- ba
#print(c(i,d,bamax))
              acc1 <- ba1
              acc2 <- ba2
              ntmax <- nt
              dmax <- d
              indmax <- colnames(datat)[1:d]
              ind1max <- ind1
              ind2max <- ind2
              outmax <- outAll
              gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
              colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
              gmaxTest <- gTest
            }
          }
          if (crit == "bat"){
            class.pred <- predict(outAll$treeAll,newdata=gTest)
            class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
      # sum(class.pred == gTest$class) / length(class.pred)
            conti <- table(class.pred,gTest$class)
            if (dim(conti)[1] == 1){                                                ## NEWNEW
              conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
            }
            ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
            ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
            ba <- (ba1 + ba2)/2
            if (ba > bamax & outAll$interpAll == FALSE){
              bamax <- ba
              acc1 <- ba1
              acc2 <- ba2
              ntmax <- nt 
              dmax <- d
              indmax <- colnames(datat)[1:d]
              ind1max <- ind1
              ind2max <- ind2
              outmax <- outAll
              gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
              colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
              gmaxTest <- gTest
            }
          }
          if (crit == "ta"){
            class.pred <- predict(outAll$treeAll,newdata=gTest)
            class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
            conti <- table(class.pred,gTest$class)
            if (dim(conti)[1] == 1){                                                ## NEWNEW
              conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
            }
            ta1 <- conti[1,1] / (conti[1,1] + conti[2,1])
            ta2 <- conti[2,2] / (conti[1,2] + conti[2,2])
            ta <- sum(class.pred == gTest$class) / length(class.pred)
            bat <- (ta1 + ta2)/2
            if (ta > tamax & outAll$interpAll == FALSE){
              bamax <- bat
              tamax <- ta
              acc1 <- ta1
              acc2 <- ta2
              ntmax <- nt 
              dmax <- d
              indmax <- colnames(datat)[1:d]
              ind1max <- ind1
              ind2max <- ind2
              outmax <- outAll
              gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
              colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
              gmaxTest <- gTest
            }
          }
          if (crit == "tal"){
            sTest <- stest[!(name %in% ind1 | name %in% ind2)]  ### NEWNEW
            splitTest <- table(gnameT)  
            for (i in 2:length(splitTest)){
              splitTest[i] <- splitTest[i] + splitTest[i-1]  #### splits according to gTest
            }
            class.pred <- predict(outAll$treeAll,newdata=gTest)
            class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
            conti <- table(class.pred,gTest$class)
            if (dim(conti)[1] == 1){                                                ## NEWNEW
              conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
            }
            ta1 <- conti[1,1] / (conti[1,1] + conti[2,1])
            ta2 <- conti[2,2] / (conti[1,2] + conti[2,2])
            bat <- (ta1 + ta2)/2
            ta <- sum(class.pred == gTest$class) / length(class.pred)
            predind <- class.pred == gTest$class
            indfal <- which(predind == FALSE)
            indfal <- indfal[!(indfal %in% splitTest)]
            indsplit <- split(indfal, cumsum(c(1, diff(indfal) != 1)))
            j <- 0
            indlong <- list()
            indstart <- vector()
            indlen <- vector()
            sumlen <- 0
            for (l in 1:length(indsplit)){
              if (length(unlist(indsplit[as.character(l)])) > ktest){ 
                j <- j + 1
                indlong[j] <- indsplit[as.character(l)]
                indstart[j] <- unlist(indlong[j])[1]
                indlen[j] <- length(unlist(indlong[j]))
                sumlen <- sumlen + indlen[j]
              }
            }
            kstart <- length(indstart)
            kustart <- length(unique(sTest[indstart]))
            if ((kustart < kumin & outAll$interpAll == FALSE) |
              (kustart == kumin & ta > tamax & outAll$interpAll == FALSE) |
              (kustart == kumin & ta == tamax & sumlen < summin & outAll$interpAll == FALSE) ){
              kumin <- kustart
              mindstart <- indstart  ## NEWNEW
              mindlong <- indlong
#print(c(kumin,mindstart))
#print(unique(sTest[mindstart]))
              elems <- unique(sTest[mindstart])
              summin <- sumlen
              tamax <- ta
              bamax <- bat
              acc1 <- ta1
              acc2 <- ta2
              ntmax <- nt 
              dmax <- d
              indmax <- colnames(datat)[1:d]
              ind1max <- ind1
              ind2max <- ind2
              outmax <- outAll
              gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
              colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
              gmaxTest <- gTest
            }
          }
        }
      }
    } 
  }  
  if (crit == "ba"){
    message("current optimal balanced accuracy: ",round(bamax,digits=4),"\n")
    tamax <- 0
    kumin <- 0
    elems <- 0
  }
  if (crit == "bat"){
    message("current optimal balanced test accuracy: ",round(bamax,digits=4),"\n")
    tamax <- 0
    kumin <- 0
    elems <- 0
  }
  if (crit == "ta"){
    message("current optimal test accuracy: ",round(tamax,digits=4),"\n")
    elems <- 0
    kumin <- 0
  }
  if (crit == "tal"){
    message("current optimal number of long misclassified parts: ",kumin)
    message("current optimal test accuracy: ",round(tamax,digits=4),"\n")
  }
}
if (Ni > 0){
  interp <- c(interp,(lMit*N*Ni*length(Pit)))
} else {
  interp <- c(interp,(lMit*N))
}
##
  result <- list(interp=interp, dmax=dmax, ntmax = ntmax, acc1 = acc1, acc2 = acc2, 
     bamax=bamax,tamax=tamax,kumin=kumin,elems=elems,mindlong=mindlong,
     ind1max = ind1max, ind2max = ind2max, indmax=indmax, modbest = outmax$treeAll, gmaxTrain = 
     gmaxTrain, gmaxTest = gmaxTest)
}