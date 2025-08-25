#' @noRd
#'
PrInDTstruc2 <- function(data,class,ctestv=NA,name,check,labs,N,Ni=NA,Pit,Eit,undersamp=TRUE,co,crit="ba",ktest=0,stest=integer(length=0),p1,p2,minsplit=20,minbucket=7,repvar,valdat,indrep=0,thr=0.5){
lMit <- length(Eit)
lDit <- length(Pit)
if (N > 0){
  Ni <- N
}
if (is.na(Ni) == TRUE){
  Ni <- 99
}
lablarge <- names(table(data$class))[1]  # class 1 = large
labsmall <- names(table(data$class))[2]  # class 2 = small
nam1max <- list()
nam2max <- list()
bamax <- 0
tamax <- 0
D <- dim(data)[2] - 1
kumin <- length(stest)
summin <- length(stest)
mindlong <- list()
interp <- 0
############
for (d in Pit){
  set.seed(7654321)  ## for single equations
message("Number of predictors: ",d)
  for (n in 1:Ni){
#  set.seed(321*n + n) # 17:4:0.99: 0.7097482, 17:4:0.95: 0.709784  # changed for what reason?
# set.seed(54321*n + n) # not used anymore
# print(n)
    datat <- data
    datat$class <- NULL
    ind <- sample(1:D)[1:d]
    datat <- datat[,ind]
    datat <- cbind(datat,data$class)
    colnames(datat)[d+1] <- "class"
#print(n)
    if (N > 0){
      for (M in Eit){
#      print(M)
#      set.seed(7654321)
        for (i in 1:N){
#          ind1 <- sample(unique(name[check %in% labs[1,]]),prob=p1)[1:M] 
#          ind2 <- sample(unique(name[check %in% labs[2,]]),prob=p2)[1:M]
          ind1 <- sample(names(p1),prob=p1)[1:M] 
          ind2 <- sample(names(p2),prob=p2)[1:M]
          gTrain <- datat[(name %in% ind1 | name %in% ind2),]
          gTest <- datat[!(name %in% ind1 | name %in% ind2),]
          sTest <- stest[!(name %in% ind1 | name %in% ind2)]
          gname <- name[(name %in% ind1 | name %in% ind2)]
          gnameT <- name[!(name %in% ind1 | name %in% ind2)]
          gcheck <- check[(name %in% ind1 | name %in% ind2)]
          if (undersamp == TRUE & crit != "tal"){
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
            gTest <- datat[!(name %in% ind1 | name %in% ind2),]
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
          outAll <- PrInDTAll(zTTrain,class,ctestv=ctestv,conf.level=co,minsplit=minsplit,minbucket=minbucket)
#print(outAll$treeAll)
          if (outAll$interpAll == FALSE){
            interp <- interp + 1
          } 
          if (crit == "ba"){
            class.pred <- predict(outAll$treeAll,newdata=valdat)
            class.pred <- relevel(class.pred,ref=levels(valdat$class)[1])
          # sum(class.pred == data$class) / length(class.pred)
            nam1 <- list()
            nam2 <- list()
            if (indrep == 2){
        #   acc1E <- 0
        #  acc2E <- 0
        #  pred <- predict(out$modbest,newdata=data)
              ch <- table(repvar,class.pred)
        # class <- datain[,names(datain)==classname]
              conti <- cbind(table(repvar,valdat$class),ch)
              no1 <- sum(conti[,1] > 0)
              no2 <- sum(conti[,2] > 0)
              for (i in 1:dim(conti)[1]){
                if (conti[i,2] > 0 & conti[i,4] <= (conti[i,2]*thr)){ 
                  nam2 <- c(nam2,rownames(ch)[i])
                }
                if (conti[i,1] > 0 & conti[i,3] < (conti[i,1]*thr)){ 
                  nam1 <- c(nam1,rownames(ch)[i])
                }
              }
              ba2 <- 1- length(nam2) / no2
              ba1 <- 1 - length(nam1) / no1  
              ba <- (ba1 + ba2)/2
              if (length(nam1) == 0) {nam1 <- "-"}
              if (length(nam2) == 0) {nam2 <- "-"}
            } else {
              conti <- table(class.pred,valdat$class)
              if (dim(conti)[1] == 1){                                                ## NEWNEW
                conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
              }
              ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
              ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
              ba <- (ba1 + ba2)/2
            }
            if (ba > bamax & outAll$interpAll == FALSE){
              bamax <- ba
              acc1 <- ba1
              acc2 <- ba2
              ntmax <- nt
              dmax <- D
              indmax <- colnames(zTTrain)[1:d]
              ind1max <- ind1
              ind2max <- ind2
              nam1max <- nam1
              nam2max <- nam2
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
              indmax <- colnames(zTTrain)[1:d]
#print(indmax)
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
              indmax <- colnames(zTTrain)[1:d]
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
            labku <- unique(sTest[indstart])
            kustart <- length(labku)
            if ((kustart < kumin & outAll$interpAll == FALSE) |
                (kustart == kumin & ta > tamax & outAll$interpAll == FALSE) |
                (kustart == kumin & ta == tamax & sumlen < summin & outAll$interpAll == FALSE) ){
              kumin <- kustart
              labkumin <- labku
              summin <- sumlen
              mindstart <- indstart  ## NEWNEW
              mindlong <- indlong
#print(kumin)
#print(labkumin)
#print(ta)
              tamax <- ta
              bamax <- bat
              acc1 <- ta1
              acc2 <- ta2
              ntmax <- nt 
              dmax <- d
              indmax <- colnames(zTTrain)[1:d]
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
    } else {  ## N = 0
#print(n)
      gTrain <- datat
      zTTrain <- gTrain
      nt <- nrow(zTTrain)
      outAll <- PrInDTAll(zTTrain,class,ctestv=ctestv,conf.level=co,minsplit=minsplit,minbucket=minbucket)
      if (outAll$interpAll == FALSE){
        interp <- interp + 1
      }
#plot(outAll$treeAll)
      class.pred <- predict(outAll$treeAll,newdata=data)
#print(class.pred)
      class.pred <- relevel(class.pred,ref=levels(data$class)[1])
      conti <- table(class.pred,data$class)
      if (dim(conti)[1] == 1){                                                ## NEWNEW
        conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
      }
#### setting standards
      sTest <- stest
      gTest <- data            
      zTname <- name
      zTcheck <- check
      if (crit == "ba"){
        class.pred <- predict(outAll$treeAll,newdata=valdat)
        class.pred <- relevel(class.pred,ref=levels(valdat$class)[1])
      # sum(class.pred == data$class) / length(class.pred)
        nam1 <- list()
        nam2 <- list()
        if (indrep == 2){
        #  acc1E <- 0
        #  acc2E <- 0
        #  pred <- predict(out$modbest,newdata=data)
          ch <- table(repvar,class.pred)
        # class <- datain[,names(datain)==classname]
          conti <- cbind(table(repvar,valdat$class),ch)
          no1 <- sum(conti[,1] > 0)
          no2 <- sum(conti[,2] > 0)
          for (i in 1:dim(conti)[1]){
            if (conti[i,2] > 0 & conti[i,4] <= (conti[i,2]*thr)){ 
              nam2 <- c(nam2,rownames(ch)[i])
            }
            if (conti[i,1] > 0 & conti[i,3] < (conti[i,1]*thr)){ 
              nam1 <- c(nam1,rownames(ch)[i])
            }
          }
          ba2 <- 1- length(nam2) / no2
          ba1 <- 1 - length(nam1) / no1  
          ba <- (ba1 + ba2)/2
          if (length(nam1) == 0) {nam1 <- "-"}
          if (length(nam2) == 0) {nam2 <- "-"}
        } else {
          conti <- table(class.pred,valdat$class)
          if (dim(conti)[1] == 1){                                                ## NEWNEW
            conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
          }
          ba1 <- conti[1,1] / (conti[1,1] + conti[2,1])
          ba2 <- conti[2,2] / (conti[1,2] + conti[2,2])
          ba <- (ba1 + ba2)/2
        }
        if (ba > bamax & outAll$interpAll == FALSE){
          bamax <- ba
          acc1 <- ba1
          acc2 <- ba2
          ntmax <- nt
          dmax <- D
          indmax <- colnames(zTTrain)[1:D]
          ind1max <- ind1
          ind2max <- ind2
          nam1max <- nam1
          nam2max <- nam2
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
          indmax <- colnames(zTTrain)[1:d]
#print(indmax)
          ind1max <- indmax
          ind2max <- indmax
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
          indmax <- colnames(zTTrain)[1:d]
          ind1max <- indmax
          ind2max <- indmax
          outmax <- outAll
          gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
          colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
          gmaxTest <- gTest
        }
      }
      if (crit == "tal"){
        splitTest <- table(sTest) ## no element selection! 
        for (i in 2:length(splitTest)){
          splitTest[i] <- splitTest[i] + splitTest[i-1]  #### splits according to sTest
        }
        class.pred <- predict(outAll$treeAll,newdata=gTest)
        class.pred <- relevel(class.pred,ref=levels(gTest$class)[1])
        conti <- table(class.pred,gTest$class)
#print(conti)
        if (dim(conti)[1] == 1){                                                ## NEWNEW
          conti <- matrix(c(conti,0,0),nrow=2,byrow=TRUE)
        }
        ta1 <- conti[1,1] / (conti[1,1] + conti[2,1])
        ta2 <- conti[2,2] / (conti[1,2] + conti[2,2])
        bat <- (ta1 + ta2)/2
        if (bat == 0.5){
          kustart <- kumin + 1  # whole small class is wrongly predicted
        } else {
          ta <- sum(class.pred == gTest$class) / length(class.pred)
          predind <- class.pred == gTest$class
          indfal <- which(predind == FALSE)
          indfal <- indfal[!(indfal %in% splitTest)]
          indsplit <- split(indfal,cumsum(c(1, diff(indfal) != 1))) 
#print(length(indsplit))
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
          labku <- unique(sTest[indstart])
          kustart <- length(labku)
          if ((kustart < kumin & outAll$interpAll == FALSE) |
              (kustart == kumin & ta > tamax & outAll$interpAll == FALSE) |
              (kustart == kumin & ta == tamax & sumlen < summin & outAll$interpAll == FALSE) ){
            kumin <- kustart
            summin <- sumlen
            labkumin <- labku
            mindstart <- indstart  ## NEWNEW
            mindlong <- indlong
#print(kumin)
#print(labkumin)
#print(ta)
            tamax <- ta
            bamax <- bat
            acc1 <- ta1
            acc2 <- ta2
            ntmax <- nt 
            dmax <- d
            indmax <- colnames(zTTrain)[1:d]
            ind1max <- indmax
            ind2max <- indmax
            outmax <- outAll
            gmaxTrain <- cbind(zTTrain,as.vector(zTname),as.vector(zTcheck))
            colnames(gmaxTrain) <- c(colnames(zTTrain),"SubStruc","check")
            gmaxTest <- data
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
    elems <- labkumin
    message("current optimal test accuracy: ",round(tamax,digits=4),"\n")
  }
#if (crit == "ba"){message("current optimal balanced accuracy: ",round(bamax,digits=4))}
#if (crit == "bat"){message("current optimal balanced test accuracy: ",round(bamax,digits=4))}
#if (crit == "ta"){message("current optimal test accuracy: ",round(tamax,digits=4))} # {print(tamax)}
#if (crit == "tal"){message("current optimal number of wrongly predicted test squencies longer than ",ktest,": ",round(kumin))} # {print(kumin)}
}
if (N > 0){
  interp <- c(interp,(lMit*N*Ni*lDit))
} else {
  interp <- c(interp,(Ni*lDit))
}
result <- list(interp=interp,dmax=dmax,ntmax=ntmax,acc1=acc1,acc2=acc2,nam1=nam1max,nam2=nam2max,
     bamax=bamax,tamax=tamax,kumin=kumin,elems=elems,mindlong=mindlong,
     ind1max=ind1max,ind2max=ind2max,indmax=indmax,modbest = outmax$treeAll,gmaxTrain=gmaxTrain,gmaxTest=gmaxTest)
}