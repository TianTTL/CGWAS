step1 <- function(cgwasenv) {
  # StatExtraction

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  naid <- unique(foreach(i=1:cgwasenv$.TRAIT_NUM,.combine="c",.inorder=F,.export='StatE1') %dopar% StatE1(i, cgwasenv))
  printTmp <- foreach(i=1:cgwasenv$.TRAIT_NUM,.inorder=F,.export='StatE2') %dopar% StatE2(i, naid, cgwasenv)

  write.table(naid, file.path(cgwasenv$.CGWAS_COLDATA_PATH, 'NA_list'), row.names = F, col.names = F)

  stopCluster(sb)
}

step2 <- function(cgwasenv) {
  # InflationCorr

  normea <- integrate(ga,qchisq(cgwasenv$.CCP,1,lower.tail=F),qchisq(1,1),cgwasenv)$value
  inseq <- ppoints(20000)

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  infm <- foreach(i=1:cgwasenv$.TRAIT_NUM,.combine="rbind",.inorder=T,.export='calinf') %dopar% calinf(i,inseq,normea,cgwasenv)

  colnames(infm) <- c("UnCorrX2","Lambda","EstInf","CorrX2","CorrX2Ccp","EstEffProp")
  write.table(signif(infm,5),file.path(cgwasenv$.CGWAS_INFCOR_PATH, "InflationStat.txt"),row.names=F,quote=F)

  stopCluster(sb)
}


step3 <- function(cgwasenv) {
  # CorEstimation

  pairma <- t(combn(cgwasenv$.TRAIT_NUM, 2))
  n <- nrow(pairma)

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  tresm <- foreach(i=1:n,.combine="rbind",.inorder=T,.export='CorE') %dopar% CorE(i, pairma, cgwasenv)
  tresm <- matrix(tresm, nrow=n)

  corm <- t(combn(cgwasenv$.TRAIT_NUM, 2))
  n <- nrow(corm)
  corm <- cbind(corm,signif(tresm, 5))
  colnames(corm) <- c("Trait1","Trait2","GeneralCor","EstBgCor","EstEffCor")
  write.table(corm,
              file.path(cgwasenv$.CGWAS_INFCOR_PATH, "EstCorrelation.txt"),
              row.names=F, quote=F)

  stopCluster(sb)
}

step4 <- function(cgwasenv) {
  # preGMA

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  normea <- integrate(ga,qchisq(cgwasenv$.CCP,1,lower.tail=F),qchisq(1,1),cgwasenv)$value

  weieffs.tmp <- read.table(file.path(cgwasenv$.CGWAS_INFCOR_PATH, "InflationStat.txt"),
                                   header=T,stringsAsFactors=F)
  weieffs <- as.numeric(weieffs.tmp[,4])
  weieffsccp <- as.numeric(weieffs.tmp[,5])

  orgeffs <- cgwasenv$.TRAIT_EFFECT_SIZE; orgeffslist <- as.list(orgeffs)
  dirclist <- rep(1,cgwasenv$.TRAIT_NUM); dirclist <- as.list(dirclist)
  coeflist <- rep(1,cgwasenv$.TRAIT_NUM); coeflist <- as.list(coeflist)

  dropm <- foreach(i=1:cgwasenv$.TRAIT_NUM,.inorder=F,.export='copyfun.tcomb.gmeta') %dopar% copyfun.tcomb.gmeta(i, cgwasenv)

  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_INFCOR_PATH, "EstCorrelation.txt"),
                                 header=T,stringsAsFactors=F)[,c(1:2,4:5,3)])
  orgcorm <- diag(cgwasenv$.TRAIT_NUM)
  n <- 0
  for(i in 2:nrow(orgcorm)){
    orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
    n <- n+nrow(orgcorm)-i+1
  }

  motjv <- abs(coridm[,3])

  n <- cgwasenv$.TRAIT_NUM+1
  roundn <- 1
  existid <- rep(1,cgwasenv$.TRAIT_NUM)
  rec <- c()
  idlist <- as.list(1:cgwasenv$.TRAIT_NUM)

  if(!all(motjv^2<0.64)){
    orgresm <- as.data.frame(foreach(i=1:cgwasenv$.TRAIT_NUM,.combine="cbind",.inorder=T,.export='readtm.gmeta') %dopar%
                               readtm.gmeta(i, cgwasenv))

    while(!all(motjv^2<0.64)){
      mtarid <- order(motjv,decreasing=T)[1]
      selidv <- coridm[mtarid,1:2]
      bgc <- coridm[mtarid,3]
      dirclist[[n]] <- as.numeric(c(dirclist[[selidv[1]]],bgc/abs(bgc)*dirclist[[selidv[2]]]))
      idlist[[n]] <- as.numeric(c(idlist[[selidv[1]]],idlist[[selidv[2]]]))
      tsolcorm <- solve(orgcorm[idlist[[n]],idlist[[n]]])

      wv <- orgeffs[idlist[[n]]]
      coef <- wv/sqrt(sum(wv^2))
      pcoef <- metaf(wv*dirclist[[n]],tsolcorm)
      pmv <- (pcoef^2)/coef^2
      newv <- orgresm[,n] <- as.matrix(orgresm[,idlist[[n]]])%*%pcoef
      pnewv <- newv^2
      weieffsccp[n] <- signif(mean(pnewv[which(pnewv>qchisq(1-cgwasenv$.CCP,1))])/normea,5)
      e5 <- signif((sum((weieffsccp[idlist[[n]]]-1)*pmv)+1),5)
      coeflist[[n]] <- coef <- signif(pcoef,5)
      orgeffs[n] <- round(sum(pmv*wv))
      weieffs[n] <- signif(mean(pnewv),5)

      trec <- paste0("Round ",roundn)
      roundn <- roundn+1
      rec <- append(rec,trec)
      print(trec)

      trec <- paste0(selidv[1]," ",selidv[2]," ",orgeffs[selidv[1]],
                     "(",weieffs[selidv[1]]," ",weieffsccp[selidv[1]],") ",
                     orgeffs[selidv[2]],
                     "(",weieffs[selidv[2]]," ",weieffsccp[selidv[2]],")")
      rec <- append(rec,trec)
      print(trec)

      trec <- paste0("Pmeta ",weieffsccp[n]," ",e5,
                     " (",weieffsccp[selidv[1]]," ",weieffsccp[selidv[2]],") ")
      rec <- append(rec,trec)
      print(trec)

      trec <- paste0(paste(idlist[[selidv[1]]],collapse="-")," ",
                     paste(idlist[[selidv[2]]],collapse="-")," -> ",
                     n," ",orgeffs[n],"(",weieffs[n]," ",weieffsccp[n],") | ",
                     paste0(coef,collapse=" ")," | ",bgc)
      rec <- append(rec,trec)
      print(trec)

      trec <- paste0(" ")
      rec <- append(rec,trec)
      print(trec)

      write.table(as.data.frame(signif(newv,6)),
                         file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(n,".gmeta")),
                         row.names=F,col.names=F,quote=F)

      existid[n] <- 1
      existid[selidv] <- 0

      newcornum <- which(existid[-n]!=0)
      if(length(newcornum)==0){
        break
      }

      if(weieffs[n]<1.02){
        newv <- newv*sqrt(1.02/weieffs[n])
      }
      newv2 <- newv^2
      newv2m <- var(newv)

      newcoridm <- foreach(i=newcornum,.combine="rbind",.inorder=F,.export=c('multiCorr')) %dopar%
        multiCorr(i,n,newv,newv2,newv2m,orgeffs,weieffs,cgwasenv)
      deleteid <- union(union(union(which(coridm[,1]==selidv[1]),which(coridm[,1]==selidv[2])),which(coridm[,2]==selidv[1])),which(coridm[,2]==selidv[2]))
      coridm <- coridm[-deleteid,]
      motjv <- motjv[-deleteid]

      if(length(nrow(newcoridm))!=0){
        motjv <- c(motjv,abs(newcoridm[,3]))
      } else{
        motjv <- c(motjv,abs(newcoridm[3]))
      }
      coridm <- rbind(coridm,newcoridm)
      n <- n+1
    }

    for(i in which(existid==0)){
      file.remove(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(i,".gmeta")))
    }

    for(n in which(existid!=0)){
      if(length(idlist[[n]])!=1){
        if(weieffsccp[n]<1){
          ep <- 1/weieffsccp[n]-1
          ev <- 1+ep
          pnewv <- ev*orgresm[,n]^2
          weieffsccp[n] <- signif(mean(pnewv[which(pnewv>qchisq(1-cgwasenv$.CCP,1))])/normea,5)
          while(weieffsccp[n]<1){
            ev <- ev+ep
            pnewv <- ev*orgresm[,n]^2
            weieffsccp[n] <- signif(mean(pnewv[which(pnewv>qchisq(1-cgwasenv$.CCP,1))])/normea,5)
          }
          weieffs[n] <- signif(mean(pnewv),5)
          orgresm[,n] <- orgresm[,n]*sqrt(ev)
          fwrite(as.data.frame(signif(orgresm[,n],6)),
                             file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(n,".gmeta")),
                             row.names=F,col.names=F,quote=F)
        }
      }
    }
  }

  write.table(rec,
              file.path(cgwasenv$.CGWAS_GMA_PATH, "Record0.txt"),
              row.names=F,col.names=F,quote=F)

  trait.name.new <- outid <- which(existid!=0)
  for(i in 1:length(outid)){
    trait.name.new[i] <- paste(cgwasenv$.TRAIT_NAME[idlist[[outid[i]]]],collapse="-")
  }
  write.table(as.data.frame(cbind(1:length(outid),outid,trait.name.new,orgeffs[outid],weieffs[outid],weieffsccp[outid])),
              file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaPreCombination"),
              row.names=F,col.names=F,quote=F)

  cnv <- c()
  for(i in 1:length(outid)){
    if(outid[i]!=i){
      file.rename(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(outid[i],".gmeta")),
                  file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(i,".gmeta")))
    }
    coridm[coridm==outid[i]] <- i
    cnv <- c(cnv,length(idlist[[outid[i]]]))
  }

  ccm <- matrix(NA,length(outid),7*max(cnv)+3)
  for(i in 1:length(outid)){
    ccm[i,1:cnv[i]] <- idlist[[outid[i]]]
    ccm[i,(max(cnv)+1):(max(cnv)+cnv[i])] <- cgwasenv$.TRAIT_DATA_NAME[idlist[[outid[i]]]]
    ccm[i,(2*max(cnv)+1):(2*max(cnv)+cnv[i])] <- dirclist[[outid[i]]]
    ccm[i,(3*max(cnv)+1):(3*max(cnv)+cnv[i])] <- coeflist[[outid[i]]]
    ccm[i,(4*max(cnv)+1):(4*max(cnv)+cnv[i])] <- orgeffs[idlist[[outid[i]]]]
    ccm[i,(5*max(cnv)+1):(5*max(cnv)+cnv[i])] <- weieffs[idlist[[outid[i]]]]
    ccm[i,(6*max(cnv)+1):(6*max(cnv)+cnv[i])] <- weieffsccp[idlist[[outid[i]]]]
    ccm[i,(7*max(cnv)+1)] <- orgeffs[outid[i]]
    ccm[i,(7*max(cnv)+2)] <- weieffs[outid[i]]
    ccm[i,(7*max(cnv)+3)] <- weieffsccp[outid[i]]
  }

  write.table(as.data.frame(ccm),
              file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaCoefficient0"),
              row.names=F,col.names=F,quote=F)

  coridm <- matrix(coridm[order(coridm[,2]),], nrow = nrow(coridm))
  coridm <- matrix(coridm[order(coridm[,1]),], nrow = nrow(coridm))
  coridm <- as.data.frame(coridm)

  colnames(coridm)[1:2] <- c("Combination1","Combination2")
  write.table(coridm,
              file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta0.corest"),
              row.names=F,quote=F)

  stopCluster(sb)
}

step5 <- function(cgwasenv) {
  # GMA

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  precomb.tmp <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaPreCombination"),
                            header=F,stringsAsFactors=F)
  trait.es.new <- as.numeric(precomb.tmp[,4])
  trait.num.new <- length(trait.es.new)
  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta0.corest"),
                                 header=T,stringsAsFactors=F))

  orgcorm <- diag(trait.num.new)
  n <- 0
  for(i in 2:nrow(orgcorm)){
    orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
    n <- n+nrow(orgcorm)-i+1
  }

  tpov <- coridm[,5]-coridm[,3]

  nn <- 0
  while(!all(eigen(orgcorm)$values>0)){
    nn <- nn+1
    coridm[,3] <- coridm[,3]+0.1*tpov
    orgcorm <- diag(trait.num.new)
    n <- 0
    for(i in 2:nrow(orgcorm)){
      orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
      n <- n+nrow(orgcorm)-i+1
    }
  }

  print(nn)

  effcorm <- diag(trait.num.new)
  n <- 0
  for(i in 2:nrow(effcorm)){
    effcorm[i:nrow(effcorm),i-1] <- effcorm[i-1,i:nrow(effcorm)] <- coridm[(n+1):(n+nrow(effcorm)-i+1),4]
    n <- n+nrow(effcorm)-i+1
  }

  sigeffcorm <- cbind(orgcorm/abs(orgcorm),effcorm/abs(effcorm))

  orgresm <- foreach(i=1:trait.num.new,.combine="cbind",.inorder=T,.export='readtm.gmeta') %dopar%
             readtm.gmeta(i, cgwasenv)

  localsep <- 2e4

  gmetam <- foreach(i=iter(orgresm,by="row",chunksize=localsep),.combine="rbind",.inorder=T,.export='metaf.GMA') %dopar%
            metaf.GMA(i, sigeffcorm, trait.es.new, orgcorm, trait.num.new, cgwasenv)

  fwrite(as.data.frame(signif(gmetam[,1],4)),
                     file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta.p"),
                     sep=" ",row.names=F,col.names=F,quote=F)
  gmetam <- gmetam[,-1]
  fwrite(as.data.frame(gmetam),
                     file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta.dirc"),
                     sep=" ",na="NA",row.names=F,col.names=F,quote=F)

  qfswm <- foreach(i=iter(orgresm,by="row",chunksize=localsep),.combine="rbind",.inorder=T,.export='qfsw') %dopar%
            qfsw(i, orgcorm, trait.num.new, cgwasenv)

  fwrite(as.data.frame(qfswm[,(trait.num.new+1):ncol(qfswm)]),
                     file.path(cgwasenv$.CGWAS_GMA_PATH, "QFsw.ord"),
                     sep=" ",na="NA",row.names=F,col.names=F,quote=F)
  qfswm <- qfswm[,-((trait.num.new+1):ncol(qfswm))]
  fwrite(as.data.frame(signif(qfswm,4)),
                     file.path(cgwasenv$.CGWAS_GMA_PATH, "QFsw.ordp"),
                     sep=" ",na="NA",row.names=F,col.names=F,quote=F)
  minp <- foreach(i=iter(qfswm,by="row",chunksize=localsep),.combine="c",.inorder=T,.export='apminfun') %dopar%
            apminfun(i)
  fwrite(as.data.frame(signif(minp,4)),
                     file.path(cgwasenv$.CGWAS_GMA_PATH, "QFsw.p"),
                     row.names=F,col.names=F,quote=F)

  stopCluster(sb)
}

step6 <- function(cgwasenv) {
  # GMA+QFSW-2

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  normea <- integrate(ga,qchisq(cgwasenv$.CCP,1,lower.tail=F),qchisq(1,1),cgwasenv)$value

  effs.tmp <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaPreCombination"),
                     header=F,stringsAsFactors=F)
  orgeffs <- as.numeric(effs.tmp[,4])
  weieffs <- as.numeric(effs.tmp[,5])
  weieffsccp <- as.numeric(effs.tmp[,6])

  trait.num.new <- nrow(effs.tmp)

  dropm <- foreach(i=1:trait.num.new,.inorder=F,.export='copyfun.gmeta.gmetad') %dopar% copyfun.gmeta.gmetad(i, cgwasenv)

  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta0.corest"),
                                 header=T,stringsAsFactors=F))

  orgcorm <- diag(trait.num.new)
  n <- 0
  for(i in 2:nrow(orgcorm)){
    orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
    n <- n+nrow(orgcorm)-i+1
  }

  tpov <- coridm[,5]-coridm[,3]

  nn <- 0
  while(!all(eigen(orgcorm)$values>0)){
    nn <- nn+1
    coridm[,3] <- coridm[,3]+0.1*tpov
    orgcorm <- diag(trait.num.new)
    n <- 0
    for(i in 2:nrow(orgcorm)){
      orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
      n <- n+nrow(orgcorm)-i+1
    }
  }

  print(nn)

  motjv <- abs(coridm[,3]-coridm[,4])

  n <- trait.num.new+1
  roundn <- 1
  existid <- rep(1,trait.num.new)
  rec <- c()
  snpn <- length(as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(1,".gmetad")),
                                 header=F))[,1])
  dm <- matrix(NA,snpn,2)

  combstrlist <- list()
  for(i in 1:trait.num.new){
    combstrlist[[i]] <- matrix(c(i,1,1,1,orgeffs[i],weieffs[i],weieffsccp[i]),1)
  }

  localsep <- 2e4

  while(!all(motjv^2<0.04)){
    mtarid <- order(motjv,decreasing=T)[1]
    selidv <- coridm[mtarid,1:2]
    bgc <- coridm[mtarid,3]
    gc <- coridm[mtarid,4]
    tsolcorm <- solve(matrix(c(1,bgc,bgc,1),2))

    dm[,1] <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(selidv[1],".gmetad"))
                              ,header=F))[,1]
    dm[,2] <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(selidv[2],".gmetad"))
                              ,header=F))[,1]
    dm2 <- dm^2
    e1 <- mean(dm2[which(dm2[,1]>qchisq(1-cgwasenv$.CCP,1)),1])/normea
    e2 <- mean(dm2[which(dm2[,2]>qchisq(1-cgwasenv$.CCP,1)),2])/normea
    if(e2>e1){
      selidv <- selidv[2:1]
      dm[,1:2] <- dm[,2:1]
      et <- e1
      e1 <- e2
      e2 <- et
    }

    wv <- orgeffs[selidv]
    pcoef <- metaf(wv*c(1,bgc/abs(bgc)),tsolcorm)
    ncoef <- metaf(wv*c(1,-bgc/abs(bgc)),tsolcorm)
    pnewv <- (dm%*%pcoef)^2
    e3 <- mean(pnewv[which(pnewv>qchisq(1-cgwasenv$.CCP,1))])/normea
    nnewv <- (dm%*%ncoef)^2
    e4 <- mean(nnewv[which(nnewv>qchisq(1-cgwasenv$.CCP,1))])/normea

    coef <- wv/sqrt(sum(wv^2))
    pmv <- (pcoef^2)/coef^2
    nmv <- (ncoef^2)/coef^2
    e5 <- (sum(c(e1-1,e2-1)*pmv)*cgwasenv$.HETEROPROP+1)
    e6 <- (sum(c(e1-1,e2-1)*nmv)*cgwasenv$.HETEROPROP+1)

    if(roundn!=1){
      trec <- paste0(" ")
      rec <- append(rec,trec)
      print(trec)
    }
    trec <- paste0("Round ",roundn)
    rec <- append(rec,trec)
    print(trec)
    roundn <- roundn+1
    trec <- paste0(selidv[1]," ",selidv[2]," ",orgeffs[selidv[1]],"(",weieffs[selidv[1]]," ",signif(e1,5),") ",orgeffs[selidv[2]],"(",weieffs[selidv[2]]," ",signif(e2,5),")")
    rec <- append(rec,trec)
    print(trec)
    if((e5<e3)&(e1<e3)&(e2<e3)&(e6<e4)&(e2<e4)&(e1<e4)){
      motjv[mtarid] <- 0
      trec <- paste0("Pmeta ",signif(e3,5),">",signif(e5,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),") Nmeta ",signif(e4,5),">",signif(e6,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),")")
      rec <- append(rec,trec)
      print(trec)
      trec <- paste0("Decision: 1 Nonmeta | ",bgc," ",gc)
      rec <- append(rec,trec)
      print(trec)
      next
    } else if((e5<e3)&(e1<e3)&(e2<e3)){
      orgeffs[n] <- round(sum(pmv*wv))
      weieffs[n] <- signif(mean(pnewv),5)
      te <- signif(e3,5)
      newv <- dm%*%pcoef
      coef <- signif(pcoef,5)
      trec <- paste0("Pmeta ",signif(e3,5),">",signif(e5,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),") Nmeta ",signif(e4,5),"<",signif(e6,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),")")
      rec <- append(rec,trec)
      print(trec)
      trec <- paste0("Decision: 2 Pmeta ",selidv[1]," ",selidv[2]," -> ",n," ",orgeffs[n],"(",weieffs[n]," ",signif(e3,5),") | ",coef[1]," ",coef[2]," | ",bgc," ",gc)
      rec <- append(rec,trec)
      print(trec)
    } else if((e6<e4)&(e2<e4)&(e1<e4)){
      orgeffs[n] <- round(sum(nmv*wv))
      weieffs[n] <- signif(mean(nnewv),5)
      te <- signif(e4,5)
      newv <- dm%*%ncoef
      coef <- signif(ncoef,5)
      trec <- paste0("Pmeta ",signif(e3,5),"<",signif(e5,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),") Nmeta ",signif(e4,5),">",signif(e6,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),")")
      rec <- append(rec,trec)
      print(trec)
      trec <- paste0("Decision: 3 Nmeta ",selidv[1]," ",selidv[2]," -> ",n," ",orgeffs[n],"(",weieffs[n]," ",signif(e4,5),") | ",coef[1]," ",coef[2]," | ",bgc," ",gc)
      rec <- append(rec,trec)
      print(trec)
    } else{
      motjv[mtarid] <- 0
      trec <- paste0("Pmeta ",signif(e3,5),"<",signif(e5,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),") Nmeta ",signif(e4,5),"<",signif(e6,5),"[",signif(cgwasenv$.HETEROPROP,5),"x] (",signif(e1,5)," ",signif(e2,5),")")
      rec <- append(rec,trec)
      print(trec)
      trec <- paste0("Decision: 4 Nonmeta | ",bgc," ",gc)
      rec <- append(rec,trec)
      print(trec)
      next
    }

    fwrite(as.data.frame(signif(newv,6)),
                       file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(n,".gmetad"))
                       ,row.names=F,col.names=F,quote=F)
    file.remove(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(selidv[1],".gmetad")))
    file.remove(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(selidv[2],".gmetad")))

    existid[n] <- 1
    existid[selidv] <- 0

    tempc1 <- combstrlist[[selidv[1]]]
    tempc2 <- combstrlist[[selidv[2]]]
    tempc1[,2] <- tempc1[,2]*coef[1]
    tempc2[,2] <- tempc2[,2]*coef[2]
    combstrlist[[n]] <- rbind(rbind(tempc1,tempc2),matrix(c(n,NA,1,1,orgeffs[n],weieffs[n],te),1))

    newcornum <- which(existid[-n]!=0)
    if(length(newcornum)==0){
      break
    }

    if(weieffs[n]<1.02){
      newv <- newv*sqrt(1.02/weieffs[n])
    }
    newv2 <- newv^2
    newv2m <- var(newv)

    newcoridm <- foreach(i=newcornum,.combine="rbind",.inorder=F,.export='multiCorr.gmetad') %dopar%
                 multiCorr.gmetad(i,n,newv,newv2,newv2m,orgeffs,weieffs,cgwasenv)
    deleteid <- union(union(union(which(coridm[,1]==selidv[1]),which(coridm[,1]==selidv[2])),which(coridm[,2]==selidv[1])),which(coridm[,2]==selidv[2]))
    coridm <- coridm[-deleteid,]
    motjv <- motjv[-deleteid]

    if(length(nrow(newcoridm))!=0){
      motjv <- c(motjv,abs(newcoridm[,3]-newcoridm[,4]))
    } else{
      motjv <- c(motjv,abs(newcoridm[3]-newcoridm[4]))
    }
    coridm <- rbind(coridm,newcoridm)
    n <- n+1
  }
  write.table(rec,
              file.path(cgwasenv$.CGWAS_GMA_PATH, "Record2.txt"),
              row.names=F,col.names=F,quote=F)

  outid <- which(existid!=0)
  write.table(as.data.frame(cbind(1:length(outid),outid)),
              file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaCombination"),
              row.names=F,col.names=F,quote=F)

  cnv <- c()
  combinfo <- c()
  for(i in 1:length(outid)){
    if(outid[i]!=i){
      file.rename(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(outid[i],".gmetad")),
                  file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(i,".gmetad")))
    }
    coridm[coridm==outid[i]] <- i
    combstrlist[[outid[i]]][,2] <- signif(combstrlist[[outid[i]]][,2],5)
    combinfo <- rbind(combinfo,c(combstrlist[[outid[i]]][nrow(combstrlist[[outid[i]]]),5:7]))
    combstrlist[[outid[i]]] <- na.omit(combstrlist[[outid[i]]])
    cnv <- c(cnv,nrow(combstrlist[[outid[i]]]))
  }

  ccm <- matrix(NA,length(outid),5*max(cnv))
  for(i in 1:length(outid)){
    ccm[i,1:cnv[i]] <- combstrlist[[outid[i]]][,1]
    ccm[i,(max(cnv)+1):(max(cnv)+cnv[i])] <- combstrlist[[outid[i]]][,2]
    ccm[i,(2*max(cnv)+1):(2*max(cnv)+cnv[i])] <- combstrlist[[outid[i]]][,5]
    ccm[i,(3*max(cnv)+1):(3*max(cnv)+cnv[i])] <- combstrlist[[outid[i]]][,6]
    ccm[i,(4*max(cnv)+1):(4*max(cnv)+cnv[i])] <- combstrlist[[outid[i]]][,7]
  }
  ccm <- cbind(ccm,combinfo)

  write.table(as.data.frame(ccm),
              file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaCoefficient2"),
              row.names=F,col.names=F,quote=F)

  coridm <- matrix(coridm[order(coridm[,2]),], nrow = nrow(coridm))
  coridm <- matrix(coridm[order(coridm[,1]),], nrow = nrow(coridm))
  coridm <- as.data.frame(coridm)

  colnames(coridm)[1:2] <- c("Combination1","Combination2")
  write.table(coridm,
              file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta2.corest"),
              row.names=F,quote=F)

  stopCluster(sb)
}

step7 <- function(cgwasenv) {
  # QFSW

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  combnum <- nrow(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaCombination"),
                             header=F,stringsAsFactors=F))

  coridm <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0("Gmeta2.corest")),
                       header=T,stringsAsFactors=F)
  gmetacorm <- diag(combnum)
  n <- 0
  if (nrow(gmetacorm) >= 2) {
    for(i in 2:nrow(gmetacorm)){
      gmetacorm[i:nrow(gmetacorm),i-1] <- gmetacorm[i-1,i:nrow(gmetacorm)] <- coridm[(n+1):(n+nrow(gmetacorm)-i+1),3]
      n <- n+nrow(gmetacorm)-i+1
    }
  }

  tpov <- coridm[,5]-coridm[,3]

  nn <- 0
  while(!all(eigen(gmetacorm)$values>0)){
    nn <- nn+1
    coridm[,3] <- coridm[,3]+0.1*tpov
    gmetacorm <- diag(combnum)
    n <- 0
    for(i in 2:nrow(gmetacorm)){
      gmetacorm[i:nrow(gmetacorm),i-1] <- gmetacorm[i-1,i:nrow(gmetacorm)] <- coridm[(n+1):(n+nrow(gmetacorm)-i+1),3]
      n <- n+nrow(gmetacorm)-i+1
    }
  }

  print(nn)

  orgresm <- foreach(i=1:combnum,.combine="cbind",.inorder=T,.export='readtm.gmetad') %dopar% readtm.gmetad(i,cgwasenv)
  orgresm <- as.matrix(orgresm)

  localsep <- 2e4

  qfswm <- foreach(i=iter(orgresm,by="row",chunksize=localsep),.combine="rbind",.inorder=T,.export='qfsw') %dopar%
           qfsw(i, gmetacorm, combnum, cgwasenv)

  fwrite(as.data.frame(qfswm[,(combnum+1):ncol(qfswm)]),
         file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaQFsw.ord"),
         sep=" ",na="NA",row.names=F,col.names=F,quote=F)
  if ((combnum+1) <= ncol(qfswm)) {
    qfswm <- as.matrix(qfswm[,-((combnum+1):ncol(qfswm))])
  }
  fwrite(as.data.frame(signif(qfswm,4)),
         file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaQFsw.ordp"),
         sep=" ",na="NA",row.names=F,col.names=F,quote=F)
  minp <- foreach(i=iter(qfswm,by="row",chunksize=localsep),.combine="c",.inorder=T,.export='apminfun') %dopar% apminfun(i)
  fwrite(as.data.frame(signif(minp,4)),
         file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaQFsw.p"),
         row.names=F,col.names=F,quote=F)

  stopCluster(sb)
}

step8 <- function(cgwasenv) {
  # Simulation

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  trait.es.new <- as.numeric(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaPreCombination"),
                                        header=F,stringsAsFactors=F)[,4])
  trait.num.new <- length(trait.es.new)
  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta0.corest"),
                                 header=T,stringsAsFactors=F))

  orgcorm <- diag(trait.num.new)
  n <- 0
  for(i in 2:nrow(orgcorm)){
    orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
    n <- n+nrow(orgcorm)-i+1
  }

  tpov <- coridm[,5]-coridm[,3]

  nn <- 0
  while(!all(eigen(orgcorm)$values>0)){
    nn <- nn+1
    coridm[,3] <- coridm[,3]+0.1*tpov
    orgcorm <- diag(trait.num.new)
    n <- 0
    for(i in 2:nrow(orgcorm)){
      orgcorm[i:nrow(orgcorm),i-1] <- orgcorm[i-1,i:nrow(orgcorm)] <- coridm[(n+1):(n+nrow(orgcorm)-i+1),3]
      n <- n+nrow(orgcorm)-i+1
    }
  }

  print(nn)

  effcorm <- diag(trait.num.new)
  n <- 0
  for(i in 2:nrow(effcorm)){
    effcorm[i:nrow(effcorm),i-1] <- effcorm[i-1,i:nrow(effcorm)] <- coridm[(n+1):(n+nrow(effcorm)-i+1),4]
    n <- n+nrow(effcorm)-i+1
  }

  sigeffcorm <- cbind(orgcorm/abs(orgcorm),effcorm/abs(effcorm))

  coridm <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0("Gmeta2.corest")),
                       header=T,stringsAsFactors=F)
  combnum <- nrow(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaCombination"),
                             header=F,stringsAsFactors=F))
  gmqfc <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaCoefficient2"),
                      header=F,stringsAsFactors=F)
  gmqfcnum <- (ncol(gmqfc)-3)/5

  swit <- FALSE
  if(combnum!=trait.num.new){
    swit <- TRUE
    corm2 <- diag(combnum)

    n <- 0
    if (2 <= nrow(corm2)) {
      for(i in 2:nrow(corm2)){
        corm2[i:nrow(corm2),i-1] <- corm2[i-1,i:nrow(corm2)] <- coridm[(n+1):(n+nrow(corm2)-i+1),3]
        n <- n+nrow(corm2)-i+1
      }
    }


    tpov <- coridm[,5]-coridm[,3]

    nn <- 0
    while(!all(eigen(corm2)$values>0)){
      nn <- nn+1
      coridm[,3] <- coridm[,3]+0.1*tpov
      corm2 <- diag(combnum)
      n <- 0
      for(i in 2:nrow(corm2)){
        corm2[i:nrow(corm2),i-1] <- corm2[i-1,i:nrow(corm2)] <- coridm[(n+1):(n+nrow(corm2)-i+1),3]
        n <- n+nrow(corm2)-i+1
      }
    }

    print(nn)
  }

  snpn <- 1e4
  simuloop <- 1e4
  thv <- sqrt(qchisq(1-cgwasenv$.CCP,1))

  system.time({
    simup <- foreach(i=1:simuloop,.combine="rbind",.inorder=F,.export='Simudata') %dopar%
             Simudata(i, swit, orgcorm, snpn, sigeffcorm, corm2, thv, trait.es.new, trait.num.new, combnum, gmqfc, gmqfcnum)
  })

  sdsd <- sd(-log10(qbeta(ppoints(simuloop),1,snpn)))

  mt <- qbeta(0.5,1,snpn)/apply(simup,2,quantile,probs=0.5)

  lgcrm <- -log10(t(simup)*mt)
  sdprop <- apply(lgcrm,1,sd)/sdsd

  crsimup <- t(10^(-((lgcrm-apply(lgcrm,1,mean))/sdprop+apply(lgcrm,1,mean)+(-log10(qbeta(0.5,1,snpn))-apply(lgcrm,1,mean))*(sdprop-1)/(sdprop))))

  mmt <- min(ncol(crsimup),qbeta(0.5,1,snpn)/quantile(apply(crsimup,1,min),probs=0.5))

  write.table(as.data.frame(mt),
              file.path(cgwasenv$.CGWAS_GMA_PATH, "MTT"),
              row.names=F,col.names=F,quote=F)
  write.table(as.data.frame(mmt),
              file.path(cgwasenv$.CGWAS_GMA_PATH, "MMT"),
              row.names=F,col.names=F,quote=F)

  stopCluster(sb)

}

step9 <- function(cgwasenv) {
  # Visualization

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  Sind <- as.data.frame(fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SnpIndex"),
                                          header=T))

  threspm <- as.numeric(unlist(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "MTT"),
                                          header=F,stringsAsFactors=F)))
  swit <- FALSE
  bonfv <- as.numeric(unlist(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH, "MMT"),
                                        header=F,stringsAsFactors=F)))
  if(length(threspm)==2){
    threspm[3] <- threspm[2]
    swit <- TRUE
  }

  gmetap <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH, "Gmeta.p"),
                                        header=F,stringsAsFactors=F))[,1]
  gmetap <- gmetap*threspm[1]
  gmetap[gmetap>1] <- 1

  qfswp <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH, "QFsw.p"),
                                       header=F,stringsAsFactors=F))[,1]
  qfswp <- qfswp*threspm[3]
  qfswp[qfswp>1] <- 1

  gmetaqfswp <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH, "GmetaQFsw.p"),
                                            header=F,stringsAsFactors=F))[,1]
  gmetaqfswp <- gmetaqfswp*threspm[2]
  gmetaqfswp[gmetaqfswp>1] <- 1

  pm3 <- cbind(gmetap,gmetaqfswp,qfswp)
  localsep <- 2e4
  pm3min <- foreach(i=iter(pm3,by="row",chunksize=localsep),.combine="c",.inorder=T,.export='apminfun') %dopar% apminfun(i)

  gmetares <- mandata(-log10(gmetap),cgwasenv$.MANCO,Sind)
  qfswres <- mandata(-log10(qfswp),cgwasenv$.MANCO,Sind)
  gmetaqfswres <- mandata(-log10(gmetaqfswp),cgwasenv$.MANCO,Sind)

  pm3minres <- mandata(-log10(pm3min),cgwasenv$.MANCO,Sind)

  orgtm <- foreach(i=1:cgwasenv$.TRAIT_NUM,.combine="cbind",.inorder=T,.export='readtm.tcomb') %dopar% readtm.tcomb(i, cgwasenv)

  orgminp <- foreach(i=iter(orgtm,by="row",chunksize=localsep),.combine="c",.inorder=T,.export='apmaxfun') %dopar% apmaxfun(i)

  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_INFCOR_PATH, "EstCorrelation.txt"),
                                 header=T,stringsAsFactors=F)[,c(1:2,4:5)])
  orgcorm <- diag(cgwasenv$.TRAIT_NUM)
  n <- 0
  for(i in 2:cgwasenv$.TRAIT_NUM){
    orgcorm[i:cgwasenv$.TRAIT_NUM,i-1] <- orgcorm[i-1,i:cgwasenv$.TRAIT_NUM] <- coridm[(n+1):(n+cgwasenv$.TRAIT_NUM-i+1),3]
    n <- n+cgwasenv$.TRAIT_NUM-i+1
  }

  orgminp <- pchisq(orgminp,1,lower.tail=F)
  orgminres <- mandata(-log10(orgminp),cgwasenv$.MANCO,Sind)


  odpm3min <- sort(pm3min)
  BHp3m <- p.adjust(odpm3min*bonfv,method="fdr")

  idpp3m <- realidp <- (nrow(Sind)/1e6)

  while((realidp<=idpp3m)&(realidp!=1)){
    idpp3m <- realidp
    n <- 0
    fpn <- 0
    while(fpn<idpp3m){
      n <- n+1
      fpn <- BHp3m[n]*n
      if(BHp3m[n]>0.05){
        break
      }
    }

    thresp3m <- mean(odpm3min[c(n-1,n)])
    if(n==1){
      fdrp3m <- BHp3m[1]
    } else{
      fdrp3m <- BHp3m[n-1]
    }

    temv <- fdrfun(pm3minres,thresp3m,Gt=5e-8/bonfv,cgwasenv)[c(1,3)]
    realidp <- temv[1]/temv[2]
    if((realidp<=idpp3m)&((sum(odpm3min<=thresp3m)*fdrp3m)<realidp)){
      idpp3m <- realidp
      break
    }
  }

  odorgminp <- sort(orgminp)
  BHorgminp <- p.adjust(odorgminp*effvar(orgcorm),method="fdr")

  idporg <- realidp <- (nrow(Sind)/1e6)

  while((realidp<=idporg)&(realidp!=1)){
    idporg <- realidp
    n <- 0
    fpn <- 0
    while(fpn<idporg){
      n <- n+1
      fpn <- BHorgminp[n]*n
      if(BHorgminp[n]>0.05){
        break
      }
    }

    thresorgmin <- mean(odorgminp[c(n-1,n)])
    if(n==1){
      fdrorg <- BHorgminp[1]
    } else{
      fdrorg <- BHorgminp[n-1]
    }

    temv <- fdrfun(orgminres,thresorgmin,Gt=5e-8/effvar(orgcorm),cgwasenv)[c(1,3)]
    realidp <- temv[1]/temv[2]
    if((realidp<=idporg)&((sum(odorgminp<=thresorgmin)*fdrorg)<realidp)){
      idporg <- realidp
      break
    }
  }

  res1 <- signif(cbind(c(5e-8/effvar(orgcorm),
                         thresorgmin,
                         fdrorg,
                         sum(odorgminp<=thresorgmin)*fdrorg/idporg,
                         fdrfun(orgminres,thresorgmin,Gt=5e-8/effvar(orgcorm),cgwasenv)),
                       c(5e-8/bonfv,
                         thresp3m,
                         fdrp3m,
                         sum(odpm3min<=thresp3m)*fdrp3m/idpp3m,
                         fdrfun(pm3minres,thresp3m,Gt=5e-8/bonfv,cgwasenv))),5)
  colnames(res1) <- c("GWASmin","CGWASmin")
  rownames(res1) <- c("Bonf-cf","Fdr-cf","Fdr","Exp-fp-idSNP","FdrSNP","SigSNP","FdrLoci","SigLoci","Fdr-only-Loci")
  write.csv(res1,
            file.path(cgwasenv$.CGWAS_RESULT_PATH, "Summary.csv"),
            quote=F)


  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "MergeMan.jpg"),
       width=7500,height=10000,res=600)
  par(mar=c(2,5,2,2))
  #layout(matrix(c(rep(c(rep(1,4),rep(2,4),rep(6,4)),4),rep(c(rep(3,3),rep(4,3),rep(5,3),rep(7,3)),3)),12,7))
  layout(matrix(c(rep(c(rep(1,4),rep(3,4),rep(6,4)),4),rep(c(rep(2,3),rep(4,3),rep(5,3),rep(7,3)),3)),12,7))
  manhattan(orgminres,orgminp,orgcorm,thresp3m,thresorgmin,pm3min,pm3,cex.axis=1,pt.cex=0.7,pt.col=c('lightblue2',rgb(102/255,190/255,1)),pch=20,St=thresorgmin,Gt=5e-8/effvar(orgcorm),Lm=cgwasenv$.MANCO,hl=2,locisep=cgwasenv$.LOCISEP)
  manhattan(pm3minres,orgminp,orgcorm,thresp3m,thresorgmin,pm3min,pm3,cex.axis=1,pt.cex=0.7,pt.col=c('lightblue2',rgb(102/255,190/255,1)),pch=20,St=thresp3m,Gt=5e-8/bonfv,Lm=cgwasenv$.MANCO,hl=4,locisep=cgwasenv$.LOCISEP)
  manhattan(pm3minres,orgminp,orgcorm,thresp3m,thresorgmin,pm3min,pm3,cex.axis=1,pt.cex=0.7,pt.col=c('lightblue2',rgb(102/255,190/255,1)),pch=20,St=thresp3m,Gt=5e-8/bonfv,Lm=cgwasenv$.MANCO,hl=1,locisep=cgwasenv$.LOCISEP)
  manhattan(gmetares,orgminp,orgcorm,thresp3m,thresorgmin,pm3min,pm3,cex.axis=1,pt.cex=0.7,pt.col=c('brown2','firebrick4'),pch=20,St=thresp3m,Gt=5e-8,Lm=cgwasenv$.MANCO,hl=3,locisep=cgwasenv$.LOCISEP)
  manhattan(gmetaqfswres,orgminp,orgcorm,thresp3m,thresorgmin,pm3min,pm3,cex.axis=1,pt.cex=0.7,pt.col=c('darkgoldenrod2','darkgoldenrod4'),pch=20,St=thresp3m,Gt=5e-8,Lm=cgwasenv$.MANCO,hl=3,locisep=cgwasenv$.LOCISEP)
  manbar(pm3minres,pm3,orgminp,orgtm,orgcorm,thresorgmin,orgminp,cex.axis=1,St=thresp3m,Gt=5e-8,Lm=cgwasenv$.MANCO,locisep=cgwasenv$.LOCISEP)
  manhattan(qfswres,orgminp,orgcorm,thresp3m,thresorgmin,pm3min,pm3,cex.axis=1,pt.cex=0.7,pt.col=c('limegreen','springgreen4'),pch=20,St=thresp3m,Gt=5e-8,Lm=cgwasenv$.MANCO,hl=3,locisep=cgwasenv$.LOCISEP)

  dev.off()

  stopCluster(sb)
}

step10 <- function(cgwasenv) {
  # Summary

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  Sind <- as.data.frame(fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH,"SnpIndex"),
                              header=T))

  threspm <- as.numeric(unlist(read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"MTT"),
                                          header=F,stringsAsFactors=F)))
  swit <- FALSE
  bonfv <- length(threspm)
  if(length(threspm)==2){
    threspm[3] <- threspm[2]
    swit <- TRUE
  }

  basedminfo <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaPreCombination"),
                           header=F,stringsAsFactors=F)
  advdminfo <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaCoefficient2"),
                          header=F,stringsAsFactors=F)

  gmetadirc <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"Gmeta.dirc"),
                                           header=F,stringsAsFactors=F))
  gmetatm <- foreach(i=1:length(basedminfo[,1]),.combine="cbind",.inorder=T,.export='readtm.gmeta') %dopar% readtm.gmeta(i,cgwasenv)
  gmetaqfswtm <- foreach(i=1:length(advdminfo[,1]),.combine="cbind",.inorder=T,.export='readtm.gmetad') %dopar% readtm.gmetad(i,cgwasenv)
  gmetaqfswtm <- as.matrix(gmetaqfswtm)

  qfsword <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"QFsw.ord"),
                                         header=F,stringsAsFactors=F))
  qfswordp <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"QFsw.ordp"),
                                          header=F,stringsAsFactors=F))

  gmetaqfsword <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaQFsw.ord"),
                                              header=F,stringsAsFactors=F))
  gmetaqfswordp <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaQFsw.ordp"),
                                               header=F,stringsAsFactors=F))

  gmetap <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"Gmeta.p"),
                                                 header=F,stringsAsFactors=F))[,1]
  qfswp <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"QFsw.p"),
                                       header=F,stringsAsFactors=F))[,1]
  gmetaqfswp <- as.matrix(fread(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaQFsw.p"),
                                            header=F,stringsAsFactors=F))[,1]
  orgpm3 <- cbind(gmetap,gmetaqfswp,qfswp)

  gmetap <- gmetap*threspm[1]
  gmetap[gmetap>1] <- 1
  qfswp <- qfswp*threspm[3]
  qfswp[qfswp>1] <- 1
  gmetaqfswp <- gmetaqfswp*threspm[2]
  gmetaqfswp[gmetaqfswp>1] <- 1

  pm3 <- cbind(gmetap,gmetaqfswp,qfswp)
  localsep <- 2e4
  pm3min <- foreach(i=iter(pm3,by="row",chunksize=localsep),.combine="c",.inorder=T,.export='apminfun') %dopar% apminfun(i)

  orgtm <- foreach(i=1:cgwasenv$.TRAIT_NUM,.combine="cbind",.inorder=T,.export='readtm.tcomb') %dopar% readtm.tcomb(i,cgwasenv)

  localsep <- 2e4
  orgminp <- foreach(i=iter(orgtm,by="row",chunksize=localsep),.combine="c",.inorder=T,.export='apmaxfun') %dopar% apmaxfun(i)

  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_INFCOR_PATH,"EstCorrelation.txt"),
                                 header=T,stringsAsFactors=F)[,c(1:2,4:5)])
  orgcorm <- diag(cgwasenv$.TRAIT_NUM)
  n <- 0
  for(i in 2:cgwasenv$.TRAIT_NUM){
    orgcorm[i:cgwasenv$.TRAIT_NUM,i-1] <- orgcorm[i-1,i:cgwasenv$.TRAIT_NUM] <- coridm[(n+1):(n+cgwasenv$.TRAIT_NUM-i+1),3]
    n <- n+cgwasenv$.TRAIT_NUM-i+1
  }

  orgminp <- pchisq(orgminp,1,lower.tail=F)

  st <- read.csv(file.path(cgwasenv$.CGWAS_RESULT_PATH,"Summary.csv"),
                 row.names=1,header=T,stringsAsFactors=F)

  orgseid <- order(orgminp)[1:st[5,1]]
  #orgseid <- order(orgminp)[1:st[5,2]]
  #tt <- Sind[orgseid,]
  #tt <- tt[order(tt[,2]),]
  #tt <- tt[order(tt[,1]),]
  #write.csv(tt,"../../Result/GWASmorehits.csv",row.names=F,quote=F)

  if(st[6,1]!=0){
    orgsig <- c(rep(1,st[6,1]),rep(2,st[5,1]-st[6,1]))
  } else{
    orgsig <- c(rep(2,st[5,1]))
  }

  forloci <- cbind(Sind[orgseid,],orgminp[orgseid])
  forloci <- forloci[order(forloci[,2]),]
  forloci <- forloci[order(forloci[,1]),]
  orglociv <- fdrfun.summary(forloci,cgwasenv)
  orginfom <- cbind(rownames(orglociv),orglociv[,-4],orgsig)

  p3mseid <- order(pm3min)[1:st[5,2]]
  if(st[6,2]!=0){
    p3msig <- c(rep(1,st[6,2]),rep(2,st[5,2]-st[6,2]))
  } else{
    p3msig <- c(rep(2,st[5,2]))
  }

  forloci <- cbind(Sind[p3mseid,],pm3min[p3mseid])
  forloci <- forloci[order(forloci[,2]),]
  forloci <- forloci[order(forloci[,1]),]
  p3mlociv <- fdrfun.summary(forloci,cgwasenv)
  p3minfom <- cbind(rownames(p3mlociv),p3mlociv[,-4],p3msig)

  allselid <- as.data.frame(union(orgseid,p3mseid))

  colnames(allselid) <- colnames(orginfom)[1] <- colnames(p3minfom)[1] <- "SummaryID"

  allinfom <- merge(merge(allselid,p3minfom,by="SummaryID",all.x=T,all.y=T,sort=F)[,-c(2:4)],
                    orginfom,by="SummaryID",all.x=T,all.y=T,sort=F)[,-c(4:6)]
  allinfom[is.na(allinfom[,3]),3] <- 3
  allinfom[is.na(allinfom[,5]),5] <- 3

  allselid <- allinfom[,1]
  op <- pchisq(orgtm[allselid,]^2,1,lower.tail=F)
  mt <- cgwasenv$.TRAIT_NAME[apply(op,1,function(a){return(order(a)[1])})]
  allinfom <- cbind(Sind[allselid,],allinfom,mt,
                    signif(as.data.frame(cbind(pm3min[allselid],pm3min[allselid]*bonfv,
                                               orgminp[allselid],orgminp[allselid]*effvar(orgcorm),
                                               orgpm3[allselid,],pm3[allselid,],orgtm[allselid,],
                                               op,gmetatm[allselid,],pchisq(gmetatm[allselid,]^2,1,lower.tail=F),
                                               gmetadirc[allselid,],qfsword[allselid,],qfswordp[allselid,],
                                               gmetaqfswtm[allselid,],pchisq(gmetaqfswtm[allselid,]^2,1,lower.tail=F),
                                               gmetaqfsword[allselid,],gmetaqfswordp[allselid,])),5))

  allinfom <- allinfom[order(allinfom[,4]),]

  colnames(allinfom) <- c("CHR","BP","SNP",
                          "SummaryID","Cp-loci","Cp-significance","Sp-loci","Sp-significance",
                          "MinTrait","Mincp-uncorr","Mincp-corr","Minsp-uncorr","Minsp-corr",
                          "Gm-uncorr","Gmqf-uncorr","Qf-uncorr","Gm-corr","Gmqf-corr","Qf-corr",
                          paste0("T_",cgwasenv$.TRAIT_NAME),paste0("P_",cgwasenv$.TRAIT_NAME),
                          paste0("T_gm_",basedminfo[,3]),paste0("P_gm_",basedminfo[,3]),
                          paste0("D_gm_",basedminfo[,3]),paste0("O_qf_",1:length(basedminfo[,1])),
                          paste0("Cp_qf_",1:length(basedminfo[,1])),paste0("T_gmqf_",1:length(advdminfo[,1])),
                          paste0("P_gmqf_",1:length(advdminfo[,1])),paste0("O_gmqf_",1:length(advdminfo[,1])),
                          paste0("Cp_gmqf_",1:length(advdminfo[,1])))

  cpod <- cpnv <- allinfom[,5]
  for(i in 1:max(allinfom[,5],na.rm=T)){
    ttid <- which(allinfom[,5]==i)
    cpod[ttid] <- rank(allinfom[ttid,10],ties.method="random")
    if(all(allinfom[ttid,8]==3)){
      cpnv[ttid] <- 1
    } else{
      cpnv[ttid] <- 2
    }
  }
  cpnv <- as.data.frame(cpnv)
  colnames(cpnv) <- "Cp-novel-loci"
  cpod <- as.data.frame(cpod)
  colnames(cpod) <- "Cp-order"

  spod <- spnv <- allinfom[,7]
  for(i in 1:max(allinfom[,7],na.rm=T)){
    ttid <- which(allinfom[,7]==i)
    spod[ttid] <- rank(allinfom[ttid,12],ties.method="random")
    if(all(allinfom[ttid,6]==3)){
      spnv[ttid] <- 1
    } else{
      spnv[ttid] <- 2
    }
  }
  spnv <- as.data.frame(spnv)
  colnames(spnv) <- "Sp-novel-loci"
  spod <- as.data.frame(spod)
  colnames(spod) <- "Sp-order"

  write.csv(cbind(allinfom[,c(1:3,5)],cpnv,cpod,allinfom[,6:7],spnv,spod,allinfom[,c(8:10,12,17:19)]),
            file.path(cgwasenv$.CGWAS_RESULT_PATH,"Summaryhits.csv"),
            row.names=F,quote=F)
  allinfom <- allinfom[,-9]
  write.csv(allinfom,
            file.path(cgwasenv$.CGWAS_RESULT_PATH,"Allhits.csv"),
            row.names=F,quote=F)

  stopCluster(sb)
}

step11 <- function(cgwasenv) {
  # Sumbarplot

  sb <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(sb)

  coridm <- as.matrix(read.table(file.path(cgwasenv$.CGWAS_INFCOR_PATH, "EstCorrelation.txt"),
                                 header=T,stringsAsFactors=F)[,c(1:2,4:5)])
  orgcorm <- diag(cgwasenv$.TRAIT_NUM)
  n <- 0
  for(i in 2:cgwasenv$.TRAIT_NUM){
    orgcorm[i:cgwasenv$.TRAIT_NUM,i-1] <- orgcorm[i-1,i:cgwasenv$.TRAIT_NUM] <- coridm[(n+1):(n+cgwasenv$.TRAIT_NUM-i+1),3]
    n <- n+cgwasenv$.TRAIT_NUM-i+1
  }

  yl <- c(seq(0,12,2),seq(15,60,5))
  yll <- c(0,2,4,6,8,10,12,15,20,30,40,60,80,100,150,200,300)

  highcorc <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaCoefficient0"),
                         header=F,stringsAsFactors=F)
  gmqfc <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaCoefficient2"),
                      header=F,stringsAsFactors=F)

  highcorinfo <- read.table(file.path(cgwasenv$.CGWAS_GMA_PATH,"GmetaPreCombination"),
                            header=F,stringsAsFactors=F)

  highcorcnum <- (ncol(highcorc)-3)/7
  gmqfcnum <- (ncol(gmqfc)-3)/5

  dm <- read.csv(file.path(cgwasenv$.CGWAS_RESULT_PATH,"Allhits.csv"),
                 header=T,stringsAsFactors=F)
  sumdm <- read.csv(file.path(cgwasenv$.CGWAS_RESULT_PATH,"Summaryhits.csv"),
                    header=T,stringsAsFactors=F)

  allid <- sort(union(which(sumdm[,6]==1),which(sumdm[,10]==1)))
  figname <- c()
  for(ii in allid){
    if(is.na(sumdm[ii,10])){
      figname <- c(figname,paste0(sumdm[ii,1],"_",round(sumdm[ii,2]/1e6,1),"_C",sumdm[ii,4],"_",sumdm[ii,3]))
    } else if(is.na(sumdm[ii,6])){
      figname <- c(figname,paste0(sumdm[ii,1],"_",round(sumdm[ii,2]/1e6,1),"_S",sumdm[ii,8],"_",sumdm[ii,3]))
    } else if(sumdm[ii,6]==sumdm[ii,10]){
      figname <- c(figname,paste0(sumdm[ii,1],"_",round(sumdm[ii,2]/1e6,1),"_C",sumdm[ii,4],"S",sumdm[ii,8],"_",sumdm[ii,3]))
    } else if(as.numeric(sumdm[ii,6])==1){
      figname <- c(figname,paste0(sumdm[ii,1],"_",round(sumdm[ii,2]/1e6,1),"_C",sumdm[ii,4],"s",sumdm[ii,8],"_",sumdm[ii,3]))
    } else{
      figname <- c(figname,paste0(sumdm[ii,1],"_",round(sumdm[ii,2]/1e6,1),"_c",sumdm[ii,4],"S",sumdm[ii,8],"_",sumdm[ii,3]))
    }
  }

  nuu <- foreach(i=1:length(allid),.inorder=F,.export="barfig") %dopar%
    barfig(i,dm,allid,figname,yl,yll,highcorc,highcorinfo,highcorcnum,gmqfc,gmqfcnum,orgcorm,cgwasenv)

  stopCluster(sb)
}
