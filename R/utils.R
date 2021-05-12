StatE1 <- function(traitid, cgwasenv) {
  df <- data.table::fread(file.path(cgwasenv$.CGWAS_DATA_DIR, cgwasenv$.TRAIT_DATA_NAME[traitid])
                          ,header=T,stringsAsFactors=F)
  df <- as.data.frame(df)
  bpm <- df[,cgwasenv$.ASSOC_COLUMN_INDEX]
  data.table::fwrite(bpm,
                     file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid,".bp")),
                     row.names=F,col.names=T,quote=F)
  return(unique(which(is.na(bpm[,1])),which(is.na(bpm[,2]))))
}

StatE2 <- function(traitid, naid, cgwasenv) {
  bpm <- data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid,".bp")),
                           header=T,stringsAsFactors=F)
  bpm <- as.data.frame(bpm)
  if(length(naid)!=0){
    bpm <- bpm[-naid,]
  }
  bpm[,1] <- sign(bpm[,1])
  bpm[bpm[,2]<(1e-300),2] <- 1e-300
  data.table::fwrite(data.frame(signif(-qnorm(bpm[,2]/2)*bpm[,1],6)),
         file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid,".ot")),
         row.names=F,col.names=F,quote=F)
  file.remove(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid,".bp")))
  if(traitid == 1) {
    df <- data.table::fread(file.path(cgwasenv$.CGWAS_DATA_DIR, cgwasenv$.TRAIT_DATA_NAME[1])
                            ,header=T,stringsAsFactors=F)
    df <- as.data.frame(df)
    bpm <- df[,cgwasenv$.SNP_COLUMN_INDEX]
    if (length(naid) != 0) {
      bpm <- bpm[-naid,]
    }
    data.table::fwrite(bpm,
                       file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SnpIndex"),
                       sep=" ",na="NA",row.names=F,quote=F)
  }
}

ga <- function(x, cgwasenv) {
  sqrt(x)*exp(-x/2)/sqrt(2*pi)/cgwasenv$.CCP
}

qcf <- function(lam,ppv,tq,meaq,inseq) {
  qcr <- rep(NA,length(lam))
  for(i in 1:length(qcr)){
    qcr[i] <- mean(pchisq(tq/lam[i],1)*(1-ppv[i])+pchisq(tq/(lam[i]+(meaq-lam[i])/ppv[i]),1)*ppv[i]-inseq)
  }
  return(qcr^2)
}

odgridse <- function(rb,intv,lamsp,tq,meaq,N,inseq) {
  while(intv>1e-4){
    cid <- order(rb)[1]
    if(cid==1){
      lamsp <- seq(lamsp[cid],lamsp[cid+1],length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1,5)] <- rb[c(cid,cid+1)]
      rb[3] <- qcf(lamsp[3],propv[3],tq,meaq,inseq)
      intv <- intv/4
    } else if(cid==5){
      lamsp <- seq(lamsp[cid-1],lamsp[cid],length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1,5)] <- rb[c(cid-1,cid)]
      rb[3] <- qcf(lamsp[3],propv[3],tq,meaq,inseq)
      intv <- intv/4
    } else{
      lamsp <- seq(lamsp[cid-1],lamsp[cid+1],length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1,3,5)] <- rb[c(cid-1,cid,cid+1)]
      intv <- intv/2
    }
    rb[c(2,4)] <- qcf(lamsp[c(2,4)],propv[c(2,4)],tq,meaq,inseq)
  }
  return(lamsp[order(rb)[1]])
}

calinf <- function(traitid,inseq,normea,cgwasenv) {
  otm <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH,paste0(traitid,".ot")),
                                     header=F))[,1]
  tm2 <- otm^2
  meaq <- mean(tm2)
  vq <- var(tm2)
  N <- vq-2*meaq^2
  tlam <- median(tm2)/qchisq(0.5,1)
  if(N<=0){
    lamv <- min(tlam,meaq)
    propm <- 0
  } else if((meaq-1)<=0){
    lamv <- meaq-sqrt(N/1e6)
    lamv <- min(lamv,tlam)
    propm <- 1/(N/3/(meaq-lamv)^2+1)
  } else{
    tq <- quantile(tm2,inseq)
    lamsp <- seq(1,1+(meaq-1)*0.9999,length.out=5)
    propv <- 1/(N/3/(meaq-lamsp)^2+1)
    trange <- 2
    if(sum(propv>(trange*(meaq-lamsp)))!=0){
      lamsp <- seq(meaq+sqrt(1/4/trange^2-N/3)-1/2/trange,1+(meaq-1)*0.9999,length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
    }
    rb <- qcf(lamsp,propv,tq,meaq,inseq)
    intv <- (meaq-1)/2
    lamv <- odgridse(rb,intv,lamsp,tq,meaq,N,inseq)
    propm <- 1/(N/3/(meaq-lamv)^2+1)
    if((propm<1e-6)|(meaq<lamv)){
      lamv <- meaq-sqrt(N/1e6)
    }
    lamv <- min(lamv,tlam)
    propm <- 1/(N/3/(meaq-lamv)^2+1)
  }
  tm3 <- tm2 <- tm2/lamv
  backupl <- lamv
  e3 <- e2 <- mean(tm2[which(tm2>qchisq(1-cgwasenv$.CCP,1))])/normea
  pv <- 1
  while(e2<1){
    lamv <- backupl
    lamv <- lamv*(e3*pv-pv+1)
    tm2 <- tm2/(e3*pv-pv+1)
    e2 <- mean(tm2[which(tm2>qchisq(1-cgwasenv$.CCP,1))])/normea
    tm2 <- tm3
    pv <- pv+1
  }
  data.table::fwrite(as.data.frame(signif(otm/sqrt(lamv),6)),file.path(cgwasenv$.CGWAS_INFCOR_PATH,paste0(traitid,".tcomb")),row.names=F,col.names=F,quote=F)
  return(c(meaq,tlam,lamv,meaq/lamv,e2,max(0,propm)))
}

CorE <- function(i, pairma, cgwasenv) {
  t1 <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_INFCOR_PATH, paste0(pairma[i,1],".tcomb")),
                        header=F))[,1]
  t12 <- t1^2
  m12 <- mean(t12)
  if(m12<1.02){
    t1 <- t1*sqrt(1.02/m12)
    t12 <- t12*1.02/m12
  }
  v12 <- var(t1)

  t2 <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_INFCOR_PATH, paste0(pairma[i,2],".tcomb")),
                                    header=F))[,1]
  t22 <- t2^2
  m22 <- mean(t22)
  if(m22<1.02){
    t2 <- t2*sqrt(1.02/m22)
    t22 <- t22*1.02/m22
  }
  v22 <- var(t2)

  wv <- c(cgwasenv$.TRAIT_EFFECT_SIZE[pairma[i,1]],cgwasenv$.TRAIT_EFFECT_SIZE[pairma[i,2]])
  wvv <- min(wv[1]/wv[2],wv[2]/wv[1])
  tmc1 <- t12+t22
  tmc2 <- 2*t1*t2
  testid <- testcor <- 0
  selid <- tmc1<qchisq(0.5,2)
  newp <- cv <- cor(t1[selid],t2[selid])
  cm <- cov(t1,t2)
  gc <- cm/sqrt(v12*v22)
  while(abs(newp-testcor)>1e-4){
    testcor <- newp
    testid <- c(testid,testcor)
    tv <- (tmc1-testcor*tmc2)/(1-testcor^2)
    selid <- tv<qchisq(0.5,2)
    cv <- c(cv,cor(t1[selid],t2[selid]))
    rcc <- summary(lm(cv~testid))$coefficients[,1]
    newp <- rcc[1]/(1-rcc[2])
    if(abs(newp)>0.995){
      newp <- gc
    }
    if(length(testid)>5){
      testid <- testid[-1]
      cv <- cv[-1]
    }
  }
  bgc <- newp*(1+(2*abs(min(v12,v22)-1))/100)
  efc <- (cm-newp)/sqrt(abs((v12-1)*(v22-1)))*(1+0.04*wvv)
  return(c(gc,bgc,efc))
}

metaf <- function(wv,cm) {
  return(as.numeric(t(wv)%*%cm/sqrt(as.vector(t(wv)%*%cm%*%wv))))
}

copyfun.tcomb.gmeta <- function(traitid, cgwasenv) {
  file.copy(file.path(cgwasenv$.CGWAS_INFCOR_PATH, paste0(traitid, ".tcomb")),
            file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid, ".gmeta")),
            overwrite = T, recursive = T)
}

readtm.gmeta <- function(traitid, cgwasenv) {
  return(as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid,".gmeta")),
                                     header=F,stringsAsFactors=F))[,1])
}

readtm.gmetad <- function(traitid, cgwasenv) {
  return(as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid,".gmetad")),
                                     header=F,stringsAsFactors=F))[,1])
}

readtm.tcomb <- function(traitid, cgwasenv) {
  return(as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_INFCOR_PATH, paste0(traitid,".tcomb")),
                                     header=F,stringsAsFactors=F))[,1])
}

multiCorr <- function(traitid,n,newv,newv2,newv2m,orgeffs,weieffs,cgwasenv) {
  oldv <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid,".gmeta")),
                                      header=F))[,1]
  if(weieffs[traitid]<1.02){
    oldv <- oldv*sqrt(1.02/weieffs[traitid])
  }
  oldv2 <- oldv^2
  oldv2m <- var(oldv)

  # EstCorr
  wv <- orgeffs[c(traitid,n)]
  wvv <- min(wv[1]/wv[2],wv[2]/wv[1])
  tmc1 <- oldv2+newv2
  tmc2 <- 2*oldv*newv
  testid <- testcor <- 0
  selid <- tmc1<qchisq(0.5,2)
  newp <- cv <- cor(oldv[selid],newv[selid])
  cm <- as.numeric(cov(oldv,newv))
  gc <- as.numeric(cm/sqrt(oldv2m*newv2m))
  while(abs(newp-testcor)>1e-4){
    testcor <- newp
    testid <- c(testid,testcor)
    tv <- (tmc1-testcor*tmc2)/(1-testcor^2)
    selid <- tv<qchisq(0.5,2)
    cv <- c(cv,cor(oldv[selid],newv[selid]))
    rcc <- summary(lm(cv~testid))$coefficients[,1]
    newp <- rcc[1]/(1-rcc[2])
    if(abs(newp)>0.995){
      newp <- gc
    }
    if(length(testid)>5){
      testid <- testid[-1]
      cv <- cv[-1]
    }
  }
  bgc <- newp*(1+(2*abs(min(oldv2m,newv2m)-1))/100)
  efc <- (cm-newp)/sqrt(abs((oldv2m-1)*(newv2m-1)))*(1+0.04*wvv)

  resv <- c(traitid,n,signif(c(bgc,efc,gc),5))
  return(resv)
}

multiCorr.gmetad <- function(traitid,n,newv,newv2,newv2m,orgeffs,weieffs,cgwasenv) {
  oldv <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid,".gmetad")),
                                      header=F))[,1]
  if(weieffs[traitid]<1.02){
    oldv <- oldv*sqrt(1.02/weieffs[traitid])
  }
  oldv2 <- oldv^2
  oldv2m <- var(oldv)

  # EstCorr
  wv <- orgeffs[c(traitid,n)]
  wvv <- min(wv[1]/wv[2],wv[2]/wv[1])
  tmc1 <- oldv2+newv2
  tmc2 <- 2*oldv*newv
  testid <- testcor <- 0
  selid <- tmc1<qchisq(0.5,2)
  newp <- cv <- cor(oldv[selid],newv[selid])
  cm <- as.numeric(cov(oldv,newv))
  gc <- as.numeric(cm/sqrt(oldv2m*newv2m))
  while(abs(newp-testcor)>1e-4){
    testcor <- newp
    testid <- c(testid,testcor)
    tv <- (tmc1-testcor*tmc2)/(1-testcor^2)
    selid <- tv<qchisq(0.5,2)
    cv <- c(cv,cor(oldv[selid],newv[selid]))
    rcc <- summary(lm(cv~testid))$coefficients[,1]
    newp <- rcc[1]/(1-rcc[2])
    if(abs(newp)>0.995){
      newp <- gc
    }
    if(length(testid)>5){
      testid <- testid[-1]
      cv <- cv[-1]
    }
  }
  bgc <- newp*(1+(2*abs(min(oldv2m,newv2m)-1))/100)
  efc <- (cm-newp)/sqrt(abs((oldv2m-1)*(newv2m-1)))*(1+0.04*wvv)

  resv <- c(traitid,n,signif(c(bgc,efc,gc),5))
  return(resv)
}

metaf.GMA <- function(tm, sigeffcorm, trait.es.new, orgcorm, trait.num.new, cgwasenv) {
  resm <- matrix(1, nrow = nrow(tm), ncol = trait.num.new + 1)
  thv <- sqrt(qchisq(1 - cgwasenv$.CCP,1))
  for(i in 1:nrow(tm)){
    tv <- tm[i,]
    rid <- which(abs(tv)>thv)
    if(length(rid)!=0){
      temcorm <- matrix(sigeffcorm[rid,c(rid,rid + trait.num.new)], nrow = length(rid))

      # metafm
      scm <- solve(orgcorm[rid,rid])
      sigwv <- t(temcorm*trait.es.new[rid])
      coefm <- t(sigwv%*%scm/as.numeric(sqrt((sigwv%*%scm*sigwv)%*%rep(1,length(rid)))))

      pv <- pchisq((t(tv[rid]) %*% coefm)^2,1,lower.tail=F)
      resm[i,1] <- min(pv)
      resm[i,2:(trait.num.new + 1)] <- sigeffcorm[,c(rid, rid + trait.num.new)[order(pv)[1]]]
    }
  }
  return(resm)
}

qfsw <- function(tm, corm, trait.num.new, cgwasenv) {
  trtvf <- function(v,id,n,cm) {
    return(pchisq(sum(t(v[id[1:n]])%*%solve(cm[id[1:n],id[1:n]])*v[id[1:n]]),n,lower.tail=F))
  }

  resm <- matrix(NA,nrow(tm),2*trait.num.new)
  thv <- sqrt(qchisq(1-cgwasenv$.CCP,1))
  for(i in 1:nrow(tm)){
    tv <- tm[i,]
    rid <- which(abs(tv)>thv)
    if(length(rid)==0){
      resm[i,1] <- 1
    } else{
      reatop <- n <- 1
      trid <- rid <- rid[order(abs(tv[rid]),decreasing=T)]
      temtop <- trtvf(tv,rid,n,corm)
      while(temtop<reatop){
        rid <- trid
        resm[i,n] <- reatop <- temtop
        if(n==length(rid)){
          break
        }
        nh <- n <- n+1
        temtop <- trtvf(tv,rid,n,corm)
        while((temtop>=reatop)&(nh<length(rid))){
          trid <- rid
          nh <- nh+1
          t1 <- trid[nh]
          trid[(n+1):nh] <- trid[n:(nh-1)]
          trid[n] <- t1
          temtop <- trtvf(tv,trid,n,corm)
        }
      }
      resm[i,(trait.num.new + 1):(trait.num.new + length(rid))] <- rid
    }
  }
  return(resm)
}

apminfun <- function(x) {
  return(apply(x,1,min,na.rm=T))
}

copyfun.gmeta.gmetad <- function(traitid, cgwasenv) {
  file.copy(file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid, ".gmeta")),
            file.path(cgwasenv$.CGWAS_GMA_PATH, paste0(traitid, ".gmetad")),
            overwrite = T, recursive = T)
}

Simudata <- function(i, swit, orgcorm, snpn, sigeffcorm, corm2, thv, trait.es.new, trait.num.new, combnum, gmqfc, gmqfcnum) {
  mkdf <- function(cm,n) {
    df <- MASS::mvrnorm(n,mu=rep(0,ncol(cm)),Sigma=cm)
    return(df)
  }

  sortrt <- function(tv,bgc,thv) {
    rid <- which(abs(tv)>thv)
    if(length(rid)==0){
      return(1)
    } else{
      reatop <- n <- 1
      trid <- rid <- rid[order(abs(tv[rid]),decreasing=T)]
      temtop <- trtvf(tv,rid,n,bgc)
      while(temtop<reatop){
        rid <- trid
        reatop <- temtop
        if(n==length(rid)){
          break
        }
        nh <- n <- n+1
        temtop <- trtvf(tv,rid,n,bgc)
        while((temtop>=reatop)&(nh<length(rid))){
          trid <- rid
          nh <- nh+1
          t1 <- trid[nh]
          trid[(n+1):nh] <- trid[n:(nh-1)]
          trid[n] <- t1
          temtop <- trtvf(tv,trid,n,bgc)
        }
      }
      return(reatop)
    }
  }

  trtvf <- function(v,id,n,cm) {
    return(pchisq(sum(t(v[id[1:n]])%*%solve(cm[id[1:n],id[1:n]])*v[id[1:n]]),n,lower.tail=F))
  }

  metaf.simulation <- function(tv,bgc,efd,thv,trait.es.new,trait.num.new) {
    rid <- which(abs(tv)>thv)
    if(length(rid)==0){
      return(1)
    } else{
      temcorm <- matrix(efd[rid, c(rid, rid + trait.num.new)], nrow = length(rid))

      # metafm
      scm <- solve(bgc[rid,rid])
      sigwv <- t(temcorm*trait.es.new[rid])
      coefm <- t(sigwv%*%scm/as.numeric(sqrt((sigwv%*%scm*sigwv)%*%rep(1,length(rid)))))

      return(min(pchisq((t(tv[rid])%*%coefm)^2, 1, lower.tail=F)))
    }
  }

  Premeta <- function(otm, snpn, combnum, gmqfc, gmqfcnum){
    tm <- matrix(NA,snpn,combnum)
    for(i in 1:combnum){
      wid <- sum(!is.na(gmqfc[i,1:gmqfcnum]))
      tm[,i] <- as.matrix(otm[,as.numeric(gmqfc[i,1:wid])])%*%as.numeric(gmqfc[i,(1:wid)+gmqfcnum])
    }
    return(tm)
  }

  Scalem <- function(m){
    for(i in 1:ncol(m)){
      m[,i] <- scale(m[,i])
    }
    return(m)
  }

  tm <- Scalem(mkdf(orgcorm,snpn))
  pm1 <- apply(tm,1,metaf.simulation, bgc = orgcorm, efd = sigeffcorm, thv = thv, trait.es.new, trait.num.new)
  pm3 <- apply(tm,1,sortrt, bgc = orgcorm, thv = thv)
  if(swit){
    tm <- Scalem(Premeta(tm, snpn, combnum, gmqfc, gmqfcnum))
    pm2 <- apply(tm,1,sortrt, bgc = corm2, thv = thv)
    return(c(min(pm1),min(pm2),min(pm3)))
  } else{
    return(c(min(pm1),min(pm3)))
  }
}

apmaxfun <- function(x) {
  return(apply(x^2,1,max,na.rm=T))
}

fdrfun <- function(temres,thr,Gt,cgwasenv) {
  highlm <- temres[which(temres[,4]>=(-log10(thr))),]
  odtid <- order(highlm[,4],decreasing=T)
  n <- 0
  sn <- 0
  fn <- 0
  while(length(odtid)!=0){
    keyid <- odtid[1]
    if(highlm[keyid,4]>=(-log10(Gt))){
      sn <- sn+1
    } else{
      fn <- fn+1
    }
    idirc <- ddirc <- highlm[keyid,2]
    oldddirc <- oldidirc <- -1
    while((oldddirc!=ddirc)|(oldidirc!=idirc)){
      oldddirc <- ddirc
      oldidirc <- idirc
      tdl <- which((highlm[,1]==highlm[odtid[1],1])&(highlm[,2]>(ddirc-cgwasenv$.LOCISEP))&(highlm[,2]<(idirc+cgwasenv$.LOCISEP)))
      ddirc <- highlm[min(tdl),2]
      idirc <- highlm[max(tdl),2]
    }
    n <- n+1
    highlm <- highlm[-tdl,]
    odtid <- order(highlm[,4],decreasing=T)
  }
  return(c(length(which(temres[,4]>=(-log10(thr)))),length(which(temres[,4]>=(-log10(Gt)))),n,sn,fn))
}

fdrfun.summary <- function(temres, cgwasenv) {
  odtid <- order(temres[,4])
  n <- 0
  locitm <- c()
  while(length(odtid)!=0){
    keyid <- odtid[1]
    idirc <- ddirc <- temres[keyid,2]
    oldddirc <- oldidirc <- -1
    while((oldddirc!=ddirc)|(oldidirc!=idirc)){
      oldddirc <- ddirc
      oldidirc <- idirc
      tdl <- which((temres[,1]==temres[keyid,1])&(temres[,2]>(ddirc-cgwasenv$.LOCISEP))&(temres[,2]<(idirc+cgwasenv$.LOCISEP)))
      ddirc <- temres[min(tdl),2]
      idirc <- temres[max(tdl),2]
    }
    n <- n+1
    locitm <- rbind(locitm,cbind(temres[tdl,],n))
    temres <- temres[-tdl,]
    odtid <- order(temres[,4])
  }
  return(locitm[order(locitm[,4]),])
}

effvar <- function(corm) {
  newevals <- eigen(abs(corm),symmetric=T)$values
  newevals[newevals<0] <- 0
  IntLinewevals <- newevals
  IntLinewevals[IntLinewevals>=1] <- 1
  IntLinewevals[IntLinewevals<1] <- 0
  NonIntLinewevals <- newevals-floor(newevals)
  return(sum(NonIntLinewevals+IntLinewevals))
}

mandata <- function(minv,mpt,Sind) {
  orgres <- cbind(Sind,minv)
  remianid <- which(orgres[,4]>(-log10(mpt)))
  orgres <- orgres[remianid,]
  return(orgres)
}

manhattan <- function(indscp,orgminp,orgcorm,
                      thresp3m,thresorgmin,pm3min,pm3,
                      pt.col=c('gray10','gray50'),pt.bg=c('gray10','gray50'),
                      pt.cex=0.7,pch=21,cex.axis=1.1,
                      suggestiveline.col='black',suggestiveline.lwd=1,suggestiveline.lty=1
                      ,genomewideline.col='black',genomewideline.lwd=1,genomewideline.lty=2,
                      St=1e-6,Gt=5e-8,locisep=2.5e5,Lm=1,hl=3) {
  d=indscp
  colnames(d)=c("CHR","BP","SNP","P")
  d$pos=NA
  Lm <- -log10(Lm)
  deleteid <- d[,4]<Lm
  if(sum(deleteid)!=0){
    d <- d[!deleteid,]
  }
  ytop=max(ceiling(max(d$P)))
  if(ytop>20){
    d$P[d$P>20]=d$P[d$P>20]/2+10
    if(ytop>40){
      d$P[d$P>30]=d$P[d$P>30]/2+15
      if(ytop>100){
        d$P[d$P>45]=d$P[d$P>45]/2.5+27
        if(ytop>200){
          d$P[d$P>55]=d$P[d$P>55]/2+27.5
          if(ytop>300){
            d$P[d$P>60]=60
          }
        }
      }
    }
  }
  ymax=ceiling(max(d$P))
  ymin=Lm

  d$index=NA
  ind=0
  for(i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }

  nchr=length(unique(d$CHR))
  if (nchr==1){
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel=paste('Chromosome',unique(d$CHR),'position')
    labs=ticks
  } else{
    ticks=rep(NA,length(unique(d$CHR))+1)
    ticks[1]=0
    for(i in 1:max(d$index)){
      d[d$index==i,]$pos=(d[d$index==i,]$BP-d[d$index==i,]$BP[1])+1+ticks[i]
      ticks[i+1]=max(d[d$index==i,]$pos)+2e7
    }
    xlabel = 'Chromosome'
    labs = unique(d$CHR)
  }

  xmax=max(d$pos)*1.03
  xmin=max(d$pos)*-0.03
  plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab=expression(-log[10](italic(p))),las=1,cex.axis=cex.axis)

  blank=rep('',length(labs))
  lowerlabs=rep('',length(labs))
  upperlabs=rep('',length(labs))

  for (i in 1:length(labs)){
    if (i %% 2 == 0){
      lowerlabs[i] = labs[i]
    } else{
      upperlabs[i] = labs[i]
    }
  }

  newtick <- c()
  for(i in 1:(length(ticks)-1)){
    newtick <- c(newtick,(ticks[i]+ticks[i+1])/2)
  }

  l=c(0,c(seq(2,12,2),seq(15,60,5)))
  ll=c(0,c(2,4,6,8,10,12,15,20,30,40,60,80,100,150,200,300))
  ll=ll[(l>=Lm)&(l<=floor(ymax))]
  l=l[(l>=Lm)&(l<=floor(ymax))]

  axis(2,las=1,at=l,labels=as.character(ll),lwd=0,lwd.ticks=1,cex.axis=cex.axis)

  d$P[which(d$P==0)] <- NA
  pt.col = rep(pt.col,max(d$CHR))[1:max(d$CHR)]
  pt.bg = rep(pt.bg,max(d$CHR))[1:max(d$CHR)]

  icol=1
  for (i in unique(d$CHR)) {
    with(d[d$CHR==i, ],points(pos, P, col=pt.col[icol],bg=pt.bg[icol],cex=pt.cex,pch=pch))
    icol=icol+1
  }

  box()

  if(hl==1){
    tid <- which(d$P>=(-log10(St)))
    highlm <- d[tid,]
    temp3m <- as.data.frame(pm3[as.numeric(rownames(highlm)),])
    temorgp <- (orgminp*effvar(orgcorm))[as.numeric(rownames(highlm))]
    orgSt <- thresorgmin*effvar(orgcorm)
    odtid <- order(highlm$P,decreasing=T)
    n <- 0
    nv <- 0
    sn <- 0
    fn <- 0
    while(length(odtid)!=0){
      keyid <- odtid[1]
      if(highlm[keyid,4]>=(-log10(Gt))){
        sn <- sn+1
      } else{
        fn <- fn+1
      }
      idirc <- ddirc <- highlm$pos[keyid]
      oldddirc <- oldidirc <- -1
      while((oldddirc!=ddirc)|(oldidirc!=idirc)){
        oldddirc <- ddirc
        oldidirc <- idirc
        tdl <- which((highlm$CHR==highlm$CHR[odtid[1]])&(highlm$pos>(ddirc-locisep))&(highlm$pos<(idirc+locisep)))
        ddirc <- highlm$pos[min(tdl)]
        idirc <- highlm$pos[max(tdl)]
      }
      tpv <- as.numeric(c(min(temp3m[tdl,1],na.rm = T),min(temp3m[tdl,2],na.rm = T),min(temp3m[tdl,3],na.rm = T)))
      colv <- c("brown2","darkgoldenrod2","limegreen")
      colv <- colv[tpv<=St]
      tpv <- tpv[tpv<=St]
      if(min(temorgp[tdl])<=orgSt){
        tpv <- c(tpv,min(temorgp[tdl]))
        colv <- c(colv,"mediumpurple3")
      } else{
        nv <- nv+1
      }
      colv <- colv[order(tpv)]
      tpv <- tpv[order(tpv)]

      for(ni in 1:length(tpv)){
        points(highlm$pos[keyid],highlm$P[keyid]+(ymax-Lm)*(0.008+length(tpv)*0.008),col=colv[ni],cex=0.9*length(tpv)+0.6-0.9*ni,lwd=2.5,pch=5,xpd=T)
      }

      n <- n+1
      highlm <- highlm[-tdl,]
      temp3m <- temp3m[-tdl,]
      temorgp <- temorgp[-tdl]
      odtid <- order(highlm$P,decreasing=T)
    }

    text(max(d$pos)*0.92,ymax,expression(paste(N[(Loci)]," = ")),adj=1,cex=1.5)
    text(max(d$pos)*0.92,ymax,paste0(n,"(",nv,")"),adj=0,cex=1.5)
    text(max(d$pos)*0.92,ymax*0.96,expression(paste(N[(Sig)],"/",N[(Fdr)]," : ")),adj=1,cex=1.5)
    text(max(d$pos)*0.92,ymax*0.96,paste0(sn,"/",fn),adj=0,cex=1.5)
  } else if(hl==2){
    tid <- which(d$P>=(-log10(St)))
    highlm <- d[tid,]
    pm3lm <- pm3min[as.numeric(rownames(highlm))]
    #    	pm3lm <- orgminp[as.numeric(rownames(highlm))]
    odtid <- order(highlm$P,decreasing=T)
    n <- 0
    nv <- 0
    sn <- 0
    fn <- 0
    while(length(odtid)!=0){
      keyid <- odtid[1]
      if(highlm[keyid,4]>=(-log10(Gt))){
        sn <- sn+1
      } else{
        fn <- fn+1
      }
      idirc <- ddirc <- highlm$pos[keyid]
      oldddirc <- oldidirc <- -1
      while((oldddirc!=ddirc)|(oldidirc!=idirc)){
        oldddirc <- ddirc
        oldidirc <- idirc
        tdl <- which((highlm$CHR==highlm$CHR[odtid[1]])&(highlm$pos>(ddirc-locisep))&(highlm$pos<(idirc+locisep)))
        ddirc <- highlm$pos[min(tdl)]
        idirc <- highlm$pos[max(tdl)]
      }

      points(highlm$pos[keyid],highlm$P[keyid]+(ymax-Lm)*0.02, col="mediumpurple3",cex=1,lwd=2.5,pch=5)
      if(min(pm3lm[tdl])>thresp3m){
        #			if(min(pm3lm[tdl])>thresorgmin){
        points(highlm$pos[keyid],highlm$P[keyid]+(ymax-Lm)*0.02, col="darkgoldenrod2",cex=1,lwd=2,pch=4)
        nv <- nv+1
      }

      pm3lm <- pm3lm[-tdl]
      highlm <- highlm[-tdl,]
      odtid <- order(highlm$P,decreasing=T)
      n <- n+1
    }

    text(max(d$pos)*0.92,ymax,expression(paste(N[(Loci)]," = ")),adj=1,cex=1.5)
    text(max(d$pos)*0.92,ymax,paste0(n,"(",nv,")"),adj=0,cex=1.5)
    text(max(d$pos)*0.92,ymax*0.96,expression(paste(N[(Sig)],"/",N[(Fdr)]," : ")),adj=1,cex=1.5)
    text(max(d$pos)*0.92,ymax*0.96,paste0(sn,"/",fn),adj=0,cex=1.5)
  } else if(hl==3){
    tid <- which(d$P>=(-log10(St)))
    highlm <- d[tid,]
    odtid <- order(highlm$P,decreasing=T)
    n <- 0
    sn <- 0
    fn <- 0
    while(length(odtid)!=0){
      keyid <- odtid[1]
      if(highlm[keyid,4]>=(-log10(Gt))){
        sn <- sn+1
      } else{
        fn <- fn+1
      }
      idirc <- ddirc <- highlm$pos[keyid]
      points(highlm$pos[keyid],highlm$P[keyid]+(ymax-Lm)*0.03, col="mediumpurple3",cex=1,lwd=2.5,pch=5)
      oldddirc <- oldidirc <- -1
      while((oldddirc!=ddirc)|(oldidirc!=idirc)){
        oldddirc <- ddirc
        oldidirc <- idirc
        tdl <- which((highlm$CHR==highlm$CHR[odtid[1]])&(highlm$pos>(ddirc-locisep))&(highlm$pos<(idirc+locisep)))
        ddirc <- highlm$pos[min(tdl)]
        idirc <- highlm$pos[max(tdl)]
      }
      highlm <- highlm[-tdl,]
      odtid <- order(highlm$P,decreasing=T)
      n <- n+1
    }
    text(max(d$pos)*0.92,ymax,expression(paste(N[(Loci)]," = ")),adj=1,cex=1.2)
    text(max(d$pos)*0.92,ymax,n,adj=0,cex=1.2)
    text(max(d$pos)*0.92,ymax*0.96,expression(paste(N[(Sig)],"/",N[(Fdr)]," : ")),adj=1,cex=1.2)
    text(max(d$pos)*0.92,ymax*0.96,paste0(sn,"/",fn),adj=0,cex=1.2)
  } else if(hl==4){
    tid <- which(d$P>=(-log10(St)))
    highlm <- d[tid,]
    temp3m <- as.data.frame(pm3[as.numeric(rownames(highlm)),])
    temorgp <- (orgminp*effvar(orgcorm))[as.numeric(rownames(highlm))]
    orgSt <- thresorgmin*effvar(orgcorm)
    odtid <- order(highlm$P,decreasing=T)
    n <- 0
    sn <- 0
    fn <- 0
    while(length(odtid)!=0){
      sni <- 0
      fni <- 0
      keyid <- odtid[1]
      if(highlm[keyid,4]>=(-log10(Gt))){
        sni <- 1
      } else{
        fni <- 1
      }
      idirc <- ddirc <- highlm$pos[keyid]
      oldddirc <- oldidirc <- -1
      while((oldddirc!=ddirc)|(oldidirc!=idirc)){
        oldddirc <- ddirc
        oldidirc <- idirc
        tdl <- which((highlm$CHR==highlm$CHR[odtid[1]])&(highlm$pos>(ddirc-locisep))&(highlm$pos<(idirc+locisep)))
        ddirc <- highlm$pos[min(tdl)]
        idirc <- highlm$pos[max(tdl)]
      }
      tpv <- as.numeric(c(min(temp3m[tdl,1]),min(temp3m[tdl,2]),min(temp3m[tdl,3])))
      colv <- c("brown2","darkgoldenrod2","limegreen")
      colv <- colv[tpv<=St]
      tpv <- tpv[tpv<=St]
      if(min(temorgp[tdl])>orgSt){
        colv <- colv[order(tpv)]
        tpv <- tpv[order(tpv)]
        for(ni in 1:length(tpv)){
          points(highlm$pos[keyid],highlm$P[keyid]+(ymax-Lm)*(0.008+length(tpv)*0.008),col=colv[ni],cex=0.75*length(tpv)+0.5-0.75*ni,lwd=2,pch=5,xpd=T)
        }
        n <- n+1
        if(sni==1){
          sn <- sn+1
        }
        if(fni==1){
          fn <- fn+1
        }
      }

      highlm <- highlm[-tdl,]
      temp3m <- temp3m[-tdl,]
      temorgp <- temorgp[-tdl]
      odtid <- order(highlm$P,decreasing=T)
    }
    text(max(d$pos)*0.92,ymax,expression(paste(N[(Loci)]," = ")),adj=1,cex=1.2)
    text(max(d$pos)*0.92,ymax,n,adj=0,cex=1.2)
    text(max(d$pos)*0.92,ymax*0.96,expression(paste(N[(Sig)],"/",N[(Fdr)]," : ")),adj=1,cex=1.2)
    text(max(d$pos)*0.92,ymax*0.96,paste0(sn,"/",fn),adj=0,cex=1.2)
  }

  text(newtick,ymin-(ymax-ymin)*0.06,lowerlabs,cex=0.8,xpd=NA)
  text(newtick,ymin-(ymax-ymin)*0.06,upperlabs,cex=0.8,xpd=NA)

  abline(h=-log10(St), col=suggestiveline.col[1],lwd=suggestiveline.lwd,lty=suggestiveline.lty)
  abline(h=-log10(Gt), col=genomewideline.col[1],lwd=genomewideline.lwd,lty=genomewideline.lty)
}

manbar <- function(d,temp3m,orgminpp,orgtm2,orgcorm,thresorgmin,orgminp,cex.axis=1.1,St=1e-6,Gt=5e-8,Lm=1,locisep=2.5e5) {
  colnames(d)=c("CHR","BP","SNP","P")
  Lm <- -log10(Lm)
  temp3mp <- temp3m[as.numeric(rownames(d)),]
  temp3m <- -log10(temp3m)[as.numeric(rownames(d)),]
  orgminprop <- -log10(orgminpp*effvar(orgcorm))[as.numeric(rownames(d))]
  orgminpp <- (orgminpp*effvar(orgcorm))[as.numeric(rownames(d))]
  orgtm2 <- matrix(orgtm2[as.numeric(rownames(d)),], nrow = nrow(d))
  orgtm2 <- pchisq(orgtm2^2,1,lower.tail=F)
  deleteid <- d[,4]<Lm
  if(sum(deleteid)!=0){
    d <- d[!deleteid,]
    orgminpp <- orgminpp[!deleteid]
    temp3m <- temp3m[!deleteid,]
    orgminprop <- orgminprop[!deleteid]
    orgtm2 <- orgtm2[!deleteid,]
  }
  ytop=max(ceiling(max(d$P)))
  if(ytop>20){
    d$P[d$P>20]=d$P[d$P>20]/2+10
    temp3m[temp3m>20]=temp3m[temp3m>20]/2+10
    orgminprop[orgminprop>20]=orgminprop[orgminprop>20]/2+10
    if(ytop>40){
      d$P[d$P>30]=d$P[d$P>30]/2+15
      temp3m[temp3m>30]=temp3m[temp3m>30]/2+15
      orgminprop[orgminprop>30]=orgminprop[orgminprop>30]/2+15
      if(ytop>100){
        d$P[d$P>45]=d$P[d$P>45]/2.5+27
        temp3m[temp3m>45]=temp3m[temp3m>45]/2.5+27
        orgminprop[orgminprop>45]=orgminprop[orgminprop>45]/2.5+27
        if(ytop>200){
          d$P[d$P>55]=d$P[d$P>55]/2+27.5
          temp3m[temp3m>55]=temp3m[temp3m>55]/2+27.5
          orgminprop[orgminprop>55]=orgminprop[orgminprop>55]/2+27.5
          if(ytop>300){
            d$P[d$P>60]=60
            temp3m[temp3m>60]=60
            orgminprop[orgminprop>60]=60
          }
        }
      }
    }
  }
  ymax=ceiling(max(d$P))

  tid <- which(d$P>=(-log10(St)))
  highlm <- d[tid,]
  temtemp3m <- matrix(temp3m[tid,], nrow = length(tid))
  orgtm2p <- matrix(orgtm2[tid,], nrow = length(tid))
  temp3morgp <- matrix(temp3mp[tid,], nrow = length(tid))
  temorgp <- orgminprop[tid]
  temorgprop <- as.numeric(apply(temp3morgp,1,min))/orgminpp[tid]
  orgSt <- thresorgmin*effvar(orgcorm)
  odtid <- order(highlm$P,decreasing=T)

  cgcp <- cgcl <- cgsp <- cgprop <- cghl <- cgty <- c()
  while(length(odtid)!=0){
    keyid <- odtid[1]
    idirc <- ddirc <- highlm$BP[keyid]
    oldddirc <- oldidirc <- -1
    while((oldddirc!=ddirc)|(oldidirc!=idirc)){
      oldddirc <- ddirc
      oldidirc <- idirc
      tdl <- which((highlm$CHR==highlm$CHR[odtid[1]])&(highlm$BP>(ddirc-locisep))&(highlm$BP<(idirc+locisep)))
      ddirc <- highlm$BP[min(tdl)]
      idirc <- highlm$BP[max(tdl)]
    }
    tpv <- as.numeric(c(max(temtemp3m[tdl,1],na.rm = T),max(temtemp3m[tdl,2],na.rm = T),max(temtemp3m[tdl,3],na.rm = T),max(temorgp[tdl])))
    colv <- c("brown2","darkgoldenrod2","limegreen","mediumpurple3")
    cgcp <- rbind(cgcp,tpv[order(tpv,decreasing=T)])
    cgcl <- rbind(cgcl,colv[order(tpv,decreasing=T)])
    cgsp <- rbind(cgsp,orgtm2p[keyid,])
    cgprop <- c(cgprop,temorgprop[keyid])
    cgty <- c(cgty,1)
    if(tpv[4]>=(-log10(orgSt))){
      cghl <- c(cghl,2)
    } else{
      cghl <- c(cghl,1)
    }

    highlm <- highlm[-tdl,]
    temtemp3m <- matrix(temtemp3m[-tdl,], nrow = nrow(temtemp3m) - length(tdl))
    temorgp <- temorgp[-tdl]
    temorgprop <- temorgprop[-tdl]
    orgtm2p <-  matrix(orgtm2p[-tdl,], nrow = nrow(orgtm2p) - length(tdl))
    odtid <- order(highlm$P,decreasing=T)
  }

  tid <- which(orgminp[as.numeric(rownames(d))]<=thresorgmin)
  highlm <- d[tid,]
  temtemp3m <- as.data.frame(temp3m[tid,])
  orgtm2p <- as.data.frame(orgtm2[tid,])
  temp3morgp <- as.data.frame(temp3mp[tid,])
  temorgp <- orgminprop[tid]
  temorgprop <- as.numeric(apply(temp3morgp,1,min))/orgminpp[tid]
  odtid <- order(highlm$P,decreasing=T)

  while(length(odtid)!=0){
    keyid <- odtid[1]
    idirc <- ddirc <- highlm$BP[keyid]
    oldddirc <- oldidirc <- -1
    while((oldddirc!=ddirc)|(oldidirc!=idirc)){
      oldddirc <- ddirc
      oldidirc <- idirc
      tdl <- which((highlm$CHR==highlm$CHR[odtid[1]])&(highlm$BP>(ddirc-locisep))&(highlm$BP<(idirc+locisep)))
      ddirc <- highlm$BP[min(tdl)]
      idirc <- highlm$BP[max(tdl)]
    }
    if(max(highlm$P[tdl])<(-log10(St))){
      tpv <- as.numeric(c(max(temtemp3m[tdl,1],na.rm = T),max(temtemp3m[tdl,2],na.rm = T),max(temtemp3m[tdl,3],na.rm = T),max(temorgp[tdl])))
      colv <- c("brown2","darkgoldenrod2","limegreen","mediumpurple3")
      cgcp <- rbind(cgcp,tpv[order(tpv,decreasing=T)])
      cgcl <- rbind(cgcl,colv[order(tpv,decreasing=T)])
      cgsp <- rbind(cgsp,orgtm2p[keyid,])
      cgprop <- c(cgprop,temorgprop[keyid])
      cgty <- c(cgty,2)
      cghl <- c(cghl,0)
    }
    highlm <- highlm[-tdl,]
    temtemp3m <- temtemp3m[-tdl,]
    temorgp <- temorgp[-tdl]
    temorgprop <- temorgprop[-tdl]
    orgtm2p <- orgtm2p[-tdl,]
    odtid <- order(highlm$P,decreasing=T)
  }

  cgcp <- matrix(cgcp[order(cgprop),], nrow = length(cgprop))
  cgcl <- matrix(cgcl[order(cgprop),], nrow = length(cgprop))
  cgsp <- as.matrix(cgsp[order(cgprop),])
  cghl <- cghl[order(cgprop)]
  cgprop <- cgprop[order(cgprop)]
  if(length(cgty)>100){
    if((sum(cgty==1)>50)&(sum(cgty==2)>50)){
      cgcp <- cgcp[-(51:(nrow(cgcp)-50)),]
      cgcl <- cgcl[-(51:(nrow(cgcl)-50)),]
      cgsp <- cgsp[-(51:(nrow(cgsp)-50)),]
      cghl <- cghl[-(51:(length(cghl)-50))]
      cgprop <- cgprop[-(51:(length(cgprop)-50))]
    } else if(sum(cgty==1)>50){
      cgcp <- cgcp[-((100-sum(cgty==2)+1):(nrow(cgcp)-sum(cgty==2))),]
      cgcl <- cgcl[-((100-sum(cgty==2)+1):(nrow(cgcl)-sum(cgty==2))),]
      cgsp <- cgsp[-((100-sum(cgty==2)+1):(nrow(cgsp)-sum(cgty==2))),]
      cghl <- cghl[-((100-sum(cgty==2)+1):(length(cghl)-sum(cgty==2)))]
      cgprop <- cgprop[-((100-sum(cgty==2)+1):(length(cgprop)-sum(cgty==2)))]
    } else if(sum(cgty==2)>50){
      cgcp <- cgcp[-((sum(cgty==1)+1):(nrow(cgcp)-100+sum(cgty==1))),]
      cgcl <- cgcl[-((sum(cgty==1)+1):(nrow(cgcl)-100+sum(cgty==1))),]
      cgsp <- cgsp[-((sum(cgty==1)+1):(nrow(cgsp)-100+sum(cgty==1))),]
      cghl <- cghl[-((sum(cgty==1)+1):(length(cghl)-100+sum(cgty==1)))]
      cgprop <- cgprop[-((sum(cgty==1)+1):(length(cgprop)-100+sum(cgty==1)))]
    }
  }

  cgcp[cgcp<Lm] <- Lm
  temm <- matrix(0,nrow(cgsp),9)
  for(i in 1:nrow(temm)){
    for(ii in 1:9){
      temm[i,ii] <- sum(cgsp[i,]<(1/(10^ii)))
    }
  }
  cgsp <- temm
  ym <- max(cgsp)
  cgsp <- Lm-3*ymax/25-cgsp/ym*ymax/3

  if(sum(cgprop>1)!=0){
    cgl <- (nrow(cgcp)+2)
    if (sum(cgprop<1) >= 1) {
      cgl <- c(1:sum(cgprop<1),(sum(cgprop<1)+3):cgl)
    }
    xmax=max(cgl)+2
  } else{
    cgl <- 1:nrow(cgcp)
    xmax=max(cgl)+4
  }

  ymin <- Lm-ymax/3-3*ymax/25
  xmin=0
  plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab=expression(paste("      ",N[GWAS],"                                                                ",-log[10](italic(p)),"                                ")),las=1,cex.axis=cex.axis)

  l=c(0,c(seq(2,12,2),seq(15,60,5)))
  ll=c(0,c(2,4,6,8,10,12,15,20,30,40,60,80,100,150,200,300))
  ll=ll[(l>=Lm)&(l<=floor(ymax))]
  l=l[(l>=Lm)&(l<=floor(ymax))]

  axis(2,las=1,at=l,labels=as.character(ll),lwd=0,lwd.ticks=1,cex.axis=cex.axis)
  axis(2,las=1,at=c(Lm-3*ymax/25-ymax/3*round(ym*0.2)/ym,Lm-3*ymax/25-ymax/3*round(ym*0.4)/ym,Lm-3*ymax/25-ymax/3*round(ym*0.6)/ym,Lm-3*ymax/25-ymax/3*round(ym*0.8)/ym,Lm-3*ymax/25-ymax/3),labels=as.character(c(round(ym*0.2),round(ym*0.4),round(ym*0.6),round(ym*0.8),ym)),lwd=0,lwd.ticks=1,cex.axis=cex.axis)

  colp <- c("grey90","grey80","grey70","grey60","grey50","grey40","grey30","grey20","grey10")

  for(i in 1:4){
    did <- which(cgcl[,i]=="mediumpurple3")
    if(length(did)!=0){
      rect(cgl[did],Lm,cgl[did]+1,cgcp[did,i],col="mediumpurple3",border="mediumpurple4",lwd=1.5)
      oid <- which(cgcl[,i]!="mediumpurple3")
      if (length(oid) > 0) {
        rect(cgl[oid],Lm,cgl[oid]+1,cgcp[oid,i],col=cgcl[oid,i],border=T,density=20,lwd=1.5)
      }
    } else{
      rect(cgl,Lm,cgl+1,cgcp[,i],col=cgcl[,i],border=T,density=20,lwd=1.5)
    }
  }

  for(i in 1:9){
    rect(cgl,cgsp[,i],cgl+1,Lm-3*ymax/25,col=colp[i],border=NA)
  }

  if(sum(cghl==2)!=0){
    did <- which(cghl==2)
    rect(cgl[did],Lm-2*ymax/25,cgl[did]+1,Lm-ymax/25,col="darkgoldenrod2",border=NA)
  }
  if(sum(cghl==1)!=0){
    did <- which(cghl==1)
    rect(cgl[did],Lm-2*ymax/25,cgl[did]+1,Lm-ymax/25,col="limegreen",border=NA)
  }

  abline(v=(sum(cgprop<1)+2))

  box()
}

barfig <- function(ii,dm,allid,figname,yl,yll,highcorc,highcorinfo,highcorcnum,gmqfc,gmqfcnum,orgcorm,cgwasenv) {

  logtran <- function(v){
    v <- -log10(v)
    ytop=max(ceiling(max(v)))
    if(ytop>20){
      v[v>20]=v[v>20]/2+10
      if(ytop>40){
        v[v>30]=v[v>30]/2+15
        if(ytop>100){
          v[v>45]=v[v>45]/2.5+27
          if(ytop>200){
            v[v>55]=v[v>55]/2+27.5
            if(ytop>300){
              v[v>60]=60
            }
          }
        }
      }
    }
    return(v)
  }

  effvar <- function(corm){
    newevals <- eigen(abs(corm),symmetric=T)$values
    newevals[newevals<0] <- 0
    IntLinewevals <- newevals
    IntLinewevals[IntLinewevals>=1] <- 1
    IntLinewevals[IntLinewevals<1] <- 0
    NonIntLinewevals <- newevals-floor(newevals)
    return(sum(NonIntLinewevals+IntLinewevals))
  }

  pdm <- dm[allid[ii],]

  jpeg(file.path(cgwasenv$.CGWAS_BARPLOTS_PATH,paste0(figname[ii],".jpg")),width=(4000+cgwasenv$.TRAIT_NUM*50),height=5000,res=600)
  par(mar=c(3.5,4,2,1))
  layout(matrix(c(4,1,3,2),2,2))
  labcex <- 0.4

  pdm2 <- as.numeric(logtran(pdm[c(16,13)]))
  tdm3 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+1:nrow(highcorinfo))])
  pdm3 <- as.numeric(logtran(pdm[c(18+2*cgwasenv$.TRAIT_NUM+nrow(highcorinfo)+1:nrow(highcorinfo))]))
  ddm3 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+2*nrow(highcorinfo)+1:nrow(highcorinfo))])

  oid <- which(pdm3<1)
  pid <- which((pdm3>=1)&(ddm3==1))
  nid <- which((pdm3>=1)&(ddm3==-1))

  ymax <- max(c(pdm2,pdm3))
  spb <- ymax/50

  plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(0,nrow(highcorinfo)+as.numeric(length(oid)!=0)+as.numeric(length(pid)!=0)+as.numeric(length(nid)!=0)+4),ylim=c(0,ymax),xlab="",ylab=expression(-log[10](italic(p))),main="Generalized Meta",las=1,cex.axis=0.7,mgp=c(2.25,1,0))
  axis(2,las=1,at=yl[yl<=floor(ymax)],labels=as.character(yll[yl<=floor(ymax)]),lwd.ticks=1,cex.axis=0.7)
  rect(1,0,2,pdm2[1],col="brown2",border="brown4",lwd=1.25)
  rect(2,0,3,pdm2[2],col="brown4",border=T,density=25,lwd=1.25)
  text(1:2+0.5,rep(-spb,4),c("G-Meta","G-Meta-uncorr"),cex=labcex,srt=45,adj=1,xpd=NA)

  n <- 4
  di <- 0
  if(length(pid)!=0){
    di <- di+sum(tdm3[pid])
  }
  if(length(nid)!=0){
    di <- di-sum(tdm3[nid])
  }

  gmeft <- c()
  if(length(pid)!=0){
    opid <- pid[order(pdm3[pid]*tdm3[pid]/abs(tdm3[pid]))]
    tdm3b <- tdm3c <- tdm3[opid]
    pd <- which(tdm3c>0)
    nd <- which(tdm3c<0)
    tdm3c[pd] <- "salmon1"
    tdm3b[pd] <- "chocolate3"
    tdm3c[nd] <- "skyblue1"
    tdm3b[nd] <- "steelblue3"
    rect(n:(n+length(pid)-1),0,(n+1):(n+length(pid)),pdm3[opid],col=tdm3c,border=tdm3b)
    if(di>0){
      if(length(nd)!=0){
        text(n+nd-0.5,pdm3[opid[nd]]+spb,"*",cex=1.2,col="grey20",xpd=NA)
      }
      gmeft <- c(gmeft,opid[pd])
    } else{
      if(length(pd)!=0){
        text(n+pd-0.5,pdm3[opid[pd]]+spb,"*",cex=1.2,col="grey20",xpd=NA)
      }
      gmeft <- c(gmeft,opid[nd])
    }
    text(n:(n+length(pid)-1)+0.5,rep(-spb,length(pid)),highcorinfo[opid,3],cex=labcex,col="chocolate3",srt=45,adj=1,xpd=NA)
    n <- n+length(pid)+1
  }
  if(length(nid)!=0){
    onid <- nid[order(pdm3[nid]*tdm3[nid]/abs(tdm3[nid]))]
    tdm3b <- tdm3c <- tdm3[onid]
    pd <- which(tdm3c>0)
    nd <- which(tdm3c<0)
    tdm3c[pd] <- "salmon1"
    tdm3b[pd] <- "chocolate3"
    tdm3c[nd] <- "skyblue1"
    tdm3b[nd] <- "steelblue3"
    rect(n:(n+length(nid)-1),0,(n+1):(n+length(nid)),pdm3[onid],col=tdm3c,border=tdm3b,lwd=1.25)
    if(di>0){
      if(length(pd)!=0){
        text(n+pd-0.5,pdm3[onid[pd]]+spb,"*",cex=1.2,col="grey20",xpd=NA)
      }
      gmeft <- c(gmeft,onid[nd])
    } else{
      if(length(nd)!=0){
        text(n+nd-0.5,pdm3[onid[nd]]+spb,"*",cex=1.2,col="grey20",xpd=NA)
      }
      gmeft <- c(gmeft,onid[pd])
    }
    text(n:(n+length(nid)-1)+0.5,rep(-spb,length(nid)),highcorinfo[onid,3],cex=labcex,col="steelblue3",srt=45,adj=1,xpd=NA)
    n <- n+length(nid)+1
  }

  if(length(oid)!=0){
    ooid <- oid[order(pdm3[oid],decreasing=T)]
    rect(n:(n+length(oid)-1),0,(n+1):(n+length(oid)),pdm3[ooid],col="grey80",border="grey20",lwd=1.25)
    text(n:(n+length(oid)-1)+0.5,rep(-spb,length(oid)),highcorinfo[ooid,3],cex=labcex,col="grey20",srt=45,adj=1,xpd=NA)
  }


  pdm2 <- as.numeric(logtran(pdm[c(18,15)]))
  odm3 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+3*nrow(highcorinfo)+1:nrow(highcorinfo))])
  odm3 <- odm3[!is.na(odm3)]
  opdm3 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+4*nrow(highcorinfo)+1:nrow(highcorinfo))])
  opdm3 <- logtran(opdm3[!is.na(opdm3)])

  ymax <- max(c(pdm2,pdm3,opdm3))
  spb <- ymax/50

  plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(0,nrow(highcorinfo)+as.numeric(length(odm3)!=nrow(highcorinfo))+5),ylim=c(0,ymax),xlab="",ylab=expression(-log[10](italic(p))),main="Quadratic function",las=1,cex.axis=0.7,mgp=c(2.25,1,0))
  axis(2,las=1,at=yl[yl<=floor(ymax)],labels=as.character(yll[yl<=floor(ymax)]),lwd.ticks=1,cex.axis=0.7)
  rect(1,0,2,pdm2[1],col="chartreuse3",border="green4",lwd=1.25)
  rect(2,0,3,pdm2[2],col="green4",border=T,density=25,lwd=1.25)
  text(1:2+0.5,rep(-spb,4),c("Sw-Qf","Sw-Qf-uncorr"),cex=labcex,srt=45,adj=1,xpd=NA)

  n <- 4
  rect(n:(n+length(opdm3)-1),0,(n+1):(n+length(opdm3)),opdm3,col="chartreuse3",border="green4",lwd=1.25)
  tdm3b <- tdm3c <- tdm3[odm3]
  pd <- which(tdm3c>0)
  nd <- which(tdm3c<0)
  tdm3c[pd] <- "salmon1"
  tdm3b[pd] <- "chocolate3"
  tdm3c[nd] <- "skyblue1"
  tdm3b[nd] <- "steelblue3"
  rect(n:(n+length(odm3)-1),0,(n+1):(n+length(odm3)),pdm3[odm3],col=tdm3c,border=tdm3b,lwd=1.25)
  text(n:(n+length(odm3)-1)+0.5,rep(-spb,length(odm3)),highcorinfo[odm3,3],cex=labcex,col="grey20",srt=45,adj=1,xpd=NA)
  text((n:(n+length(odm3)-1))[1:length(opdm3)]+0.5,rep(-spb,length(opdm3)),highcorinfo[odm3[1:length(opdm3)],3],cex=labcex,col="green4",srt=45,adj=1,xpd=NA)
  qfeft <- odm3[1:length(opdm3)]
  n <- n+length(odm3)+1

  if(length(setdiff(1:nrow(highcorinfo),odm3))!=0){
    oid <- setdiff(1:nrow(highcorinfo),odm3)
    ooid <- oid[order(pdm3[oid],decreasing=T)]
    rect(n:(n+length(oid)-1),0,(n+1):(n+length(oid)),pdm3[ooid],col="grey80",border="grey20",lwd=1.25)
    text(n:(n+length(oid)-1)+0.5,rep(-spb,length(oid)),highcorinfo[ooid,3],cex=labcex,col="grey20",srt=45,adj=1,xpd=NA)
  }


  pdm2 <- as.numeric(logtran(pdm[c(17,14)]))
  tdm5 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+5*nrow(highcorinfo)+1:nrow(gmqfc))])
  pdm5 <- as.numeric(logtran(pdm[c(18+2*cgwasenv$.TRAIT_NUM+5*nrow(highcorinfo)+nrow(gmqfc)+1:nrow(gmqfc))]))
  odm4 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+5*nrow(highcorinfo)+2*nrow(gmqfc)+1:nrow(gmqfc))])
  odm4 <- odm4[!is.na(odm4)]
  opdm4 <- as.numeric(pdm[c(18+2*cgwasenv$.TRAIT_NUM+5*nrow(highcorinfo)+3*nrow(gmqfc)+1:nrow(gmqfc))])
  opdm4 <- logtran(opdm4[!is.na(opdm4)])

  ymax <- max(c(pdm2,pdm3,pdm5,opdm4))
  spb <- ymax/50

  plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(0,nrow(highcorinfo)+as.numeric(length(odm4)!=nrow(gmqfc))+5),ylim=c(0,ymax),xlab="",ylab=expression(-log[10](italic(p))),main="Generalized Meta + Quadratic function",las=1,cex.axis=0.7,mgp=c(2.25,1,0))
  axis(2,las=1,at=yl[yl<=floor(ymax)],labels=as.character(yll[yl<=floor(ymax)]),lwd.ticks=1,cex.axis=0.7)
  rect(1,0,2,pdm2[1],col="darkgoldenrod2",border="darkgoldenrod4",lwd=1.25)
  rect(2,0,3,pdm2[2],col="darkgoldenrod4",border=T,density=25,lwd=1.25)
  text(1:2+0.5,rep(-spb,4),c("G-Meta + Sw-Qf","G-Meta + Sw-Qf - uncorr"),cex=labcex,srt=45,adj=1,xpd=NA)

  tdm5b <- tdm5c <- tdm5[odm4]
  pd <- which(tdm5c>0)
  nd <- which(tdm5c<0)
  tdm5c[pd] <- "salmon1"
  tdm5b[pd] <- "chocolate3"
  tdm5c[nd] <- "skyblue1"
  tdm5b[nd] <- "steelblue3"

  n <- 4
  gmqfeft <- c()
  if (length(odm4) >0) {
    for(i in 1:length(odm4)){
      wid <- sum(!is.na(gmqfc[odm4[i],1:gmqfcnum]))
      temselid <- as.numeric(gmqfc[odm4[i],1:wid])
      temo <- order(pdm3[temselid],decreasing=T)
      temselid <- temselid[temo]
      temtdm <- as.numeric(tdm3[temselid]*(gmqfc[odm4[i],(1:wid)+gmqfcnum][temo]))
      temcol <- rep("black",wid)
      if(!is.na(opdm4[i])){
        rect(n,0,n+wid,opdm4[i],col="darkgoldenrod2",border="darkgoldenrod4",lwd=1.25)
        gmqfeft <- c(gmqfeft,temselid[(pdm3[temselid]>=1)&((temtdm*tdm5[odm4[i]])>0)])
        temcol[(pdm3[temselid]>=1)&((temtdm*tdm5[odm4[i]])>0)] <- "darkgoldenrod3"
      }
      rect(n,0,n+wid,pdm5[odm4[i]],col=tdm5c[i],border=tdm5b[i],lwd=1.25)
      if(wid!=1){
        pd <- which(temtdm>0)
        nd <- which(temtdm<0)
        temtdm[pd] <- "chocolate3"
        temtdm[nd] <- "steelblue3"
        rect(n:(n+wid-1),0,(n+1):(n+wid),pdm3[temselid],col=temtdm,border=T,density=25)
        text(n:(n+wid-1)+0.5,-spb,highcorinfo[temselid,3],col=temcol,cex=labcex,srt=45,adj=1,xpd=NA)
      } else{
        text(n+0.5,-spb,highcorinfo[gmqfc[odm4[i],1],3],col=temcol,cex=labcex,srt=45,adj=1,xpd=NA)
      }
      n <- n+wid
    }
  }

  if(length(setdiff(1:nrow(gmqfc),odm4))!=0){
    oid <- setdiff(1:nrow(gmqfc),odm4)
    ooid <- oid[order(pdm5[oid],decreasing=T)]

    n <- n+1
    for(i in 1:length(ooid)){
      wid <- sum(!is.na(gmqfc[ooid[i],1:gmqfcnum]))
      rect(n,0,n+wid,pdm5[ooid[i]],col="grey80",border="grey20",lwd=1.25)
      if(wid!=1){
        temselid <- as.numeric(gmqfc[ooid[i],1:wid])
        temo <- order(pdm3[temselid],decreasing=T)
        temselid <- temselid[temo]
        temtdm <- as.numeric(tdm3[temselid]*(gmqfc[ooid[i],(1:wid)+gmqfcnum][temo]))
        pd <- which(temtdm>0)
        zd <- which(temtdm==0)
        nd <- which(temtdm<0)
        temtdm[pd] <- "chocolate3"
        temtdm[zd] <- "grey20"
        temtdm[nd] <- "steelblue3"
        rect(n:(n+wid-1),0,(n+1):(n+wid),pdm3[temselid],col=temtdm,border=T,density=25)
        text(n:(n+wid-1)+0.5,-spb,highcorinfo[temselid,3],cex=labcex,srt=45,adj=1,xpd=NA)
        temymax <- max(c(pdm3[temselid],pdm5[ooid[i]]))
        lines(c(n,n),c(temymax+spb,temymax+1.5*spb),lwd=1.25,col="grey20")
        lines(c(n,n+wid),c(temymax+1.5*spb,temymax+1.5*spb),lwd=1.25,col="grey20")
        lines(c(n+wid,n+wid),c(temymax+spb,temymax+1.5*spb),lwd=1.25,col="grey20")
      } else{
        text(n+0.5,-spb,highcorinfo[gmqfc[ooid[i],1],3],cex=labcex,srt=45,adj=1,xpd=NA)
      }
      n <- n+wid
    }
  }


  pdm2 <- as.numeric(logtran(pdm[c(16:18,12)]))
  tdm4 <- as.numeric(pdm[c(18+1:cgwasenv$.TRAIT_NUM)])
  pdm4 <- as.numeric(logtran(pdm[c(18+cgwasenv$.TRAIT_NUM+1:cgwasenv$.TRAIT_NUM)]))

  opid <- order(pdm3,decreasing=T)

  ymax <- max(c(pdm2,pdm3,pdm4))
  spb <- ymax/50

  spa <- ymax/100
  spaa <- 0.05
  hpm <- (4000+50*cgwasenv$.TRAIT_NUM)*ymax/(cgwasenv$.TRAIT_NUM+7)/4500
  plot(0,col=F,xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(0,cgwasenv$.TRAIT_NUM+7),ylim=c(-spb-spa*3/2-4*hpm,ymax),xlab="",ylab=expression(-log[10](italic(p))),main=figname[ii],las=1,cex.axis=0.7,mgp=c(2.25,1,0))

  axis(2,las=1,at=yl[yl<=floor(ymax)],labels=as.character(yll[yl<=floor(ymax)]),lwd.ticks=1,cex.axis=0.7)
  rect(1:4,0,2:5,pdm2[1:4],col=c("brown2","darkgoldenrod2","chartreuse3","mediumpurple3"),border=c("brown4","darkgoldenrod4","green4","mediumpurple4"),lwd=1.25)
  text(1:4+0.5,rep(-spb,4),c("G-Meta","G-Meta + Sw-Qf","Sw-Qf","Min-GWAS"),cex=labcex,srt=45,adj=1,xpd=NA)

  tdm3b <- tdm3 <- tdm3[opid]
  pd <- which(tdm3>0)
  zd <- which(tdm3==0)
  nd <- which(tdm3<0)
  tdm3[pd] <- "salmon1"
  tdm3b[pd] <- "chocolate3"
  tdm3[zd] <- "grey20"
  tdm3b[zd] <- "grey20"
  tdm3[nd] <- "skyblue1"
  tdm3b[nd] <- "steelblue3"

  n <- 6
  for(i in 1:length(tdm3)){
    wid <- sum(!is.na(highcorc[opid[i],1:highcorcnum]))
    rect(n,0,n+wid,pdm3[opid[i]],col=tdm3[i],border=tdm3b[i],lwd=1.25)
    if(wid!=1){
      temselid <- as.numeric(highcorc[opid[i],1:wid])
      temselid <- temselid[order(pdm4[temselid],decreasing=T)]
      temtdm <- tdm4[temselid]
      pd <- which(temtdm>0)
      zd <- which(temtdm==0)
      nd <- which(temtdm<0)
      temtdm[pd] <- "chocolate3"
      temtdm[zd] <- "grey20"
      temtdm[nd] <- "steelblue3"
      rect(n:(n+wid-1),0,(n+1):(n+wid),pdm4[temselid],col=temtdm,border=T,density=25)
      temymax <- max(c(pdm4[temselid],pdm3[opid[i]]))
      lines(c(n,n),c(temymax+spb,temymax+1.5*spb),lwd=1.25,col="grey20")
      lines(c(n,n+wid),c(temymax+1.5*spb,temymax+1.5*spb),lwd=1.25,col="grey20")
      lines(c(n+wid,n+wid),c(temymax+spb,temymax+1.5*spb),lwd=1.25,col="grey20")
      text(n:(n+wid-1)+0.5,-2*spb-3*spa/2-4*hpm,cgwasenv$.TRAIT_NAME[temselid],cex=labcex,srt=45,adj=1,xpd=NA)

      temtdm <- pdm4[temselid]
      pd <- which(temtdm>=(-log10(0.05/effvar(orgcorm))))
      nd <- which(temtdm<(-log10(0.05/effvar(orgcorm))))
      temtdm[pd] <- "mediumpurple3"
      temtdm[nd] <- "grey80"
      rect(n:(n+wid-1)+spaa,-spb,(n+1):(n+wid)-spaa,-spb-hpm,col=temtdm,border=NA)

      if(length(which(gmeft==opid[i]))!=0){
        temtdm <- pdm4[temselid]
        pd <- which(pdm4[temselid]>=1)
        nd <- which(pdm4[temselid]<1)
        temtdm[pd] <- "brown2"
        temtdm[nd] <- "grey80"
        rect(n:(n+wid-1)+spaa,-spb-spa/2-hpm,(n+1):(n+wid)-spaa,-spb-spa/2-2*hpm,col=temtdm,border=NA)
      } else{
        rect(n:(n+wid-1)+spaa,-spb-spa/2-hpm,(n+1):(n+wid)-spaa,-spb-spa/2-2*hpm,col="grey80",border=NA)
      }
      if(length(which(gmqfeft==opid[i]))!=0){
        temtdm <- pdm4[temselid]
        pd <- which(pdm4[temselid]>=1)
        nd <- which(pdm4[temselid]<1)
        temtdm[pd] <- "darkgoldenrod2"
        temtdm[nd] <- "grey80"
        rect(n:(n+wid-1)+spaa,-spb-2*spa/2-2*hpm,(n+1):(n+wid)-spaa,-spb-2*spa/2-3*hpm,col=temtdm,border=NA)
      } else{
        rect(n:(n+wid-1)+spaa,-spb-2*spa/2-2*hpm,(n+1):(n+wid)-spaa,-spb-2*spa/2-3*hpm,col="grey80",border=NA)
      }
      if(length(which(qfeft==opid[i]))!=0){
        temtdm <- pdm4[temselid]
        pd <- which(pdm4[temselid]>=1)
        nd <- which(pdm4[temselid]<1)
        temtdm[pd] <- "chartreuse3"
        temtdm[nd] <- "grey80"
        rect(n:(n+wid-1)+spaa,-spb-3*spa/2-3*hpm,(n+1):(n+wid)-spaa,-spb-3*spa/2-4*hpm,col=temtdm,border=NA)
      } else{
        rect(n:(n+wid-1)+spaa,-spb-3*spa/2-3*hpm,(n+1):(n+wid)-spaa,-spb-3*spa/2-4*hpm,col="grey80",border=NA)
      }

    } else{
      text(n+0.5,-2*spb-3*spa/2-4*hpm,highcorinfo[opid[i],3],cex=labcex,srt=45,adj=1,xpd=NA)
      if(pdm3[opid[i]]>=(-log10(0.05/effvar(orgcorm)))){
        rect(n+spaa,-spb,n+1-spaa,-spb-hpm,col="mediumpurple3",border=NA)
      } else{
        rect(n+spaa,-spb,n+1-spaa,-spb-hpm,col="grey80",border=NA)
      }
      if(length(which(gmeft==opid[i]))!=0){
        rect(n+spaa,-spb-spa/2-hpm,n+1-spaa,-spb-spa/2-2*hpm,col="brown2",border=NA)
      } else{
        rect(n+spaa,-spb-spa/2-hpm,n+1-spaa,-spb-spa/2-2*hpm,col="grey80",border=NA)
      }
      if(length(which(gmqfeft==opid[i]))!=0){
        rect(n+spaa,-spb-2*spa/2-2*hpm,n+1-spaa,-spb-2*spa/2-3*hpm,col="darkgoldenrod2",border=NA)
      } else{
        rect(n+spaa,-spb-2*spa/2-2*hpm,n+1-spaa,-spb-2*spa/2-3*hpm,col="grey80",border=NA)
      }
      if(length(which(qfeft==opid[i]))!=0){
        rect(n+spaa,-spb-3*spa/2-3*hpm,n+1-spaa,-spb-3*spa/2-4*hpm,col="chartreuse3",border=NA)
      } else{
        rect(n+spaa,-spb-3*spa/2-3*hpm,n+1-spaa,-spb-3*spa/2-4*hpm,col="grey80",border=NA)
      }
    }
    n <- n+wid
  }

  abline(h=-log10(5e-8),col="grey20",lwd=0.8,lty=2)
  abline(h=-log10(0.05/effvar(orgcorm)),col="grey20",lwd=0.8,lty=2)

  dev.off()
}
