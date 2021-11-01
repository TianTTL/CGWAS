StatE1 <- function(traitid, cgwasenv) {
  df <- data.table::fread(cgwasenv$.GWAS_FILE_PATH[traitid],
                          header = T, stringsAsFactors = F, nThread = 1)
  df <- as.data.frame(df)
  bpm <- df[,cgwasenv$.ASSOC_COLUMN_INDEX[4:5]]
  data.table::fwrite(bpm,
                     file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid, ".bp")),
                     row.names = F, col.names = T, quote = F, nThread = 1)
  return(unique(which(is.na(bpm[,1])), which(is.na(bpm[,2]))))
}

StatE2 <- function(traitid, naid, cgwasenv) {
  bpm <- data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid, ".bp")),
                           header = T, stringsAsFactors = F, nThread = 1)
  bpm <- as.data.frame(bpm)
  if(length(naid) != 0) {
    bpm <- bpm[-naid,]
  }
  data.table::fwrite(as.data.frame(signif(bpm[,1], 7)),
         file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[traitid], ".beta")),
         row.names = F, col.names = F, quote = F, nThread = 1)
  bpm[,1] <- sign(bpm[,1])
  bpm[bpm[,2] < (1e-300),2] <- 1e-300
  data.table::fwrite(as.data.frame(signif(sqrt(qchisq(bpm[,2], 1, lower.tail=F))*bpm[,1], 7)),
         file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[traitid], ".stat")),
         row.names = F, col.names = F, quote = F, nThread = 1)
  file.remove(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(traitid, ".bp")))

  if (traitid == 1) {
    df <- data.table::fread(cgwasenv$.GWAS_FILE_PATH[1],
                            header = T, stringsAsFactors = F, nThread = 1)
    df <- as.data.frame(df)
    df.snp <- df[,cgwasenv$.ASSOC_COLUMN_INDEX[1:3]]
    if (length(naid) != 0) {
      df.snp <- df.snp[-naid,]
    }
    data.table::fwrite(df.snp,
                       file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SnpIndex"),
                       sep = " ", na = "NA", row.names = F, quote = F, nThread = 1)
  }
  return(nrow(bpm))
}

qcf <- function(lam, ppv, tq, meaq, inseq) {
  qcr <- rep(NA, length(lam))
  for(i in 1:length(qcr)) {
    if(ppv[i]==0) {
      qcr[i] <- mean((pchisq(tq/lam[i], 1)*(1-ppv[i])-inseq)^2)
    } else{
      qcr[i] <- mean((pchisq(tq/lam[i], 1)*(1-ppv[i])+pchisq(tq/(lam[i]+(meaq-lam[i])/ppv[i]), 1)*ppv[i]-inseq)^2)
    }
  }
  return(qcr)
}

odgridse <- function(rb, intv, lamsp, tq, meaq, N, inseq) {
  while(intv>1e-4) {
    cid <- order(rb)[1]
    if(cid==1) {
      lamsp <- seq(lamsp[cid], lamsp[cid+1], length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1, 5)] <- rb[c(cid, cid+1)]
      rb[3] <- qcf(lamsp[3], propv[3], tq, meaq, inseq)
      intv <- intv/4
    } else if(cid==5) {
      lamsp <- seq(lamsp[cid-1], lamsp[cid], length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1, 5)] <- rb[c(cid-1, cid)]
      rb[3] <- qcf(lamsp[3], propv[3], tq, meaq, inseq)
      intv <- intv/4
    } else{
      lamsp <- seq(lamsp[cid-1], lamsp[cid+1], length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1, 3, 5)] <- rb[c(cid-1, cid, cid+1)]
      intv <- intv/2
    }
    rb[c(2, 4)] <- qcf(lamsp[c(2, 4)], propv[c(2, 4)], tq, meaq, inseq)
  }
  return(lamsp[order(rb)[1]])
}

calinf <- function(tm2, inseq, mp) {
  meaq <- mean(tm2)
  vq <- mean(tm2^2)-meaq^2
  N <- vq-2*meaq^2
  tlam <- median(tm2)/qchisq(0.5, 1)
  if((N<=0)|((meaq-1)<=0)) {
    lamv <- max(1, min(meaq, tlam))
    propv <- 0
  } else{
    tq <- quantile(tm2, inseq)
    lamsp <- seq(max(meaq-sqrt(N/3/(1/mp-1)), 0.99), meaq+0.01, length.out=5)
    propv <- 1/(N/3/(meaq-lamsp)^2+1)
    rb <- qcf(lamsp, propv, tq, meaq, inseq)
    intv <- (lamsp[5]-lamsp[1])/2
    lamv <- odgridse(rb, intv, lamsp, tq, meaq, N, inseq)
    propv <- 1/(N/3/(meaq-lamv)^2+1)
    if(lamv>meaq) {
      lamsp <- seq(max(meaq-sqrt(N/3/(1/mp-1)), 0.99), meaq, length.out=5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb <- qcf(lamsp, propv, tq, meaq, inseq)
      intv <- (lamsp[5]-lamsp[1])/2
      lamv <- odgridse(rb, intv, lamsp, tq, meaq, N, inseq)
      propv <- 1/(N/3/(meaq-lamv)^2+1)
    }
  }
  return(c(lamv, propv))
}

Coresti.ICE <- function(t1, t2, id1, id2, resinfm, ranidlist, cgwasenv) {
  mt2 <- resinfm[c(id1, id2), 4]
  ax1 <- t1^2
  ax2 <- t2^2
  atmc1 <- ax1+ax2
  atmc2 <- t1*t2
  newpv <- c()
  chisq.cv <- qchisq(0.5, 2)
  for(i in 1:nrow(ranidlist)) {
    x1 <- ax1[ranidlist[i,]]
    x2 <- ax2[ranidlist[i,]]
    tmc1 <- atmc1[ranidlist[i,]]
    tmc2 <- atmc2[ranidlist[i,]]
    testid <- testcor <- 0
    selid <- tmc1 < chisq.cv
    cv <- newp <- mean(tmc2[selid])/sqrt(mean(x1[selid])*mean(x2[selid]))
    if(abs(newp-testcor)>1e-4) {
      testcor <- newp
      testid <- c(testid, testcor)
      tv <- (tmc1-2*testcor*tmc2)/(1-testcor^2)
      selid <- tv < chisq.cv
      newp <- mean(tmc2[selid])/sqrt(mean(x1[selid])*mean(x2[selid]))
      while(abs(newp-testcor)>1e-4) {
        cv <- c(cv, newp)
        rcc <- summary(lm(cv~testid))$coefficients[,1]
        testcor <- rcc[1]/(1-rcc[2])
        testid <- c(testid, testcor)
        tv <- (tmc1-2*testcor*tmc2)/(1-testcor^2)
        selid <- tv < chisq.cv
        newp <- mean(tmc2[selid])/sqrt(mean(x1[selid])*mean(x2[selid]))
        if(length(testid)>5) {
          testid <- testid[-1]
          cv <- cv[-1]
        }
      }
    }
    newpv <- c(newpv, newp)
  }
  newp <- mean(newpv)

  mxx <- mean(atmc2)
  efc <- mxx-newp
  if(efc==0) {
    efcr <- 0
  } else if(prod(mt2-1)!=0) {
    efcr <- efc/sqrt(prod(mt2-1))
    if(efcr>1) {
      efcr <- 1
    }
    if(efcr<(-1)) {
      efcr <- -1
    }
  } else{
    efcr <- efc/abs(efc)
  }

  tv <- (atmc1-2*newp*atmc2)/(1-newp^2)
  selid <- tv < qchisq(cgwasenv$.P_THRD, 2, lower.tail=F)
  meff <- mt2-c(mean(ax1[selid]), mean(ax2[selid]))
  meff[meff<0] <- 0
  mcef <- mxx-mean(atmc2[selid])
  if(mcef==0){
    mcefr <- 0
    if((meff[1]==0)&(meff[2]==0)){
      meff[1:2] <- 1e-6
    }
  } else if(prod(meff)!=0){
    mcefr <- mcef/sqrt(prod(meff))
    if(mcefr>1){
      mcefr <- 1
    }
    if(mcefr<(-1)){
      mcefr <- -1
    }
  } else{
    mcefr <- mcef/abs(mcef)
  }

  return(c(mxx/sqrt(prod(mt2)), newp, efc, efcr, meff, mcef, mcefr))
}

CorE.ICE <- function(pgpid, upnum, n, pairma, resinfm, ranidlist, cgwasenv) {
  corm <- c()
  tempsv <- c(0, 0)
  for(i in (upnum*(pgpid-1)+1):min(upnum*pgpid, n)) {
    if(pairma[i, 1]!=tempsv[1]) {
      t1 <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[pairma[i, 1]], ".efstat")),
                                        header=F, nThread = 1))[,1]
      tempsv[1] <- pairma[i,1]
    }
    if(pairma[i, 2]!=tempsv[2]) {
      t2 <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[pairma[i, 2]], ".efstat")),
                                        header=F, nThread = 1))[,1]
      tempsv[2] <- pairma[i,2]
    }
    corm <- rbind(corm, Coresti.ICE(t1, t2, pairma[i, 1], pairma[i, 2], resinfm, ranidlist, cgwasenv))
  }
  return(corm)
}

ridl.ICE <- function(repn, minsnpn, cgwasenv) {
  tm <- (as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[1], ".stat")),
                                     header=F, nThread = 1))[,1])^2
  bn <- 1
  while((length(tm)/bn)>minsnpn) {
    bn <- bn+1
  }
  sl <- floor(length(tm)/bn)
  idm <- matrix(NA, repn*bn, sl)
  for(i in 1:repn) {
    idv <- sample(1:length(tm), length(tm))
    for(ii in 1:bn) {
      idm[(i-1)*bn+ii,] <- seq(ii, length(tm), bn)[1:sl]
    }
  }
  return(idm)
}

Coresti.ebico <- function(t1, t2, id1, id2, ranidlist, cgwasenv) {
  mt2 <- as.numeric(c(id1, id2)) + 1
  ax1 <- t1^2
  ax2 <- t2^2
  atmc1 <- ax1+ax2
  atmc2 <- t1*t2
  newpv <- c()
  chisq.cv <- qchisq(0.5, 2)
  for(i in 1:nrow(ranidlist)) {
    x1 <- ax1[ranidlist[i,]]
    x2 <- ax2[ranidlist[i,]]
    tmc1 <- atmc1[ranidlist[i,]]
    tmc2 <- atmc2[ranidlist[i,]]
    testid <- testcor <- 0
    selid <- tmc1 < chisq.cv
    cv <- newp <- mean(tmc2[selid])/sqrt(mean(x1[selid])*mean(x2[selid]))
    if(abs(newp-testcor)>1e-4) {
      testcor <- newp
      testid <- c(testid, testcor)
      tv <- (tmc1-2*testcor*tmc2)/(1-testcor^2)
      selid <- tv < chisq.cv
      newp <- mean(tmc2[selid])/sqrt(mean(x1[selid])*mean(x2[selid]))
      while(abs(newp-testcor)>1e-4) {
        cv <- c(cv, newp)
        rcc <- summary(lm(cv~testid))$coefficients[,1]
        testcor <- rcc[1]/(1-rcc[2])
        testid <- c(testid, testcor)
        tv <- (tmc1-2*testcor*tmc2)/(1-testcor^2)
        selid <- tv < chisq.cv
        newp <- mean(tmc2[selid])/sqrt(mean(x1[selid])*mean(x2[selid]))
        if(length(testid)>5) {
          testid <- testid[-1]
          cv <- cv[-1]
        }
      }
    }
    newpv <- c(newpv, newp)
  }
  newp <- mean(newpv)

  mxx <- mean(atmc2)
  efc <- mxx-newp
  if(efc==0) {
    efcr <- 0
  } else if(prod(mt2-1)!=0) {
    efcr <- efc/sqrt(prod(mt2-1))
    if(efcr>1) {
      efcr <- 1
    }
    if(efcr<(-1)) {
      efcr <- -1
    }
  } else{
    efcr <- efc/abs(efc)
  }

  tv <- (atmc1-2*newp*atmc2)/(1-newp^2)
  selid <- tv<qchisq(cgwasenv$.P_THRD, 2, lower.tail=F)
  meff <- mt2-c(mean(ax1[selid]), mean(ax2[selid]))
  meff[meff<0] <- 0
  mcef <- mxx-mean(atmc2[selid])
  if(mcef==0){
    mcefr <- 0
    if((meff[1]==0)&(meff[2]==0)){
      meff[1:2] <- 1e-6
    }
  } else if(prod(meff)!=0){
    mcefr <- mcef/sqrt(prod(meff))
    if(mcefr>1){
      mcefr <- 1
    }
    if(mcefr<(-1)){
      mcefr <- -1
    }
  } else{
    mcefr <- mcef/abs(mcef)
  }

  return(c(mt2-1, mxx/sqrt(prod(mt2)), newp, efc, efcr, meff, mcef, mcefr))
}

CorE.ebico <- function(traitid, newt, TNm, x2m, did, ranidlist, cgwasenv) {
  t1 <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH,
                                              paste0(TNm[traitid, 1], ".efstat")),
                                    header=F, nThread = 1))[,1]
  corm <- Coresti.ebico(t1, newt, TNm[traitid, 2], x2m[did], ranidlist, cgwasenv)
  return(corm)
}

ebicocof <- function(cv, ssw, tm, thresc) {
  solbm <- solve(matrix(c(1, cv[4], cv[4], 1), 2))
  em <- matrix(c(cv[1], cv[5], cv[5], cv[2]), 2)
  fem <- matrix(c(cv[1], sqrt(cv[1]*cv[2]), sqrt(cv[1]*cv[2]), cv[2]), 2)
  esm <- matrix(c(cv[7], cv[9], cv[9], cv[8]), 2)
  bsm <- sqrt(c(ssw[1], ssw[2]))
  d011 <- t(c(1, 1)*bsm)
  cbef11 <- t(d011%*%solbm)/as.numeric(sqrt(d011%*%solbm%*%t(d011)))
  d021 <- t(c(1, -1)*bsm)
  cbef21 <- t(d021%*%solbm)/as.numeric(sqrt(d021%*%solbm%*%t(d021)))
  if(length(which((tm%*%cbef11)^2>thresc))>=length(which((tm%*%cbef21)^2>thresc))) {
    d01 <- d011/abs(d011)
    cbef1 <- cbef11
    bbef1 <- t(d011%*%solbm)
  } else{
    d01 <- d021/abs(d021)
    cbef1 <- cbef21
    bbef1 <- t(d021%*%solbm)
  }
  d011 <- t(c(1, 1)/bsm)%*%fem
  cbef11 <- t(d011%*%solbm)/as.numeric(sqrt(d011%*%solbm%*%t(d011)))
  d021 <- t(c(1, -1)/bsm)%*%fem
  cbef21 <- t(d021%*%solbm)/as.numeric(sqrt(d021%*%solbm%*%t(d021)))
  if(length(which((tm%*%cbef11)^2>thresc))>=length(which((tm%*%cbef21)^2>thresc))) {
    d02 <- d011/abs(d011)
    cbef2 <- cbef11
    bbef2 <- t(d011%*%solbm)
  } else{
    d02 <- d021/abs(d021)
    cbef2 <- cbef21
    bbef2 <- t(d021%*%solbm)
  }
  sid <- order(c(sum((tm%*%cbef2)^2>thresc), sum((tm%*%cbef1)^2>thresc)))[2]
  d0 <- cbind(t(d02), t(d01))[,sid]
  cbef <- cbind(cbef2, cbef1)[,sid]
  bbef <- cbind(bbef2, bbef1)[,sid]

  d11 <- t(c(1, 1)/bsm)%*%em
  coef1 <- t(d11%*%solbm)/as.numeric(sqrt(d11%*%solbm%*%t(d11)))
  d12 <- t(c(1, -1)/bsm)%*%em
  coef2 <- t(d12%*%solbm)/as.numeric(sqrt(d12%*%solbm%*%t(d12)))
  if(length(which((tm%*%coef1)^2>thresc))>=length(which((tm%*%coef2)^2>thresc))){
    boef01 <- t(d11%*%solbm)
    coef01 <- coef1
    d101 <- d11/abs(d11)
  } else{
    boef01 <- t(d12%*%solbm)
    coef01 <- coef2
    d101 <- d12/abs(d12)
  }
  d1 <- t(d101)
  coef <- coef01
  boef <- boef01

  d21 <- t(c(1, 1)/bsm)%*%esm
  coefs1 <- t(d21%*%solbm)/as.numeric(sqrt(d21%*%solbm%*%t(d21)))
  d22 <- t(c(1, -1)/bsm)%*%esm
  coefs2 <- t(d22%*%solbm)/as.numeric(sqrt(d22%*%solbm%*%t(d22)))
  if(length(which((tm%*%coefs1)^2>thresc))>=length(which((tm%*%coefs2)^2>thresc))){
    boefs01 <- t(d21%*%solbm)
    coefs01 <- coefs1
    d201 <- d21/abs(d21)
  } else{
    boefs01 <- t(d22%*%solbm)
    coefs01 <- coefs2
    d201 <- d22/abs(d22)
  }
  d2 <- t(d201)
  coefs <- coefs01
  boefs <- boefs01

  return(cbind(coef, coefs, cbef, boef, boefs, bbef, d1, d2, d0))
}

Essfun <- function(tid, mafv, cgwasenv) {
  t2m <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(tid, ".efstat")),
                                         header=F, stringsAsFactors=F, nThread = 1))[,1]
  b2m <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(tid, ".beta")),
                                         header=F, stringsAsFactors=F, nThread = 1))[,1]
  s2m <- t2m/b2m
  if(cgwasenv$.MAF_FILE_EXIST) {
    mse <- median(s2m^2/2/mafv/(1-mafv), na.rm=T)
  } else{
    mse <- median(s2m^2, na.rm=T)
  }
  return(mse)
}

ridl.ebico <- function(repn, minsnpn, cgwasenv) {
  tm <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[1], ".efstat")),
                                         header=F, stringsAsFactors=F, nThread = 1))[,1]
  bn <- 1
  while((length(tm)/bn)>minsnpn) {
    bn <- bn+1
  }
  sl <- floor(length(tm)/bn)
  idm <- matrix(NA, repn*bn, sl)
  for(i in 1:repn) {
    idv <- sample(1:length(tm), length(tm))
    for(ii in 1:bn) {
      idm[(i-1)*bn+ii,] <- seq(ii, length(tm), bn)[1:sl]
    }
  }
  return(idm)
}

gcrom <- function(rv, n) {
  corm <- diag(n)
  n <- 0
  for(i in 2:nrow(corm)) {
    corm[i:nrow(corm), i-1] <- corm[i-1, i:nrow(corm)] <- rv[(n+1):(n+nrow(corm)-i+1)]
    n <- n+nrow(corm)-i+1
  }
  return(corm)
}

ptc <- function(bgm, stm) {
  nn <- 0
  cm <- 0.01*(stm-bgm)
  while(!all(eigen(bgm)$values>0)) {
    nn <- nn+1
    bgm <- bgm+cm
  }
  while((!all(abs(solve(bgm))<10))&(nn<100)) {
    nn <- nn+1
    bgm <- bgm+cm
  }
  return(list(bgm, nn, max(abs(solve(bgm)))))
}

mkdf <- function(n, cm) {
  df <- MASS::mvrnorm(n, mu = rep(0, ncol(cm)), Sigma = cm)
  return(df)
}

trtvf <- function(tv, cm) {
  if (length(tv) == 1) {
    tmp <- as.numeric(1 / cm * tv * tv)
  } else {
    tmp <- as.numeric(t(tv) %*% chol2inv(chol(cm)) %*% tv)
  }
  return(pchisq(tmp,
                length(tv),
                lower.tail = F))
}

swtrtsimu <- function(snpn, cmg, cm, ACstatm, maxcn, cutoff.thv, cgwasenv) {
  tm <- mkdf(snpn, cmg)
  ebcm <- matrix(NA, nrow(tm), nrow(ACstatm))
  for(i in 1:nrow(ACstatm)){
    compn <- match(as.character(ACstatm[i,
                                        (4:(3+maxcn))
                                          [!is.na(ACstatm[i, 4:(3+maxcn)])]]),
                   cgwasenv$.TRAIT_NAME)
    if(length(compn)==1){
      ebcm[,i] <- tm[, compn]
    } else{
      ebcm[,i] <- tm[, compn] %*% as.numeric(ACstatm[i,
                                                     ((4+2*maxcn):(3+3*maxcn))
                                                       [!is.na(ACstatm[i, (4+2*maxcn):(3+3*maxcn)])]])
    }
  }
  return(swtrt(ebcm, cm, cutoff.thv))
}

swtrt <- function(tm, cm, cutoff.thv) {
  resv <- matrix(NA, nrow(tm), length(cutoff.thv)+1)
  nresv <- length(cutoff.thv)
  fnresv <- length(cutoff.thv) + 1
  for(i in 1:nrow(tm)) {
    tv <- tm[i,]
    tv2 <- tv^2
    cm.tmp <- cm
    rid <- tv2 > cutoff.thv[1]
    ridN <- sum(rid)
    if(ridN == 0) {
      resv[i, 1:nresv] <- NA
      resv[i, fnresv] <- pchisq(max(tv2), 1, lower.tail=F)
      next
    } else if(ridN==length(tv)) {
      resv[i, 1] <- trtvf(tv, cm.tmp)
    } else {
      tv <- tv[rid]
      tv2 <- tv2[rid]
      cm.tmp <- cm.tmp[rid, rid]
      resv[i, 1] <- trtvf(tv, cm.tmp)
    }
    for(j in 2:length(cutoff.thv)) {
      rid <- tv2 > cutoff.thv[j]
      ridN <- sum(rid)
      if(ridN == 0) {
        resv[i, j:nresv] <- NA
        break
      } else if(ridN==length(tv)) {
        resv[i, j] <- resv[i, j-1]
      } else {
        tv <- tv[rid]
        tv2 <- tv2[rid]
        cm.tmp <- cm.tmp[rid, rid]
        resv[i, j] <- trtvf(tv, cm.tmp)
      }
    }
    resv[i, fnresv] <- pchisq(max(tv2), 1, lower.tail=F)
  }
  return(resv)
}

appmin <- function(pm) {
  return(do.call(pmin, c(as.data.frame(pm), na.rm = TRUE)))
}

locisearch <- function(sp, spos, fp, cgwasenv) {
  n <- 0
  lclist <- list()
  if(length(sp)!=0) {
    keyid <- order(sp)[1]
    while(sp[keyid]<=fp) {
      idirc <- ddirc <- spos[keyid,2]
      oldddirc <- oldidirc <- -1
      while((oldddirc!=ddirc)|(oldidirc!=idirc)) {
        oldddirc <- ddirc
        oldidirc <- idirc
        tdl <- which((spos[,1]==spos[keyid, 1])&(spos[,2]>(ddirc-cgwasenv$.LOCI_INTER))&(spos[,2]<(idirc+cgwasenv$.LOCI_INTER)))
        ddirc <- spos[min(tdl),2]
        idirc <- spos[max(tdl),2]
      }
      n <- n+1
      lclist[[n]] <- tdl
      spos[tdl,] <- NA
      sp[tdl] <- NA
      keyid <- order(sp)[1]
      if(is.na(sp[keyid])) {
        break
      }
    }
  }
  return(lclist)
}

fdrf <- function(opv, Sind, cgwasenv) {
  os <- which(opv<1e-4)
  stos <- opv[os][order(opv[os])]
  fdrvo <- c()
  for(i in 1:length(stos)) {
    fdrvo <- c(fdrvo, stos[i]*nrow(Sind)/i)
  }
  mfdrvo <- min(fdrvo)

  efplcno <- lcno <- fdrsno <- fdrpo <- c()
  for(i in 1:length(cgwasenv$.FDR_SET)) {
    if(mfdrvo<cgwasenv$.FDR_SET[i]) {
      id <- max(which(fdrvo<cgwasenv$.FDR_SET[i]))
      if(id!=length(fdrvo)) {
        fdrpo <- c(fdrpo, mean(stos[c(id, id+1)]))
        fdrsno <- c(fdrsno, id)
        oso <- os[which(opv[os]<=fdrpo[i])]
        lco <- locisearch(opv[oso], Sind[oso, 1:2], fdrpo[i], cgwasenv)
        lcno <- c(lcno, length(lco))
        efplcno <- c(efplcno, lcno[i]*cgwasenv$.FDR_SET[i])
      } else{
        fdrpo <- c(fdrpo, 1e-4)
        fdrsno <- c(fdrsno, NA)
        lcno <- c(lcno, NA)
        efplcno <- c(efplcno, NA)
      }
    } else{
      fdrpo <- c(fdrpo, NA)
      fdrsno <- c(fdrsno, 0)
      lcno <- c(lcno, 0)
      efplcno <- c(efplcno, NA)
    }
  }
  return(cbind(fdrpo, fdrsno, lcno, efplcno))
}

efvn <- function(cm, snpn) {
  tm <- scale(mkdf(snpn, cm))
  tm2 <- tm^2
  tm2.max <- do.call(pmax, as.data.frame(tm2))
  ttop <- pchisq(tm2.max, 1, lower.tail = F)
  return(ttop)
}

minv <- function(tm) {
  tm2 <- tm^2
  tm2.max <- do.call(pmax, as.data.frame(tm2))
  ttop <- pchisq(tm2.max, 1, lower.tail = F)
  return(ttop)
}

calnna <- function(a, ppn) {
  a <- a[order(a)]
  tq <- sum(!is.na(a))
  tq005 <- tq * 0.05
  return(log(1 - mean(ppn[c(floor(tq005),
                            ceiling(tq005))]),
             1 - mean(a[c(floor(tq005),
                          ceiling(tq005))])))
}

tpcor.m <- function(m, efn) {
  for (i in 1: length(efn)) {
    v <- m[, i]
    m[which(v > 1e-12), i] <- 1 - (1 - v[which(v > 1e-12)]) ^ efn[i]
    m[which(v <= 1e-12), i] <- v[which(v <= 1e-12)] * efn[i]
  }
  return(m)
}

tpcor.v <- function(v, efn) {
  vt <- v
  v[which(vt > 1e-12)] <- 1 - (1 - vt[which(vt > 1e-12)]) ^ efn
  v[which(vt <= 1e-12)] <- vt[which(vt <= 1e-12)] * efn
  return(v)
}

tpcor2 <- function(v, md){
  lgv <- -log10(v)
  intp <- md[[length(md)]]
  v[lgv<=intp[1]] <- predict(md[[1]], intp[1])
  for(i in 2:length(md)){
    v[(lgv>intp[i-1])&(lgv<=intp[i])] <- predict(md[[i-1]], lgv[(lgv>intp[i-1])&(lgv<=intp[i])])
  }
  v[lgv>intp[length(intp)]] <- predict(md[[length(md)-1]], intp[length(intp)])
  if(length(md)>2){
    for(i in 2:(length(md)-1)){
      v[(lgv>(intp[i]-0.5))&(lgv<=(intp[i]+0.5))] <-
        ((intp[i] - lgv[(lgv>(intp[i]-0.5)) & (lgv<=(intp[i]+0.5))] + 0.5) *
          predict(md[[i-1]], lgv[(lgv > (intp[i]-0.5)) & (lgv <= (intp[i]+0.5))]) +
          (lgv[(lgv>(intp[i]-0.5)) & (lgv<=(intp[i]+0.5))] - intp[i] + 0.5) *
          predict(md[[i]], lgv[(lgv > (intp[i]-0.5)) & (lgv <= (intp[i]+0.5))]))
    }
  }
  return(v)
}

nullcorrection <- function(isimup, nam, cgwasenv) {
  options(warn=-1) # turn off warning temporarily
  sn <- cgwasenv$.SIMUL_N * cgwasenv$.SIMUL_SNP_N
  rsn <- sn/cgwasenv$.IND_SNP_N
  swwv <- isimup[order(isimup)]

  sid <- rsn*c(1:(cgwasenv$.IND_SNP_N-1))
  y <- swwv[sid]
  x <- sid/sn
  yy <- x/y

  rth <- seq(yy[cgwasenv$.IND_SNP_N-1], yy[cgwasenv$.IND_SNP_N/10], length.out=(1/0.025+1))
  rth <- rth[-length(rth)]
  ii <- 1
  sid <- c()
  for(i in length(yy):1){
    if(yy[i]>=rth[ii]){
      ii <- ii+1
      sid <- c(sid, i)
    }
    if(ii>length(rth)){
      break
    }
  }

  sid <- sort(unique(c(sid, round(10^seq(0, -log10(cgwasenv$.IND_SNP_N), -0.025)*(cgwasenv$.IND_SNP_N-1)))), decreasing=T)
  y <- y[sid]
  x <- x[sid]
  write.table(cbind(x, y),
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, paste0("Null", nam, "correction")),
              row.names=F, col.names=c("ExpectedQuantile&P", "ObservedP"), quote=F)

  yy <- x/y
  xx <- -log10(x)

  smn <- sum(sid<10)
  weight <- seq(0, 1, length.out=(smn+2))[2:(smn+1)]

  tv <- c()
  for(i in smn:1){
    tv <- c(tv, sum(yy[(length(yy)-smn+1):(length(yy)-i+1)]*weight[1:(smn-i+1)])/sum(weight[1:(smn-i+1)]))
  }

  for(i in cgwasenv$.LOESS_SPAN_V){
    drawnc(xx, yy, tv, i, nam, cgwasenv)
  }

  intp <- c(min(xx), -log10(cgwasenv$.LOESS_INTER_P), max(xx))

  jpeg(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, paste0("Null", nam, "correction.jpg")),
       width=3000, height=3000, res=600)
  par(mar=c(3, 3, 1, 1))
  plot(xx[(length(yy)-length(tv)+1):length(yy)], yy[(length(yy)-length(tv)+1):length(yy)], pch=20, col="grey75", xlim=c(0, max(xx)), ylim=c(0, max(c(yy, tv))))
  yy[(length(yy)-length(tv)+1):length(yy)] <- tv
  md <- list()
  for(i in 1:length(cgwasenv$.LOESS_SPAN_V)){
    md[[i]] <- loess(yy~xx, span=cgwasenv$.LOESS_SPAN_V[i])
  }
  md[[length(cgwasenv$.LOESS_SPAN_V)+1]] <- intp
  points(xx, yy, pch=20)
  lines(c(seq(min(xx), max(xx), 0.001), max(xx)), tpcor2(10^(-c(seq(min(xx), max(xx), 0.001), max(xx))), md), lwd=1.5, col="red")
  dev.off()

  ssid <- c(1:1e4, seq(1e4, 1e5, 10), seq(1e5, cgwasenv$.IND_SNP_N, 100))
  swwm <- matrix(isimup[sample(1:sn, sn)], ncol=rsn)
  for(i in 1:rsn){
    swwm[,i] <- swwm[order(swwm[,i]),i]
  }
  infm <- t(apply(swwm[ssid,], 1, quantile, probs=c(0.05, 0.5, 0.95)))

  tpq <- seq(1/cgwasenv$.IND_SNP_N, 1-1/cgwasenv$.IND_SNP_N, length.out=cgwasenv$.IND_SNP_N)
  swwm <- matrix(isimup[sample(1:sn, sn)], ncol=rsn)
  for(i in 1:rsn){
    swwm[,i] <- swwm[order(swwm[,i]),i]
  }
  ct <- tpcor2(tpq, md)
  swwm <- swwm*ct
  swwm[swwm>1] <- 1
  nulm <- t(apply(swwm[ssid,], 1, quantile, probs=c(0.05, 0.5, 0.95)))

  nulpm <- cbind(qbeta(0.05, ssid, cgwasenv$.IND_SNP_N+1-ssid), qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid), qbeta(0.95, ssid, cgwasenv$.IND_SNP_N+1-ssid))

  jpeg(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, paste0("Null", nam, "distribution.jpg")),
       width=3000, height=3000, res=600)
  par(mar=c(3.5, 3.5, 1, 1))

  plot(NA, xlim=c(0, max(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)))), ylim=c(0, max(c(-log10(infm[1, 1]), -log10(nulpm[1, 1]), -log10(nulm[1, 1])))), xlab=expression(paste("Expected   ", -Log[10](italic(p)))), ylab=expression(paste("Observed   ", -Log[10](italic(p)))), mgp=c(2, 0.7, 0), las=1, cex.axis=0.85, tck=-0.015)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(infm[,1]), lwd=1.8, col="grey75", lty=3)
  points(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(infm[,2]), pch=20, col="grey75", cex=0.6)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(infm[,3]), lwd=1.8, col="grey75", lty=3)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulpm[,1]), lwd=1.8, col="black", lty=3)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulpm[,3]), lwd=1.8, col="black", lty=3)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulm[,1]), lwd=1.8, col="#FD9001", lty=3)
  points(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulm[,2]), pch=20, col="#FD9001", cex=0.6)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulm[,3]), lwd=1.8, col="#FD9001", lty=3)
  lines(c(-log10(qbeta(0.5, cgwasenv$.IND_SNP_N+1, 1)), -log10(qbeta(0.5, 1, cgwasenv$.IND_SNP_N+1))), c(-log10(qbeta(0.5, cgwasenv$.IND_SNP_N+1, 1)), -log10(qbeta(0.5, 1, cgwasenv$.IND_SNP_N+1))), col="black")
  dev.off()

  options(warn=0) # turn on warning
  return(md)
}

drawnc <- function(xx, yy, tv, i, nam, cgwasenv){
  jpeg(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, paste0("Null", nam, "correctionSpan", i, ".jpg")),
       width=3000, height=3000, res=600)
  par(mar=c(3, 3, 1, 1))
  plot(xx[(length(yy)-length(tv)+1):length(yy)], yy[(length(yy)-length(tv)+1):length(yy)], pch=20, col="grey75", xlim=c(0, max(xx)), ylim=c(0, max(c(yy, tv))))
  yy[(length(yy)-length(tv)+1):length(yy)] <- tv
  md <- loess(yy~xx, span=i)
  points(xx, yy, pch=20)
  lines(c(seq(min(xx), max(xx), 0.001), max(xx)), predict(md, c(seq(min(xx), max(xx), 0.001), max(xx))), lwd=1.2, col="red")
  dev.off()
}

ffdrv <- function(pv, fv, Sind) {
  os <- which(pv<1e-4)
  stos <- pv[os][order(pv[os])]
  fdrvo <- c()
  for(i in 1:length(stos)) {
    fdrvo <- c(fdrvo, stos[i]*nrow(Sind)/i)
  }
  if(min(fdrvo)<fv) {
    id <- max(which(fdrvo<=fv))
    if(mean(stos[c(id, id+1)])>5e-8) {
      return(mean(stos[c(id, id+1)]))
    } else{
      return(5e-8)
    }
  } else{
    return(5e-8)
  }
}

ptc.summay <- function(bgm, stm) {
  nn <- 0
  cm <- 0.01*(stm-bgm)
  while(!all(eigen(bgm)$values>0)) {
    nn <- nn+1
    bgm <- bgm+cm
  }
  return(list(bgm, nn))
}

ttp <- function(tm) {
  return(pchisq(tm^2, 1, lower.tail=F))
}

tpcor.summary <- function(v, efn) {
  vt <- v
  vt[which(v > 1e-12)] <- 1 - (1 - v[which(v > 1e-12)]) ^ efn
  vt[which(v <= 1e-12)] <- v[which(v <= 1e-12)] * efn
  return(vt)
}

manhattan <- function(op,np,osid,nsid,lco,lcn,fpo,fpn,
                      Sind,newtick,
                      pcol=c("#009DE8", "#233FAA"), pcex=0.35, ph=20) {
  d <- cbind(Sind[,c(1, 4)], -log10(op), -log10(np))
  colnames(d) <- c("CHR", "pos", "OP", "NP")

  ytop <- ceiling(max(d[,3:4]))
  if(ytop>30) {
    d$OP[d$OP>30]=d$OP[d$OP>30]/2+15
    d$NP[d$NP>30]=d$NP[d$NP>30]/2+15
    if(ytop>70) {
      d$OP[d$OP>50]=d$OP[d$OP>50]/1.5+50/3
      d$NP[d$NP>50]=d$NP[d$NP>50]/1.5+50/3
      if(ytop>100) {
        d$OP[d$OP>60]=d$OP[d$OP>60]*0.3+42
        d$NP[d$NP>60]=d$NP[d$NP>60]*0.3+42
        if(ytop>300) {
          d$OP[d$OP>80]=80
          d$NP[d$NP>80]=80
        }
      }
    }
  }
  ymax <- ceiling(max(d[,3:4]))

  xmax <- max(d$pos)*1.03
  xmin <- max(d$pos)*-0.03

  l <- c(0, 5, 10, 15, 20, 25, seq(30, 80, 10))
  ll <- c(0, 5, 10, 15, 20, 25, 30, 50, 70, 100, 200, 300)
  tid <- which(l<ymax)
  ll <- ll[c(tid, length(tid)+1)]
  l <- l[c(tid, length(tid)+1)]

  ymax <- max(l)
    plot(0, col=F, xaxt='n',yaxt='n',bty='n',xaxs='i',xlim=c(xmin, xmax), ylim=c(-ymax, 1.2*ymax), xlab="", ylab=expression(-log[10](italic(p))), mgp=c(2, 0.7, 0), cex.lab=1.1)

  axis(2, las=1, at=l+0.2*ymax, labels=as.character(ll), tck=-0.02, mgp=c(2, 0.7, 0), lwd.ticks=1, cex.axis=0.85)
  axis(2, las=1, at=-l, labels=as.character(ll), tck=-0.02, mgp=c(2, 0.7, 0), lwd.ticks=1, cex.axis=0.85)

  pcol <- rep(pcol, max(d$CHR))[1:max(d$CHR)]

  d1 <- which((d$NP>0)&(d$NP<=1))
  d2 <- which((d$NP>1)&(d$NP<=2))
  d3 <- which((d$NP>2)&(d$NP<=3))
  d4 <- which(d$NP>3)
  tsid <- c(sample(d1, 1e5), sample(d2, min(length(d2), 1e5)), sample(d3, min(length(d3), 1e5)), d4)
  dd <- d[tsid[order(tsid)],]
  icol <- 1
  dd$NP <- dd$NP+0.2*ymax
  for (i in unique(dd$CHR)) {
      with(dd[dd$CHR==i,], points(pos, NP, col=pcol[icol], cex=pcex, pch=ph))
          icol=icol+1
  }
  d1 <- which((d$OP>0)&(d$OP<=1))
  d2 <- which((d$OP>1)&(d$OP<=2))
  d3 <- which((d$OP>2)&(d$OP<=3))
  d4 <- which(d$OP>3)
  tsid <- c(sample(d1, 1e5), sample(d2, min(length(d2), 1e5)), sample(d3, min(length(d3), 1e5)), d4)
  dd <- d[tsid[order(tsid)],]
  icol <- 1
  dd$OP <- -dd$OP
  for (i in unique(dd$CHR)) {
      with(dd[dd$CHR==i,], points(pos, OP, col=pcol[icol], cex=pcex, pch=ph))
          icol=icol+1
  }

  locistatm <- matrix(rep(0, 4), 2)
  if((length(lco)!=0)&(length(lcn)!=0)) {
    for(i in 1:length(lco)) {
      if(length(intersect(osid[lco[[i]]], nsid[unlist(lcn)]))==0) {
        points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method="random")==1)], c(2, 3)]*c(1, -1)-c(0, 0.024*ymax), col="#FF6305", cex=1.1, pch=18)
        points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method="random")==1)], c(2, 4)]+c(0, 0.2*ymax), col="#FF6305", cex=2*pcex, pch=20)
        locistatm[1,1] <- locistatm[1,1]+1
      } else{
        points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method="random")==1)], c(2, 3)]*c(1, -1)-c(0, 0.017*ymax), col="#C0CD28", cex=0.8, lwd=2, pch=4)
        locistatm[2,1] <- locistatm[2,1]+1
      }
    }
    for(i in 1:length(lcn)) {
      if(length(intersect(nsid[lcn[[i]]], osid[unlist(lco)]))==0) {
        points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method="random")==1)], c(2, 4)]+c(0, 0.224*ymax), col="#FF6305", cex=1.1, pch=18)
        points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method="random")==1)], c(2, 3)]*c(1, -1), col="#FF6305", cex=2*pcex, pch=20)
        locistatm[2,2] <- locistatm[2,2]+1
      } else{
        points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method="random")==1)], c(2, 4)]+c(0, 0.217*ymax), col="#C0CD28", cex=0.8, lwd=2, pch=4)
        locistatm[1,2] <- locistatm[1,2]+1
      }
    }
  } else if(length(lco)!=0) {
    for(i in 1:length(lco)) {
      points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method="random")==1)], c(2, 3)]*c(1, -1)-c(0, 0.024*ymax), col="#FF6305", cex=1.1, pch=18)
      points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method="random")==1)], c(2, 4)]+c(0, 0.2*ymax), col="#FF6305", cex=2*pcex, pch=20)
      locistatm[1,1] <- locistatm[1,1]+1
    }
  } else if(length(lcn)!=0) {
    for(i in 1:length(lcn)) {
      points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method="random")==1)], c(2, 4)]+c(0, 0.224*ymax), col="#FF6305", cex=1.1, pch=18)
      points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method="random")==1)], c(2, 3)]*c(1, -1), col="#FF6305", cex=2*pcex, pch=20)
      locistatm[2,2] <- locistatm[2,2]+1
    }
  }

  rect(xmin, 0, xmax, 0.2*ymax, col="white", border=NA)
  text(newtick, 0.16*ymax, c(1:18, "19", 20, "21", 22), cex=0.7, xpd=NA)
  text(newtick[8], 0.07*ymax, "Chromosome", cex=1.1, xpd=NA)
  abline(h=-log10(5e-8)+0.2*ymax, lty=2)
  abline(h=-log10(fpn)+0.2*ymax, lty=1)
  abline(h=-log10(1)+0.2*ymax, lty=1)
  abline(h=log10(5e-8), lty=2)
  abline(h=log10(fpo), lty=1)
  abline(h=log10(1), lty=1)
  return(locistatm)
}

manpos <- function(idm) {
  chn <- unique(idm$CHR)
  pos <- idm$BP
  newtick <- c()
  incm <- 0
  for(i in chn) {
    temid <- which(idm$CHR==i)
    temBP <- idm$BP[temid]
    decm <- temBP[1]
    pos[temid] <- pos[temid]-decm+1+incm
    incm <- incm+temBP[length(temBP)]-temBP[1]+1e7
    newtick <- c(newtick, (pos[temid[1]]+pos[temid[length(temid)]])/2)
  }
  return(list(pos, newtick))
}
