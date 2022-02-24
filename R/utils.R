# Send R Output to a File
logOutput <- function(..., cgwasenv) {
  sink(file.path(cgwasenv$.CGWAS_RESULT_PATH, 'LogFile'), append = TRUE, split = TRUE)
  cat(paste0(...))
  sink()
}

paramOutput <- function(cgwasenv) {
  logOutput("\n========== C-GWAS basic parameters ==========\n\n",
            cgwasenv = cgwasenv)
  logOutput("GWAS result file path : ", dirname(cgwasenv$.GWAS_FILE_PATH[1]), '\n',
            cgwasenv = cgwasenv)
  logOutput("SNP information file path : ", paste(cgwasenv$.SNP_FILE_PATH, collapse = ", "), '\n',
            cgwasenv = cgwasenv)
  logOutput("GWAS name : ", paste(cgwasenv$.TRAIT_NAME, collapse = ", "), '\n',
            cgwasenv = cgwasenv)
  logOutput("Exclude NA : ", cgwasenv$.EXCLUDE_NA, '\n',
            cgwasenv = cgwasenv)
  logOutput("MAF path : ", cgwasenv$.MAF_FILE_PATH, '\n',
            cgwasenv = cgwasenv)
  logOutput("Output path : ", cgwasenv$.CGWAS_DIR, '\n',
            cgwasenv = cgwasenv)
  logOutput("Keep i-EbICoW output : ", cgwasenv$.KEEP_EbICoW, '\n',
            cgwasenv = cgwasenv)
  logOutput("Parallel number : ", cgwasenv$.PARAL_NUM, '\n',
            cgwasenv = cgwasenv)
  logOutput("\n========== C-GWAS advanced parameters ==========\n\n",
            cgwasenv = cgwasenv)
  paramAdvL <- list(cgwasenv$.IND_SNP_N,
                    cgwasenv$.SIMUL_DEP,
                    cgwasenv$.P_THRD_STUDY,
                    cgwasenv$.P_THRD_SUGST,
                    cgwasenv$.P_MAIN_EFFECT,
                    cgwasenv$.MIN_EbICo_POWER_INC,
                    cgwasenv$.MIN_CORR_DIFF,
                    cgwasenv$.HIGH_CORR2_RES,
                    cgwasenv$.SAMPLE_SIZE_INC,
                    cgwasenv$.TWT_STRAT_CUT,
                    cgwasenv$.LOESS_INTER_P,
                    cgwasenv$.LOESS_SPAN_V,
                    cgwasenv$.LOCI_INTER)
  paramAdvDefL <- list(1e6,
                       100,
                       0.05/cgwasenv$.IND_SNP_N,
                       1/cgwasenv$.IND_SNP_N,
                       3/cgwasenv$.IND_SNP_N,
                       1,
                       0.05,
                       0.5,
                       0.5,
                       10^(seq(0, 1-ceiling(log10(cgwasenv$.IND_SNP_N)), -1/3))[-1],
                       c(0.05, 0.001),
                       c(0.03, 0.1, 0.75),
                       2.5e5)
  paramAdvId <- c("Independent SNP number",
                  "Quantile simulation depth",
                  "Study-wide significant threshold",
                  "Suggestive significant threshold",
                  "Main effect threshold",
                  "Min EbICoW power increase ratio",
                  "Min |psi-pi| difference",
                  "Maximum squared psi restriction",
                  "Min Ess increase ratio",
                  "TWT stratification cutoff",
                  "LOESS interval quantile",
                  "LOESS span vector",
                  "Loci interval")
  userDefIdx <- c()
  for (i in 1:length(paramAdvL)) {
    if (!all(paramAdvL[[i]]==paramAdvDefL[[i]])) {
      userDefIdx <- c(userDefIdx, i)
    }
  }
  logOutput("User define ", length(userDefIdx), "/13 parameters\n",
            cgwasenv = cgwasenv)
  if (length(userDefIdx)!=0) {
    logOutput("\n", cgwasenv = cgwasenv)
    for (i in userDefIdx) {
      logOutput(paramAdvId[i], " : ", paste(paramAdvL[[i]], collapse = ","), '\n',
                cgwasenv = cgwasenv)
    }
    logOutput("\n", cgwasenv = cgwasenv)
  }
  logOutput("Other ", 13 - length(userDefIdx), "/13 parameters are in default\n\n\n",
            cgwasenv = cgwasenv)
}

StatE1 <- function(traitid, cgwasenv) {
  df <- data.table::fread(cgwasenv$.GWAS_FILE_PATH[traitid],
                          header = T, stringsAsFactors = F, nThread = 1)
  df <- as.data.frame(df)
  bpm <- df
  return(unique(which(is.na(bpm[,1])), which(is.na(bpm[,2]))))
}

StatE2 <- function(traitid, naidList, naid, snpPosOrder, cgwasenv) {
  bpm <- data.table::fread(cgwasenv$.GWAS_FILE_PATH[traitid],
                           header = T, stringsAsFactors = F, nThread = 1)
  bpm <- as.matrix(bpm)
  if (cgwasenv$.EXCLUDE_NA) {
    if (length(naid) != 0) {
      bpm <- bpm[-naid,]
    }
  } else {
    if (length(naidList[[traitid]]) != 0) {
      bpm <- bpm[naidList[[traitid]], 1] <- 0
      bpm <- bpm[naidList[[traitid]], 2] <- 1
    }
  }
  bpm <- bpm[snpPosOrder,]

  data.table::fwrite(as.data.frame(signif(bpm[,1], 7)),
         file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(cgwasenv$.TRAIT_NAME[traitid], ".beta")),
         row.names = F, col.names = F, quote = F, nThread = 1)
  bpm[bpm[,2] < (1e-300), 2] <- 1e-300
  data.table::fwrite(as.data.frame(signif(sqrt(qchisq(bpm[,2], 1, lower.tail = F))*sign(bpm[,1]), 7)),
         file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(cgwasenv$.TRAIT_NAME[traitid], ".stat")),
         row.names = F, col.names = F, quote = F, nThread = 1)

  return(nrow(bpm))
}

qcf <- function(lam, ppv, tq, meaq, inseq) {
  qcr <- rep(NA, length(lam))
  for (i in 1:length(qcr)) {
    if (ppv[i]==0) {
      qcr[i] <- mean((pchisq(tq/lam[i], 1)*(1-ppv[i])-inseq)^2)
    } else{
      qcr[i] <- mean((pchisq(tq/lam[i], 1)*(1-ppv[i])+
                      pchisq(tq/(lam[i]+(meaq-lam[i])/ppv[i]), 1)*ppv[i]-inseq)^2)
    }
  }
  return(qcr)
}

odgridse <- function(rb, intv, lamsp, tq, meaq, N, inseq) {
  while (intv>1e-4) {
    cid <- order(rb)[1]
    if (cid==1) {
      lamsp <- seq(lamsp[cid], lamsp[cid+1], length.out = 5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1, 5)] <- rb[c(cid, cid+1)]
      rb[3] <- qcf(lamsp[3], propv[3], tq, meaq, inseq)
      intv <- intv/4
    } else if (cid==5) {
      lamsp <- seq(lamsp[cid-1], lamsp[cid], length.out = 5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb[c(1, 5)] <- rb[c(cid-1, cid)]
      rb[3] <- qcf(lamsp[3], propv[3], tq, meaq, inseq)
      intv <- intv/4
    } else{
      lamsp <- seq(lamsp[cid-1], lamsp[cid+1], length.out = 5)
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
  if ((N<=0)|((meaq-1)<=0)) {
    lamv <- max(1, min(meaq, tlam))
    propv <- 0
  } else {
    tq <- quantile(tm2, inseq)
    lamsp <- seq(max(meaq-sqrt(N/3/(1/mp-1)), 0.99), meaq+0.01, length.out = 5)
    propv <- 1/(N/3/(meaq-lamsp)^2+1)
    rb <- qcf(lamsp, propv, tq, meaq, inseq)
    intv <- (lamsp[5]-lamsp[1])/2
    lamv <- max(c(odgridse(rb, intv, lamsp, tq, meaq, N, inseq), 1))
    propv <- 1/(N/3/(meaq-lamv)^2+1)
    if (lamv > min(meaq, tlam)) {
      lamsp <- seq(max(meaq-sqrt(N/3/(1/mp-1)), 0.99), min(meaq, tlam), length.out = 5)
      propv <- 1/(N/3/(meaq-lamsp)^2+1)
      rb <- qcf(lamsp, propv, tq, meaq, inseq)
      intv <- (lamsp[5]-lamsp[1])/2
      lamv <- odgridse(rb, intv, lamsp, tq, meaq, N, inseq)
      lamv <- max(c(odgridse(rb, intv, lamsp, tq, meaq, N, inseq), 1))
      propv <- 1/(N/3/(meaq-lamv)^2+1)
    }
  }
  return(c(lamv, propv))
}

meanByHand <- function(x) {
  return(sum(x) / length(x))
}

CorEsti.ICE <- function(t1, t2, id1, id2, resinfm, maSpltRw, maSpltN, cgwasenv) {
  mt2 <- resinfm[c(id1, id2), 4]
  cor.reg <- c()
  chisq.cv <- qchisq(0.5, 2)

  t1m <- matrix(t1[1:maSpltN], nrow = maSpltRw, byrow = F)
  t2m <- matrix(t2[1:maSpltN], nrow = maSpltRw, byrow = F)
  t1msq <- t1m ^2
  t2msq <- t2m ^2
  t1msqT2msqAdd <- t1msq + t2msq
  t1mT2mProd <- t1m * t2m
  for (i in 1:maSpltRw) {
    t1vsq = t1msq[i,]
    t2vsq = t2msq[i,]
    t1t2SqAdd = t1msqT2msqAdd[i,]
    t1t2Prod = t1mT2mProd[i,]

    cor2v <- cor2 <- 0
    idxSlct <- t1t2SqAdd < chisq.cv
    cor1v <- cor1 <- sum(t1t2Prod[idxSlct])/sqrt(sum(t1vsq[idxSlct])*sum(t2vsq[idxSlct]))
    if (abs(cor1-cor2)>1e-4) {
      cor2 <- cor1
      cor2v <- c(cor2v, cor2)
      wald.t <- (t1t2SqAdd-2*cor2*t1t2Prod)/(1-cor2^2)
      idxSlct <- wald.t < chisq.cv
      cor1 <- sum(t1t2Prod[idxSlct])/sqrt(sum(t1vsq[idxSlct])*sum(t2vsq[idxSlct]))
      while (abs(cor1-cor2)>1e-4) {
        cor1v <- c(cor1v, cor1)
        rcc <- .lm.fit(cbind(1, cor2v), cor1v)$coefficients
        cor2 <- rcc[1]/(1-rcc[2])
        cor2v <- c(cor2v, cor2)
        wald.t <- (t1t2SqAdd-2*cor2*t1t2Prod)/(1-cor2^2)
        idxSlct <- wald.t < chisq.cv
        cor1 <- sum(t1t2Prod[idxSlct])/sqrt(sum(t1vsq[idxSlct])*sum(t2vsq[idxSlct]))
        if (length(cor2v)>5) {
          cor2v <- cor2v[-1]
          cor1v <- cor1v[-1]
        }
      }
    }
    cor.reg <- c(cor.reg, cor1)
  }
  cor.mean <- meanByHand(cor.reg)

  mxx <- meanByHand(t1mT2mProd)
  efc <- mxx - cor.mean
  if (efc == 0) {
    efcr <- 0
  } else {
    if (prod(mt2-1) != 0) {
      efcr <- efc/sqrt(prod(mt2-1))
      efcr <- min(efcr, 1)
      efcr <- max(efcr, -1)
    } else {
      efcr <- sign(efc)
    }
  }

  tv <- (t1msqT2msqAdd-2*cor.mean*t1mT2mProd)/(1-cor.mean^2)
  idxSlct <- tv < qchisq(cgwasenv$.P_MAIN_EFFECT, 2, lower.tail = F)
  meff <- mt2-c(meanByHand(t1msq[idxSlct]), meanByHand(t2msq[idxSlct]))
  meff[meff < 0] <- 0
  mcef <- mxx-meanByHand(t1mT2mProd[idxSlct])
  if (mcef == 0) {
    mcefr <- 0
    if ((meff[1] == 0) & (meff[2] == 0)) {
      meff[1:2] <- 1/cgwasenv$.IND_SNP_N
    }
  } else {
    if (prod(meff) != 0) {
      mcefr <- mcef/sqrt(prod(meff))
      mcefr <- min(mcefr, 1)
      mcefr <- max(mcefr, -1)
    } else {
      mcefr <- sign(mcef)
    }
  }

  return(c(mxx/sqrt(prod(mt2)), cor.mean, efc, efcr, meff, mcef, mcefr))
}


CorE.ICE <- function(id1, id2, resinfm, maSpltRw, maSpltN, cgwasenv) {
  t1 <- as.matrix(data.table::fread(
    file.path(cgwasenv$.CGWAS_iEbICoW_PATH,
              paste0(cgwasenv$.TRAIT_NAME[id1], ".stat")),
    header = F, nThread = 1))[,1]
  t2 <- as.matrix(data.table::fread(
    file.path(cgwasenv$.CGWAS_iEbICoW_PATH,
              paste0(cgwasenv$.TRAIT_NAME[id2], ".stat")),
    header = F, nThread = 1))[,1]
  corm <- CorEsti.ICE(t1, t2, id1, id2, resinfm, maSpltRw, maSpltN, cgwasenv)
  return(corm)
}

CorEsti.ebico <- function(t1, t2, id1, id2, maSpltRw, maSpltN, cgwasenv) {
  mt2 <- as.numeric(c(id1, id2)) + 1
  cor.reg <- c()
  chisq.cv <- qchisq(0.5, 2)

  t1m <- matrix(t1[1:maSpltN], nrow = maSpltRw, byrow = F)
  t2m <- matrix(t2[1:maSpltN], nrow = maSpltRw, byrow = F)
  t1msq <- t1m ^2
  t2msq <- t2m ^2
  t1msqT2msqAdd <- t1msq + t2msq
  t1mT2mProd <- t1m * t2m
  for (i in 1:maSpltRw) {
    t1vsq = t1msq[i,]
    t2vsq = t2msq[i,]
    t1t2SqAdd = t1msqT2msqAdd[i,]
    t1t2Prod = t1mT2mProd[i,]

    cor2v <- cor2 <- 0
    idxSlct <- t1t2SqAdd < chisq.cv
    cor1v <- cor1 <- sum(t1t2Prod[idxSlct])/sqrt(sum(t1vsq[idxSlct])*sum(t2vsq[idxSlct]))
    if (abs(cor1-cor2)>1e-4) {
      cor2 <- cor1
      cor2v <- c(cor2v, cor2)
      wald.t <- (t1t2SqAdd-2*cor2*t1t2Prod)/(1-cor2^2)
      idxSlct <- wald.t < chisq.cv
      cor1 <- sum(t1t2Prod[idxSlct])/sqrt(sum(t1vsq[idxSlct])*sum(t2vsq[idxSlct]))
      while (abs(cor1-cor2)>1e-4) {
        cor1v <- c(cor1v, cor1)
        rcc <- .lm.fit(cbind(1, cor2v), cor1v)$coefficients
        cor2 <- rcc[1]/(1-rcc[2])
        cor2v <- c(cor2v, cor2)
        wald.t <- (t1t2SqAdd-2*cor2*t1t2Prod)/(1-cor2^2)
        idxSlct <- wald.t < chisq.cv
        cor1 <- sum(t1t2Prod[idxSlct])/sqrt(sum(t1vsq[idxSlct])*sum(t2vsq[idxSlct]))
        if (length(cor2v)>5) {
          cor2v <- cor2v[-1]
          cor1v <- cor1v[-1]
        }
      }
    }
    cor.reg <- c(cor.reg, cor1)
  }
  cor.mean <- meanByHand(cor.reg)

  mxx <- meanByHand(t1mT2mProd)
  efc <- mxx - cor.mean
  if (efc == 0) {
    efcr <- 0
  } else {
    if (prod(mt2-1) != 0) {
      efcr <- efc/sqrt(prod(mt2-1))
      efcr <- min(efcr, 1)
      efcr <- max(efcr, -1)
    } else {
      efcr <- sign(efc)
    }
  }

  tv <- (t1msqT2msqAdd-2*cor.mean*t1mT2mProd)/(1-cor.mean^2)
  idxSlct <- tv < qchisq(cgwasenv$.P_MAIN_EFFECT, 2, lower.tail = F)
  meff <- mt2-c(meanByHand(t1msq[idxSlct]), meanByHand(t2msq[idxSlct]))
  meff[meff < 0] <- 0
  mcef <- mxx-meanByHand(t1mT2mProd[idxSlct])
  if (mcef == 0) {
    mcefr <- 0
    if ((meff[1] == 0) & (meff[2] == 0)) {
      meff[1:2] <- 1/cgwasenv$.IND_SNP_N
    }
  } else {
    if (prod(meff) != 0) {
      mcefr <- mcef/sqrt(prod(meff))
      mcefr <- min(mcefr, 1)
      mcefr <- max(mcefr, -1)
    } else {
      mcefr <- sign(mcef)
    }
  }

  return(c(mt2-1, mxx/sqrt(prod(mt2)), cor.mean, efc, efcr, meff, mcef, mcefr))
}

CorE.ebico <- function(traitid, newt, TNm, x2m, maSpltRw, maSpltN, cgwasenv) {
  t1 <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH,
                                              paste0(TNm[traitid, 1], ".stat")),
                                    header = F, nThread = 1))[,1]
  corm <- CorEsti.ebico(t1, newt, TNm[traitid, 2], x2m, maSpltRw, maSpltN, cgwasenv)
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
  if (length(which((tm%*%cbef11)^2>thresc))>=length(which((tm%*%cbef21)^2>thresc))) {
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
  if (length(which((tm%*%cbef11)^2>thresc))>=length(which((tm%*%cbef21)^2>thresc))) {
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
  if (length(which((tm%*%coef1)^2>thresc))>=length(which((tm%*%coef2)^2>thresc))) {
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
  if (length(which((tm%*%coefs1)^2>thresc))>=length(which((tm%*%coefs2)^2>thresc))) {
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
  t2m <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(tid, ".stat")),
                                         header = F, stringsAsFactors = F, nThread = 1))[,1]
  b2m <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(tid, ".beta")),
                                         header = F, stringsAsFactors = F, nThread = 1))[,1]
  s2m <- t2m/b2m
  if (cgwasenv$.MAF_FILE_EXIST) {
    mse <- median(s2m^2/2/mafv/(1-mafv), na.rm = T)
  } else{
    mse <- median(s2m^2, na.rm = T)
  }
  return(mse)
}

gcrom <- function(rv, n) {
  corm <- diag(n)
  n <- 0
  for (i in 2:nrow(corm)) {
    corm[i:nrow(corm), i-1] <- corm[i-1, i:nrow(corm)] <- rv[(n+1):(n+nrow(corm)-i+1)]
    n <- n+nrow(corm)-i+1
  }
  return(corm)
}

ptc <- function(bgm, stm) {
  p1 <- eigen(bgm)
  ptv <- p1$values
  ptv[ptv<1/length(ptv)] <- 1/length(ptv)
  ptv <- ptv/(sum(ptv)/length(ptv))
  newcorm <- p1$vectors %*% diag(ptv) %*% t(p1$vectors)
  bgm <- t(newcorm/sqrt(diag(newcorm)))/sqrt(diag(newcorm))
  nn <- 0
  cm <- 0.01*(stm-bgm)
  while (!all(eigen(bgm)$values>0)) {
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
  for (i in 1:nrow(ACstatm)) {
    compn <- match(as.character(ACstatm[i,
                                        (4:(3+maxcn))
                                          [!is.na(ACstatm[i, 4:(3+maxcn)])]]),
                   cgwasenv$.TRAIT_NAME)
    if (length(compn)==1) {
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
  for (i in 1:nrow(tm)) {
    tv <- tm[i,]
    tv2 <- tv^2
    cm.tmp <- cm
    rid <- tv2 > cutoff.thv[1]
    ridN <- sum(rid)
    if (ridN == 0) {
      resv[i, 1:nresv] <- NA
      resv[i, fnresv] <- pchisq(max(tv2), 1, lower.tail = F)
      next
    } else if (ridN==length(tv)) {
      resv[i, 1] <- trtvf(tv, cm.tmp)
    } else {
      tv <- tv[rid]
      tv2 <- tv2[rid]
      cm.tmp <- cm.tmp[rid, rid]
      resv[i, 1] <- trtvf(tv, cm.tmp)
    }
    for (j in 2:length(cutoff.thv)) {
      rid <- tv2 > cutoff.thv[j]
      ridN <- sum(rid)
      if (ridN == 0) {
        resv[i, j:nresv] <- NA
        break
      } else if (ridN==length(tv)) {
        resv[i, j] <- resv[i, j-1]
      } else {
        tv <- tv[rid]
        tv2 <- tv2[rid]
        cm.tmp <- cm.tmp[rid, rid]
        resv[i, j] <- trtvf(tv, cm.tmp)
      }
    }
    resv[i, fnresv] <- pchisq(max(tv2), 1, lower.tail = F)
  }
  return(resv)
}

appmin <- function(pm) {
  return(do.call(pmin, c(as.data.frame(pm), na.rm = TRUE)))
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

tpcor.m.simumin <- function(m, efn) {
  for (i in 1: length(efn)) {
    v <- m[, i]
    m[which(v > 1e-12), i] <- 1 - (1 - v[which(v > 1e-12)]) ^ efn[i]
    m[which(v <= 1e-12), i] <- v[which(v <= 1e-12)] * efn[i]
  }
  m.min <- do.call(pmin, c(as.data.frame(m), na.rm = TRUE))
  return(m.min)
}

tpcor2 <- function(v, md) {
  lgv <- -log10(v)
  intp <- md[[length(md)]]
  v[lgv<=intp[1]] <- predict(md[[1]], intp[1])
  for (i in 2:length(md)) {
    v[(lgv>intp[i-1])&(lgv<=intp[i])] <- predict(md[[i-1]], lgv[(lgv>intp[i-1])&(lgv<=intp[i])])
  }
  v[lgv>intp[length(intp)]] <- predict(md[[length(md)-1]], intp[length(intp)])
  if (length(md)>2) {
    for (i in 2:(length(md)-1)) {
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
  options(warn = -1) # turn off warning temporarily
  sn <- cgwasenv$.SIMUL_N * cgwasenv$.SIMUL_SNP_N
  rsn <- sn/cgwasenv$.IND_SNP_N
  swwv <- isimup[order(isimup)]

  sid <- rsn*c(1:(cgwasenv$.IND_SNP_N-1))
  y <- swwv[sid]
  x <- sid/sn
  yy <- x/y

  rth <- seq(yy[cgwasenv$.IND_SNP_N-1], yy[cgwasenv$.IND_SNP_N/10], length.out = (1/0.025+1))
  rth <- rth[-length(rth)]
  j <- 1
  sid <- c()
  for (i in length(yy):1) {
    if (yy[i]>=rth[j]) {
      j <- j+1
      sid <- c(sid, i)
    }
    if (j>length(rth)) {
      break
    }
  }

  sid <- sort(unique(c(sid,
                       round(10^seq(0, -log10(cgwasenv$.IND_SNP_N), -0.025)*(cgwasenv$.IND_SNP_N-1)))),
              decreasing = T)
  y <- y[sid]
  x <- x[sid]
  write.table(cbind(x, y),
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, paste0("Null", nam, "correction.txt")),
              row.names = F, col.names = c("ExpectedQuantile", "ObservedP"), quote = F)

  yy <- x/y
  xx <- -log10(x)

  smn <- sum(sid<10)
  weight <- seq(0, 1, length.out = (smn+2))[2:(smn+1)]

  tv <- c()
  for (i in smn:1) {
    tv <- c(tv, sum(yy[(length(yy)-smn+1):(length(yy)-i+1)]*weight[1:(smn-i+1)])/sum(weight[1:(smn-i+1)]))
  }

  intp <- c(min(xx), -log10(cgwasenv$.LOESS_INTER_P), max(xx))

  yy[(length(yy)-length(tv)+1):length(yy)] <- tv
  md <- list()
  for (i in 1:length(cgwasenv$.LOESS_SPAN_V)) {
    md[[i]] <- loess(yy~xx, span = cgwasenv$.LOESS_SPAN_V[i])
  }
  md[[length(cgwasenv$.LOESS_SPAN_V)+1]] <- intp

  if (cgwasenv$.IND_SNP_N < 1e5) {
    ssid <- c(1:1e4, seq(1e4, cgwasenv$.IND_SNP_N, 10))
  } else {
    ssid <- c(1:1e4, seq(1e4, 1e5, 10), seq(1e5, cgwasenv$.IND_SNP_N, 100))
  }
  swwm <- matrix(isimup[sample(1:sn, sn)], ncol = rsn)
  for (i in 1:rsn) {
    swwm[,i] <- swwm[order(swwm[,i]), i]
  }
  infm <- t(apply(swwm[ssid,], 1, quantile, probs = c(0.05, 0.5, 0.95)))

  tpq <- seq(1/cgwasenv$.IND_SNP_N, 1-1/cgwasenv$.IND_SNP_N, length.out = cgwasenv$.IND_SNP_N)
  swwm <- matrix(isimup[sample(1:sn, sn)], ncol = rsn)
  for (i in 1:rsn) {
    swwm[,i] <- swwm[order(swwm[,i]), i]
  }
  ct <- tpcor2(tpq, md)
  swwm <- swwm*ct
  swwm[swwm>1] <- 1
  nulm <- t(apply(swwm[ssid,], 1, quantile, probs = c(0.05, 0.5, 0.95)))

  nulpm <- cbind(qbeta(0.05, ssid, cgwasenv$.IND_SNP_N+1-ssid),
                 qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid),
                 qbeta(0.95, ssid, cgwasenv$.IND_SNP_N+1-ssid))

  jpeg(file.path(cgwasenv$.CGWAS_DETAIL_PATH, paste0("Null", nam, "distribution.jpg")),
       width = 3000, height = 3000, res = 600)
  par(mar = c(3.5, 3.5, 1, 1))

  plot(NA,
       xlim = c(0, max(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)))),
       ylim = c(0, max(c(-log10(infm[1, 1]), -log10(nulpm[1, 1]), -log10(nulm[1, 1])))),
       xlab = expression(paste("Expected   ", -Log[10](italic(p)))),
       ylab = expression(paste("Observed   ", -Log[10](italic(p)))),
       mgp = c(2, 0.7, 0), las = 1, cex.axis = 0.85, tck = -0.015)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(infm[,1]), lwd = 1.8, col = "grey75", lty = 3)
  points(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(infm[,2]), pch = 20, col = "grey75", cex = 0.6)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(infm[,3]), lwd = 1.8, col = "grey75", lty = 3)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulpm[,1]), lwd = 1.8, col = "black", lty = 3)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulpm[,3]), lwd = 1.8, col = "black", lty = 3)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulm[,1]), lwd = 1.8, col = "#FD9001", lty = 3)
  points(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulm[,2]), pch = 20, col = "#FD9001", cex = 0.6)
  lines(-log10(qbeta(0.5, ssid, cgwasenv$.IND_SNP_N+1-ssid)), -log10(nulm[,3]), lwd = 1.8, col = "#FD9001", lty = 3)
  lines(c(-log10(qbeta(0.5, cgwasenv$.IND_SNP_N+1, 1)), -log10(qbeta(0.5, 1, cgwasenv$.IND_SNP_N+1))),
        c(-log10(qbeta(0.5, cgwasenv$.IND_SNP_N+1, 1)), -log10(qbeta(0.5, 1, cgwasenv$.IND_SNP_N+1))),
        col = "black")
  dev.off()

  options(warn = -1) # turn off warning
  return(md)
}

locisearch <- function(sp, spos, fp, cgwasenv) {
  n <- 0
  lclist <- list()
  if (length(sp)!=0) {
    keyid <- order(sp)[1]
    while (sp[keyid]<=fp) {
      idirc <- ddirc <- spos[keyid, 2]
      oldddirc <- oldidirc <- -1
      while ((oldddirc!=ddirc)|(oldidirc!=idirc)) {
        oldddirc <- ddirc
        oldidirc <- idirc
        tdl <- which ((spos[,1]==spos[keyid, 1]) &
                        (spos[,2]>(ddirc-cgwasenv$.LOCI_INTER)) &
                        (spos[,2]<(idirc+cgwasenv$.LOCI_INTER)))
        ddirc <- spos[min(tdl), 2]
        idirc <- spos[max(tdl), 2]
      }
      n <- n+1
      lclist[[n]] <- tdl
      spos[tdl,] <- NA
      sp[tdl] <- NA
      keyid <- order(sp)[1]
      if (is.na(sp[keyid])) {
        break
      }
    }
  }
  return(lclist)
}

manhattan <- function(op, np, osid, nsid, lco, lcn, fpo, fpn, gpo, gpn,
                      Sind, newtick,
                      pcol = c("#009DE8", "#233FAA"), pcex = 0.35, ph = 20) {
  d <- cbind(Sind[, c(1, 4)], -log10(op), -log10(np))
  colnames(d) <- c("CHR", "pos", "OP", "NP")

  ytop <- ceiling(max(d[,3:4]))
  if (ytop>30) {
    d$OP[d$OP>30] = d$OP[d$OP>30]/2+15
    d$NP[d$NP>30] = d$NP[d$NP>30]/2+15
    if (ytop>70) {
      d$OP[d$OP>50] = d$OP[d$OP>50]/1.5+50/3
      d$NP[d$NP>50] = d$NP[d$NP>50]/1.5+50/3
      if (ytop>100) {
        d$OP[d$OP>60] = d$OP[d$OP>60]*0.3+42
        d$NP[d$NP>60] = d$NP[d$NP>60]*0.3+42
        if (ytop>300) {
          d$OP[d$OP>80] = 80
          d$NP[d$NP>80] = 80
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
    plot(0, col = F, xaxt = 'n', yaxt = 'n', bty = 'n', xaxs = 'i',
         xlim = c(xmin, xmax), ylim = c(-ymax, 1.2*ymax),
         xlab = "", ylab = expression(-log[10](italic(p))),
         mgp = c(2, 0.7, 0), cex.lab = 1.1)

  axis(2, las = 1, at = l+0.2*ymax,
       labels = as.character(ll), tck = -0.02, mgp = c(2, 0.7, 0),
       lwd.ticks = 1, cex.axis = 0.85)
  axis(2, las = 1, at = -l,
       labels = as.character(ll), tck = -0.02, mgp = c(2, 0.7, 0),
       lwd.ticks = 1, cex.axis = 0.85)

  pcol <- rep(pcol, max(d$CHR))[1:max(d$CHR)]

  d1 <- which((d$NP>0)&(d$NP<=1))
  d2 <- which((d$NP>1)&(d$NP<=2))
  d3 <- which((d$NP>2)&(d$NP<=3))
  d4 <- which(d$NP>3)
  tsid <- c(sample(d1, min(length(d1), 1e5)),
            sample(d2, min(length(d2), 1e5)),
            sample(d3, min(length(d3), 1e5)), d4)
  dd <- d[tsid[order(tsid)],]
  icol <- 1
  dd$NP <- dd$NP+0.2*ymax
  for (i in unique(dd$CHR)) {
      with(dd[dd$CHR==i,], points(pos, NP, col = pcol[icol], cex = pcex, pch = ph))
          icol = icol+1
  }
  d1 <- which((d$OP>0)&(d$OP<=1))
  d2 <- which((d$OP>1)&(d$OP<=2))
  d3 <- which((d$OP>2)&(d$OP<=3))
  d4 <- which(d$OP>3)
  tsid <- c(sample(d1, min(length(d1), 1e5)), sample(d2, min(length(d2), 1e5)), sample(d3, min(length(d3), 1e5)), d4)
  dd <- d[tsid[order(tsid)],]
  icol <- 1
  dd$OP <- -dd$OP
  for (i in unique(dd$CHR)) {
      with(dd[dd$CHR==i,], points(pos, OP, col = pcol[icol], cex = pcex, pch = ph))
          icol = icol+1
  }

  locistatm <- matrix(rep(0, 4), 2)
  if ((length(lco)!=0)&(length(lcn)!=0)) {
    for (i in 1:length(lco)) {
      if (length(intersect(osid[lco[[i]]], nsid[unlist(lcn)]))==0) {
        points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method = "random")==1)], c(2, 3)]*c(1, -1)-c(0, 0.024*ymax),
               col = "#FF6305", cex = 1.1, pch = 18)
        points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method = "random")==1)], c(2, 4)]+c(0, 0.2*ymax),
               col = "#FF6305", cex = 2*pcex, pch = 20)
        locistatm[1, 1] <- locistatm[1, 1]+1
      } else{
        points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method = "random")==1)], c(2, 3)]*c(1, -1)-c(0, 0.017*ymax),
               col = "#C0CD28", cex = 0.8, lwd = 2, pch = 4)
        locistatm[2, 1] <- locistatm[2, 1]+1
      }
    }
    for (i in 1:length(lcn)) {
      if (length(intersect(nsid[lcn[[i]]], osid[unlist(lco)]))==0) {
        points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method = "random")==1)], c(2, 4)]+c(0, 0.224*ymax),
               col = "#FF6305", cex = 1.1, pch = 18)
        points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method = "random")==1)], c(2, 3)]*c(1, -1),
               col = "#FF6305", cex = 2*pcex, pch = 20)
        locistatm[2, 2] <- locistatm[2, 2]+1
      } else{
        points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method = "random")==1)], c(2, 4)]+c(0, 0.217*ymax),
               col = "#C0CD28", cex = 0.8, lwd = 2, pch = 4)
        locistatm[1, 2] <- locistatm[1, 2]+1
      }
    }
  } else if (length(lco)!=0) {
    for (i in 1:length(lco)) {
      points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method = "random")==1)], c(2, 3)]*c(1, -1)-c(0, 0.024*ymax),
             col = "#FF6305", cex = 1.1, pch = 18)
      points(d[osid[lco[[i]]][which(rank(op[osid[lco[[i]]]], ties.method = "random")==1)], c(2, 4)]+c(0, 0.2*ymax),
             col = "#FF6305", cex = 2*pcex, pch = 20)
      locistatm[1, 1] <- locistatm[1, 1]+1
    }
  } else if (length(lcn)!=0) {
    for (i in 1:length(lcn)) {
      points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method = "random")==1)], c(2, 4)]+c(0, 0.224*ymax),
             col = "#FF6305", cex = 1.1, pch = 18)
      points(d[nsid[lcn[[i]]][which(rank(np[nsid[lcn[[i]]]], ties.method = "random")==1)], c(2, 3)]*c(1, -1),
             col = "#FF6305", cex = 2*pcex, pch = 20)
      locistatm[2, 2] <- locistatm[2, 2]+1
    }
  }

  rect(xmin, 0, xmax, 0.2*ymax, col = "white", border = NA)
  text(newtick, 0.16*ymax, unique(d$CHR), cex = 0.7, xpd = NA)
  text((newtick[1]+newtick[length(newtick)])/2, 0.07*ymax, "Chromosome", cex = 1.1, xpd = NA)
  abline(h = -log10(gpn)+0.2*ymax, lty = 2)
  abline(h = -log10(fpn)+0.2*ymax, lty = 1)
  abline(h = -log10(1)+0.2*ymax, lty = 1)
  abline(h = log10(gpo), lty = 2)
  abline(h = log10(fpo), lty = 1)
  abline(h = log10(1), lty = 1)
  return(locistatm)
}

manpos <- function(idm) {
  chn <- unique(idm[,1])
  pos <- idm[,2]
  newtick <- c()
  incm <- 0
  for (i in chn) {
    temid <- which(idm[,1]==i)
    temBP <- idm[temid,2]
    decm <- temBP[1]
    pos[temid] <- pos[temid]-decm+1+incm
    incm <- incm+temBP[length(temBP)]-temBP[1]+1e7
    newtick <- c(newtick, (pos[temid[1]]+pos[temid[length(temid)]])/2)
  }
  return(list(pos, newtick))
}
