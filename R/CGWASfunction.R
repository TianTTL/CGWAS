# Data formating
step1 <- function(cgwasenv) {
  logOutput("========== C-GWAS step 1 : Data formating ==========\n\n", cgwasenv = cgwasenv)
  logOutput("Read and format ",length(cgwasenv$.TRAIT_NAME)," GWASs ..\n", cgwasenv = cgwasenv)

  if (cgwasenv$.TRAIT_NUM <= 1)  {
    stop("Error: Traits number must larger than 1.")
  } # check equality among all element in snp.N

  threadNCur <- min(cgwasenv$.PARAL_NUM, cgwasenv$.TRAIT_NUM + 1)
  cl <- makeCluster(threadNCur)
  registerDoParallel(cl)
  # globalVariables(c("i"))
  i <- 1 # assign parallel control variants

  # count NA value
  naidList <- foreach(i = seq_len(cgwasenv$.TRAIT_NUM), .inorder=F) %dopar%
                StatE1(i, cgwasenv)
  naid <- unique(unlist(naidList))
  naid.Len <- length(naid)
  if(naid.Len!=0){
    if(cgwasenv$.EXCLUDE_NA) {
      logOutput(naid.Len," SNPs with NA excluded\n", cgwasenv = cgwasenv)
    } else {
      logOutput(naid.Len," SNPs with NA set in 0 effect\n", cgwasenv = cgwasenv)
    }
  } else {
    logOutput("0 SNP with NA found\n", cgwasenv = cgwasenv)
  }

  if (cgwasenv$.MAF_FILE_EXIST) {
    mafv <- as.data.frame(data.table::fread(cgwasenv$.MAF_FILE_PATH,
                                            header=T, stringsAsFactors=F,
                                            nThread = 1))[,1]
    if (cgwasenv$.EXCLUDE_NA & naid.Len != 0) {
      mafv <- mafv[-naid]
    }
    write.table(mafv, file.path(cgwasenv$.CGWAS_iEbICoW_PATH, "MAF"), row.names=F, col.names=F, quote=F)
  }
  snp.N <- foreach(i = seq_len(cgwasenv$.TRAIT_NUM), .combine="c", .inorder = F) %dopar%
    StatE2(i, naidList, naid, cgwasenv)

  # check equality among all element in snp.N
  if (cgwasenv$.MAF_FILE_EXIST) {snp.N <- c(snp.N, length(mafv))}
  if (length(unique(snp.N)) > 1)  {
    stop("Error: SNP number in each files are not equal.")
  }
  cgwasenv$.SNP_N <- as.integer(unique(snp.N))

  stopCluster(cl)

  logOutput("\nC-GWAS step 1 completed\n", cgwasenv = cgwasenv)
}

# GetI
step2 <- function(cgwasenv) {
  logOutput("========== C-GWAS step 2 : GetI ==========\n\n", cgwasenv = cgwasenv)
  logOutput("Estimating inflation for ", cgwasenv$.TRAIT_NUM, " GWASs ..\n", cgwasenv = cgwasenv)

  inseq <- ppoints(2000)

  minsnpn = 1e5
  maSpltRw = ceiling(cgwasenv$.SNP_N / minsnpn)
  maSpltN = floor(cgwasenv$.SNP_N / maSpltRw) * maSpltRw

  threadNCur <- min(cgwasenv$.PARAL_NUM, maSpltRw)
  cl <- makeCluster(threadNCur)
  registerDoParallel(cl)
  # globalVariables(c('i', 'j'))
  i <- j <- 1 # assign parallel control variants

  resinfm <- matrix(NA, cgwasenv$.TRAIT_NUM, 4)

  for(i in seq_len(cgwasenv$.TRAIT_NUM)) {
    tm <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(cgwasenv$.TRAIT_NAME[i], ".stat")),
                                      header = F,
                                      nThread = cgwasenv$.PARAL_NUM))[,1]
    tm2 <- tm^2
    resinfm[i,1] <- mean(tm2)
    resinfm[i,2] <- median(tm2)/qchisq(0.5, 1)
    mp <- 0.1
    rantm <- matrix(tm2[1:maSpltN], nrow = maSpltRw, byrow = F)
    infm <- foreach(j=iter(rantm, by="row", chunksize=1),.combine="rbind",.inorder=T) %dopar%
              calinf(j, inseq, mp)
    if (maSpltRw == 1) { infm <- t(as.matrix(infm)) }
    infm[infm[,2] > mp, 2] <- mp
    while((sum(infm[,2] == mp) > (maSpltRw*0.1)) & (mp!=0.5)) {
      mp <- mp+0.1
      infm <- foreach(j=iter(rantm, by="row", chunksize=1),.combine="rbind",.inorder=T) %dopar%
                calinf(j, inseq, mp)
      if (maSpltRw == 1) { infm <- t(as.matrix(infm)) }
      infm[infm[,2] > mp, 2] <- mp
    }
    if(mp!=0.5) {
      resinfm[i, 3] <- mean(infm[infm[, 2]!=mp, 1])
    } else{
      resinfm[i, 3] <- mean(infm[, 1])
    }
    if((resinfm[i, 1]-resinfm[i, 3]) < 1e-3) {
      resinfm[i, 4] <- 1.001
    } else{
      resinfm[i, 4] <- 1+resinfm[i, 1]-resinfm[i, 3]
    }
    data.table::fwrite(as.data.frame(signif(tm/sqrt(resinfm[i, 1]/resinfm[i, 4]), 7)),
                       file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(cgwasenv$.TRAIT_NAME[i], ".stat")),
                       row.names = F, col.names = F, quote = F)
  }
  logOutput("Adjusted GWASs test statistics written to Details/i-EbICoW/\n", cgwasenv = cgwasenv)

  colnames(resinfm) <- c("RawMeanX2","GClambda","EstInf","AdjMeanX2")
  write.table(cbind(TraitName = cgwasenv$.TRAIT_NAME, signif(resinfm, 7)),
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "SummaryGetI.txt"),
              col.names = T, row.names = F, quote = F, sep = "\t")
  logOutput("Summary of GetI written to Details/SummaryGetI.txt\n", cgwasenv = cgwasenv)

  stopCluster(cl)

  logOutput("\nC-GWAS step 2 completed\n", cgwasenv = cgwasenv)
}

# GetPsi
step3 <- function(cgwasenv) {
  logOutput("========== C-GWAS step 3 : GetPsi ==========\n\n", cgwasenv = cgwasenv)
  logOutput("Estimating background correlation for ",
            nrow(pairma <- t(combn(seq_len(cgwasenv$.TRAIT_NUM), 2))),
            " GWAS pairs ..\n", cgwasenv = cgwasenv)

  minsnpn = 1e5
  maSpltRw = ceiling(cgwasenv$.SNP_N / minsnpn)
  maSpltN = floor(cgwasenv$.SNP_N / maSpltRw) * maSpltRw
  resinfm <- as.matrix(
               read.table(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "SummaryGetI.txt"),
                          header = T)[,-1])

  threadNCur <- min(cgwasenv$.PARAL_NUM, nrow(pairma))
  cl <- makeCluster(threadNCur)
  # globalVariables(c('i'))
  i <- 1 # assign parallel control variants

  tresm <- foreach(i = 1:nrow(pairma),
                   .combine = "rbind",
                   .inorder = T) %dopar%
    CorE.ICE(pairma[i, 1], pairma[i, 2], resinfm, maSpltRw, maSpltN, cgwasenv)
  tresm <- signif(tresm, 7)

  corm <- cbind(cgwasenv$.TRAIT_NAME[pairma[,1]],
                cgwasenv$.TRAIT_NAME[pairma[,2]],
                signif(resinfm[pairma[,1], 4]-1, 7),
                signif(resinfm[pairma[,2], 4]-1, 7),
                tresm)

  colnames(corm) <- c("GWAS1","GWAS2",
                      "T1Eff","T2Eff",
                      "StatCor","Psi","EffCov","allPi","T1sEff","T2sEff","EffsCov","sigPi")
  write.table(corm,
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "BCCorrelationStat.txt"),
              row.names = F, quote = F, sep = "\t")
  write.table(corm[,c(1,2,5,6,8,12)],
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "SummaryGetPsi.txt"),
              row.names = F, quote = F, sep = "\t")
  logOutput("Summary of GetPsi written to Details/SummaryGetPsi.txt\n", cgwasenv = cgwasenv)

  stopCluster(cl)

  logOutput("\nC-GWAS step 3 completed\n", cgwasenv = cgwasenv)
}

# EbICo
step4 <- function(cgwasenv) {
  logOutput("========== C-GWAS step 4 : i-EbICoW ==========\n\n", cgwasenv = cgwasenv)
  logOutput("Iteration preparing ..\n", cgwasenv = cgwasenv)

  threadNCur <- min(cgwasenv$.PARAL_NUM, cgwasenv$.TRAIT_NUM)
  cl <- makeCluster(threadNCur)
  registerDoParallel(cl)
  # globalVariables(c('i'))
  i <- 1 # assign parallel control variants

  minsnpn = 1e5
  maSpltRw = ceiling(cgwasenv$.SNP_N / minsnpn)
  maSpltN = floor(cgwasenv$.SNP_N / maSpltRw) * maSpltRw

  thresc <- qchisq(cgwasenv$.P_MAIN_EFFECT, 1, lower.tail=F)
  thresw <- qchisq(cgwasenv$.P_MAIN_EFFECT, 2, lower.tail=F)

  if (cgwasenv$.MAF_FILE_EXIST) {
    mafv <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, "MAF"),
                                            header=F, stringsAsFactors=F,
                                            nThread = cgwasenv$.PARAL_NUM))[,1]
  }

  corm.tmp <- read.table(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "BCCorrelationStat.txt"),
                       header=T, stringsAsFactors=F)
  coridm <- as.matrix(corm.tmp[,1:2])
  cordatm <- as.matrix(corm.tmp[,3:12])

  n <- cgwasenv$.TRAIT_NUM+1
  AllEff <- as.numeric(c(cordatm[-(1:(n-2))^2/2+(n-0.5)*(1:(n-2))-n+2, 1],
                         cordatm[-(n-2)^2/2+(n-0.5)*(n-2)-n+2, 2]))
  Ess <- unlist(foreach(i = cgwasenv$.TRAIT_NAME, .inorder=T) %dopar%
                  Essfun(i, mafv, cgwasenv))
  TNm <- cbind(data.frame(TraitName = cgwasenv$.TRAIT_NAME, stringsAsFactors=F), AllEff, Ess)
  idlist <- as.list(cgwasenv$.TRAIT_NAME)
  bflist <- as.list(rep(1, cgwasenv$.TRAIT_NUM))
  tflist <- as.list(rep(1, cgwasenv$.TRAIT_NUM))
  rnlist <- as.list(rep(NA, cgwasenv$.TRAIT_NUM))

  bkcordatm <- c()
  bknum <- 1
  logfm <- c()
  rn <- 1

  logOutput("Iteration start ..\n\n", cgwasenv = cgwasenv)

  while(!(all(abs(cordatm[,6] - cordatm[,4]) + abs(cordatm[,10] - cordatm[,4]) < cgwasenv$.MIN_CORR_DIFF) &
          all(abs(cordatm[,4]) < cgwasenv$.HIGH_CORR_RES))) {
    mtarid <- order(abs(cordatm[,6]-cordatm[,4]+cordatm[,10]-cordatm[,4]), decreasing=T)[1]
    selidv <- coridm[mtarid,]
    newname <- paste0(selidv[1], "_", selidv[2])
    t2m <- cbind(
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[1], ".stat")),
                                         header=F, stringsAsFactors=F,
                                         nThread = cgwasenv$.PARAL_NUM))[,1],
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[2], ".stat")),
                                         header=F, stringsAsFactors=F,
                                         nThread = cgwasenv$.PARAL_NUM))[,1]
           )
    mse <- c(TNm[which(TNm[,1]==selidv[1]), 3], TNm[which(TNm[,1]==selidv[2]), 3])

    ttm <- t2m[which((((t2m %*%
                        solve(matrix(c(1, cordatm[mtarid, 4], cordatm[mtarid, 4], 1), 2))) * t2m) %*%
                        c(1, 1))
                     > qchisq(0.01, 2, lower.tail=F)),]
    fm <- ebicocof(cordatm[mtarid,], mse, ttm, thresc)
    cofm <- fm[,1:3]
    temsm <- (ttm %*% cofm)^2
    temwm <- ((ttm %*% solve(matrix(c(1, cordatm[mtarid, 4], cordatm[mtarid, 4], 1), 2)))*ttm) %*% c(1, 1)

    efin <- 2 -
            0.0725615289639716 * abs(cordatm[mtarid, 4]) -
            0.196105591438997 * cordatm[mtarid, 4] ^ 2 -
            0.192036017551772 * cordatm[mtarid, 4] ^ 4 -
            0.134862184571824 * cordatm[mtarid, 4] ^ 8 -
            0.0969596906274313 * cordatm[mtarid, 4] ^ 16 -
            0.260927307332186 * cordatm[mtarid, 4] ^ 64 +
            0.232255565305109 * cordatm[mtarid, 4] ^ 128 -
            0.27795687154795 * cordatm[mtarid, 4] ^ 256
    ge1 <- which((ttm[,1]^2) > qchisq(cgwasenv$.P_MAIN_EFFECT / efin, 1, lower.tail=F))
    ge2 <- which((ttm[,2]^2) > qchisq(cgwasenv$.P_MAIN_EFFECT / efin, 1, lower.tail=F))
    ge1c <- which((ttm[,1]^2) > thresc)
    ge2c <- which((ttm[,2]^2) > thresc)
    idg <- union(ge1, ge2)
    idgc <- union(ge1c, ge2c)
    idog <- intersect(ge1c, ge2c)
    ide <- which(temsm[,1]>thresc)
    ids <- which(temsm[,2]>thresc)
    ida <- which(temsm[,3]>thresc)
    idw <- which(temwm>thresw)
    gen1 <- length(ge1c)
    gen2 <- length(ge2c)
    gen <- max(c(length(idg), gen1, gen2))
    een <- length(ide)
    esn <- length(ids)
    ean <- length(ida)
    ewn <- length(idw)
    jee <- length(intersect(idgc, ide))
    jes <- length(intersect(idgc, ids))
    jea <- length(intersect(idgc, ida))
    jeeo <- length(intersect(idog, ide))
    jeso <- length(intersect(idog, ids))
    jeao <- length(intersect(idog, ida))

    logOutput(rn, " Round    ", sum((abs(cordatm[,6]-cordatm[,4])+abs(cordatm[,10]-cordatm[,4])>cgwasenv$.MIN_CORR_DIFF)|(abs(cordatm[,4])>cgwasenv$.HIGH_CORR_RES))-1, " Left\n", cgwasenv = cgwasenv)
    logOutput(substring(selidv[1],1,70),"  ",signif(cordatm[mtarid,1],3),"  ",gen1, '\n', cgwasenv = cgwasenv)
    logOutput(substring(selidv[2],1,70),"  ",signif(cordatm[mtarid,2],3),"  ",gen2, '\n', cgwasenv = cgwasenv)
    logOutput(substring(newname,1,70),"  ",signif(cordatm[mtarid,4],3), '\n', cgwasenv = cgwasenv)
    logOutput("using AEC  (",signif(cordatm[mtarid,6],3),")  ",signif(cofm[1,1],3),"  ",signif(cofm[2,1],3)," -> ",een,"/",gen,"(",signif(een/gen*100,3),"%)\n", cgwasenv = cgwasenv)
    logOutput("using SEC  (",signif(cordatm[mtarid,10],3),")  ",signif(cofm[1,2],3),"  ",signif(cofm[2,2],3)," -> ",esn,"/",gen,"(",signif(esn/gen*100,3),"%)\n", cgwasenv = cgwasenv)
    logOutput("using STB  (",signif(fm[2,9],3),")  ",signif(cofm[1,3],3),"  ",signif(cofm[2,3],3)," -> ",ean,"/",gen,"(",signif(ean/gen*100,3),"%)\n", cgwasenv = cgwasenv)
    logOutput("using WDT -> ", ewn, "/", gen, "(", signif(ewn/gen*100, 3), "%)\n", cgwasenv = cgwasenv)

    rn <- rn+1

    if ((een>=esn)&(jeeo>=jeso)&(jee>=jes)) {
      if((een>=ean)&(jeeo>=jeao)&(jee>=jea)) {
        did <- 1
      } else if((ean>=een)&(jeao>=jeeo)&(jea>=jee)) {
        did <- 3
      } else{
        did <- which(c(een, esn, ean)==max(c(een, esn, ean)))[1]
      }
    } else if((esn>=ean)&(jeso>=jeao)&(jes>=jea)) {
      did <- 2
    } else if((ean>=esn)&(jeao>=jeso)&(jea>=jes)) {
      did <- 3
    } else{
      did <- which(c(een, esn, ean)==max(c(een, esn, ean)))[1]
    }

    if ((c(een, esn, ean)[did] < max(c(gen*cgwasenv$.MIN_EbICo_POWER_INC, ewn*cgwasenv$.MIN_EbICo_POWER_INC))) &
        (abs(cordatm[mtarid, 4])<cgwasenv$.HIGH_CORR_RES)) {
      logOutput("Decision -> do not combine\n\n", cgwasenv = cgwasenv)
      bkcordatm <- rbind(bkcordatm, cordatm[mtarid,])
      logfm <- rbind(logfm,
                     c(rn-1, selidv, newname, NA, FALSE,
                       signif(c(cordatm[mtarid,c(1,2)], NA,
                                cordatm[mtarid,c(4,6,10)],
                                gen1,gen2,gen,een,esn,ean,ewn,jee,jes,jea,jeeo,jeso,jeao),
                              7),
                       rep(NA, 7)))
      cordatm[mtarid,4] <- cordatm[mtarid,6] <- cordatm[mtarid,10] <- 0
      cordatm[mtarid,3] <- bknum+1
      bknum <- bknum+1
      next
    }

    newt <- t2m %*% cofm[,did]
    x2m <- mean(newt^2)-1
    if (x2m < 0.001) {
      newt <- newt*sqrt(1.001/(x2m+1))
      x2m <- mean(newt^2)-1
    }

    bofm <- fm[,4:6]
    dofm <- fm[,7:9]

    b2m <- cbind(
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[1], ".beta")),
                                         header=F, stringsAsFactors=F,
                                         nThread = cgwasenv$.PARAL_NUM))[,1],
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[2], ".beta")),
                                         header=F, stringsAsFactors=F,
                                         nThread = cgwasenv$.PARAL_NUM))[,1]
           )
    bofv <- bofm[,did]*sqrt(mse)/sum((bofm[,did]*dofm[,did])*sqrt(mse))
    newb <- b2m %*% bofv

    if(cgwasenv$.MAF_FILE_EXIST) {
      newmse <- median((newt/newb)^2/2/mafv/(1-mafv), na.rm=T)
    } else{
      newmse <- median((newt/newb)^2, na.rm=T)
    }

    logOutput("Equivalent sample size  ", round(mse[1]), " ", round(mse[2]), " -> ", round(newmse), '\n', cgwasenv = cgwasenv)
    logOutput("Mean beta weights  ", signif(bofv[1], 3), " ", signif(bofv[2], 3), "   total ", sum((bofv*dofm[,did])), '\n', cgwasenv = cgwasenv)
    if((newmse<max(c(mse[1]+mse[2]*cgwasenv$.SAMPLE_SIZE_INC, mse[2]+mse[1]*cgwasenv$.SAMPLE_SIZE_INC))) &
       (abs(cordatm[mtarid, 4])<cgwasenv$.HIGH_CORR_RES)) {
      logOutput("Decision -> do not combine\n\n", cgwasenv = cgwasenv)
      bkcordatm <- rbind(bkcordatm, cordatm[mtarid,])
      logfm <- rbind(logfm,
                     c(rn-1, selidv, newname, c("AEC","SEC","STB")[did], FALSE,
                       signif(c(cordatm[mtarid,c(1,2)], x2m, cordatm[mtarid,c(4,6,10)],
                                gen1,gen2,gen,een,esn,ean,ewn,jee,jes,jea,jeeo,jeso,jeao),
                              7),
                       round(c(mse,newmse)),
                       signif(bofv,7),
                       signif(cofm[,did],7)))
      cordatm[mtarid,4] <- cordatm[mtarid,6] <- cordatm[mtarid,10] <- 0
      cordatm[mtarid,3] <- bknum+1
      bknum <- bknum+1
      next
    } else{
      logOutput("Decision -> use ", c("AEC", "SEC", "STB")[did], '\n\n', cgwasenv = cgwasenv)
    }

    logfm <- rbind(logfm,
                   c(rn-1, selidv, newname, c("AEC","SEC","STB")[did], TRUE,
                     signif(c(cordatm[mtarid,c(1,2)], x2m, cordatm[mtarid,c(4,6,10)],
                              gen1,gen2,gen,een,esn,ean,ewn,jee,jes,jea,jeeo,jeso,jeao),
                            7),
                     round(c(mse,newmse)),
                     signif(bofv,7),
                     signif(cofm[,did],7)))

    if (nrow(coridm) < 2) break

    selid1 <- which(TNm[,1]==selidv[1])
    selid2 <- which(TNm[,1]==selidv[2])
    TNm[selid1,] <- NA
    TNm[selid2,] <- NA
    newcordatm <- foreach(i=which(!is.na(TNm[,1])), .combine="rbind", .inorder=T) %dopar%
                    CorE.ebico(i, newt, TNm, x2m, maSpltRw, maSpltN, cgwasenv)
    newcoridm <- cbind(na.omit(TNm)[,1], newname)
    deleteid <- union(union(union(
                                  which(coridm[,1]==selidv[1]), which(coridm[,1]==selidv[2])),
                            which(coridm[,2]==selidv[1])),
                      which(coridm[,2]==selidv[2]))
    coridm <- coridm[-deleteid,]
    cordatm <- cordatm[-deleteid,]

    coridm <- rbind(coridm, newcoridm)
    cordatm <- rbind(cordatm, newcordatm)
    TNm <- rbind(TNm, TNm[nrow(TNm),])
    TNm[nrow(TNm),1] <- newname
    TNm[nrow(TNm),2] <- x2m
    TNm[nrow(TNm),3] <- newmse
    idlist[[n]] <- c(idlist[[selid1]], idlist[[selid2]])
    idlist[[selid1]] <- idlist[[selid2]] <- NA
    bflist[[n]] <- c(bflist[[selid1]]*bofv[1], bflist[[selid2]]*bofv[2])
    bflist[[selid1]] <- bflist[[selid2]] <- NA
    tflist[[n]] <- c(tflist[[selid1]]*cofm[1, did], tflist[[selid2]]*cofm[2, did])
    tflist[[selid1]] <- tflist[[selid2]] <- NA
    rnlist[[n]] <- c(rnlist[[selid1]], rnlist[[selid2]], rn-1)
    rnlist[[n]] <- rnlist[[n]][!is.na(rnlist[[n]])]
    rnlist[[selid1]] <- rnlist[[selid2]] <- NA
    n <- n+1

    if (selid1 > cgwasenv$.TRAIT_NUM) {
      file.remove(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[1], ".stat")))
      file.remove(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[1], ".beta")))
    }
    if (selid2 > cgwasenv$.TRAIT_NUM) {
      file.remove(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[2], ".stat")))
      file.remove(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(selidv[2], ".beta")))
    }

    data.table::fwrite(as.data.frame(signif(newt, 7)),
                       file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(newname, ".stat")),
                       row.names=F,col.names=F,quote=F)
    data.table::fwrite(as.data.frame(signif(newb, 7)),
                       file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(newname, ".beta")),
                       row.names=F,col.names=F,quote=F)
  }

  logOutput("Iteration end\n\n", cgwasenv = cgwasenv)

  compm <- bfm <- tfm <- cdm <- rnm <- as.data.frame(rep(NA, sum(!is.na(TNm[,1]))))
  tn <- 1
  for(i in 1:(n-1)) {
    if(!is.na(idlist[[i]][1])) {
      if(length(idlist[[i]])>ncol(compm)) {
        compm <- cbind(compm, t(rep(NA, length(idlist[[i]])-ncol(compm))))
        bfm <- cbind(bfm, t(rep(NA, length(bflist[[i]])-ncol(bfm))))
        tfm <- cbind(tfm, t(rep(NA, length(tflist[[i]])-ncol(tfm))))
        rnm <- cbind(rnm, t(rep(NA, length(rnlist[[i]])-ncol(rnm))))
      }
      compm[tn,1:length(idlist[[i]])] <- idlist[[i]]
      bfm[tn,1:length(bflist[[i]])] <- bflist[[i]]
      tfm[tn,1:length(tflist[[i]])] <- tflist[[i]]
      rnm[tn,1:length(rnlist[[i]])] <- rnlist[[i]]
      tn <- tn+1
    }
  }

  colnames(compm) <- paste0("GWAS",1:ncol(compm))
  colnames(bfm) <- paste0("bw",1:ncol(bfm))
  colnames(tfm) <- paste0("tw",1:ncol(tfm))
  colnames(rnm) <- paste0("r", 1:ncol(rnm))

  TNm <- na.omit(TNm)
  TNm[,2] <- signif(TNm[,2], 7)

  for(i in 1:nrow(cordatm)) {
    if(cordatm[i, 3]>=2) {
      cordatm[i,] <- bkcordatm[cordatm[i,3]-1,]
    }
  }

  orderid <- c()
  for(i in 1:(nrow(TNm)-1)) {
    for(j in (i+1):nrow(TNm)) {
      orderid <- c(orderid, which((coridm[,1]==TNm[i, 1])&(coridm[,2]==TNm[j, 1])))
    }
  }

  colnames(logfm) <- c("Round","GWAS1","GWAS2","NewGWAS","CovType","Evaluate",
                       "Eff1X2","Eff2X2","NewEffX2","Psi","allPi","sigPi",
                       "sigSNP1N","sigSNP2N","sigSNPminN","sigSNPAECN",
                       "sigSNPSECN","sigSNPSTBN","sigSNPPWDN","sigSNPAECuniN",
                       "sigSNPSECuniN","sigSNPSTBuniN","sigSNPAECinterN",
                       "sigSNPSECinterN","sigSNPSTBinterN",
                       "Ess1","Ess2","EssNew","bw1","bw2","tw1","tw2")

  write.table(cbind(TNm, compm, signif(bfm, 7), signif(tfm, 7), rnm),
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "EbICoW.txt"),
              row.names=F,quote=F,sep="\t")
  logOutput("Summary of EbICoW GWASs written to Details/EbICoW.txt\n", cgwasenv = cgwasenv)

  accorm <- cbind(matrix(coridm[orderid,], ncol=ncol(coridm)), matrix(signif(cordatm[orderid,], 7), ncol=ncol(cordatm)))
  colnames(accorm) <- c("EbICoWGWAS1","EbICoWGWAS2","T1Eff","T2Eff",
                        "StatCor","Psi","EffCov","allPi","T1sEff","T2sEff",
                        "EffsCov","sigPi")
  write.table(accorm,
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "ACCorrelationStat.txt"),
              row.names=F,quote=F,sep="\t")
  write.table(accorm[,c(1,2,5,6,8,12)],
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "SummaryEbICoWGetPsi.txt"),
              row.names=F,quote=F,sep="\t")
  write.table(logfm,
              file.path(cgwasenv$.CGWAS_DETAIL_PATH, "Summaryi-EbICoW.txt"),
              row.names=F,quote=F,sep="\t")
  logOutput("Summary of GetPsi for EbICoW GWASs written to Details/SummaryEbICoWGetPsi.txt\n", cgwasenv = cgwasenv)
  logOutput("Summary of i-EbICoW iteration written to Details/Summaryi-EbICoW.txt\n", cgwasenv = cgwasenv)

  stopCluster(cl)

  logOutput("\nC-GWAS step 4 completed\n", cgwasenv = cgwasenv)
}

# SWaT
step5 <- function(cgwasenv) {
  logOutput("========== C-GWAS step 5 : Truncated Wald Test ==========\n\n", cgwasenv = cgwasenv)
  logOutput("Simulation preparing ..\n", cgwasenv = cgwasenv)

  cl <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(cl)
  # globalVariables(c('i', 'j'))
  i <- j <- 1 # assign parallel control variants

  chunkSize <- 5e4

  Sind <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, "SnpIndex"),
                                          header=T,
                                          nThread = cgwasenv$.PARAL_NUM))
  ACstatm <- read.table(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "EbICoW.txt"),
                        header=T, stringsAsFactors=F)
  maxcn <- (ncol(ACstatm)-2)/4

  corm.tmp <- read.table(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "ACCorrelationStat.txt"),
                         header=T,stringsAsFactors=F)
  statcorm.ac <- gcrom(corm.tmp[,5],length(ACstatm[,1]))
  bgcorm.ac <- gcrom(corm.tmp[,6],length(ACstatm[,1]))
  posm <- ptc(bgcorm.ac, statcorm.ac)
  file.remove(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "ACCorrelationStat.txt"))

  corm.tmp <- read.table(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "BCCorrelationStat.txt"),
                         header=T, stringsAsFactors=F)
  statcorm.bc <- gcrom(corm.tmp[,5],cgwasenv$.TRAIT_NUM)
  bgcorm.bc <- gcrom(corm.tmp[,6],cgwasenv$.TRAIT_NUM)
  posmg <- ptc(bgcorm.bc, statcorm.bc)
  file.remove(file.path(cgwasenv$.CGWAS_DETAIL_PATH, "BCCorrelationStat.txt"))

  logOutput("Max absoluted value of solved EbICoWPsiMatrix : ",
            signif(posm[[3]],7), '\n', cgwasenv = cgwasenv)
  logOutput("Max absoluted value of solved GWASPsiMatrix : ",
            signif(posmg[[3]],7), '\n\n', cgwasenv = cgwasenv)
  logOutput("C-GWAS Simulation start .. (totally ",
               cgwasenv$.SIMUL_DEP * cgwasenv$.IND_SNP_N,
               " simulated independent SNPs)\n", cgwasenv = cgwasenv)

  cutoff.thv <- qchisq(cgwasenv$.TWT_STRAT_CUT, 1, lower.tail = F)
  simumin <- matrix(nrow = 0, ncol = length(cutoff.thv) + 1)
  for (i in 1:ceiling(cgwasenv$.SIMUL_N / 100)) {
    simumin.tmp <- foreach(j = ((i-1)*100+1) : min(i*100, cgwasenv$.SIMUL_N),
                           .inorder = F,
                           .combine = "rbind") %dopar%
                    swtrtsimu(cgwasenv$.SIMUL_SNP_N, posmg[[1]], posm[[1]],
                              ACstatm, maxcn, cutoff.thv, cgwasenv)
    simumin <- rbind(simumin, simumin.tmp); rm(simumin.tmp); gc()
  }
  stopCluster(cl)

  ppn <- ppoints(cgwasenv$.SIMUL_N * cgwasenv$.SIMUL_SNP_N)
  trtc <- c()
  for (i in 1:ncol(simumin)) {
    trtc <- c(trtc, calnna(simumin[,i], ppn)); gc()
  }

  cl <- makeCluster(min(cgwasenv$.PARAL_NUM, 4)) # 4 threads is enough
  registerDoParallel(cl)
  sww <- foreach(i = iter(simumin, by="row", chunksize = chunkSize),
                 .combine = "c",
                 .inorder = T) %dopar%
          tpcor.m.simumin(i, trtc)
  stopCluster(cl)
  rm(simumin); gc()

  logOutput("C-GWAS GetCoef start ..\n", cgwasenv = cgwasenv)

  md <- nullcorrection(sww, "CGWAS", cgwasenv)

  logOutput("C-GWAS Qcoef obtained\n", cgwasenv = cgwasenv)
  logOutput("LOESS model samples written to Details/NullCGWAScorrection.txt\n", cgwasenv = cgwasenv)
  logOutput("Simulated C-GWAS correction effect plotted to Details/NullCGWASdistribution.jpg\n\n", cgwasenv = cgwasenv)
  logOutput("Equivalent independent number of test at study-wide significance : ",
               signif(tpcor2(sum(sww < cgwasenv$.P_THRD_STUDY)/length(sww), md),7), '\n', cgwasenv = cgwasenv)
  logOutput("Equivalent independent number of test at nominal significance : ",
               signif(tpcor2(0.05, md),7), '\n', cgwasenv = cgwasenv)

  cl <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(cl)

  logOutput("Applying the correction to C-GWAS raw p-value .. \n\n", cgwasenv = cgwasenv)
  compn <- as.character(ACstatm[,1])
  tm <- foreach(i = compn, .combine="cbind", .inorder=T) %dopar%
          as.data.frame(data.table::fread(
            file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(i, ".stat")),
                      header = F, stringsAsFactors = F, nThread = 1))[,1]
  owp <- foreach(i = iter(as.matrix(tm), by = "row", chunksize = chunkSize),
                 .combine="rbind",
                 .inorder=T) %dopar%
           swtrt(i, posm[[1]], cutoff.thv)

  owp <- foreach(i = iter(owp, by="row", chunksize = chunkSize),
                 .inorder = T,
                 .combine = "rbind") %dopar%
           tpcor.m(i, trtc)

  nntp <- ntp <- appmin(owp)
  tq <- seq(1/cgwasenv$.IND_SNP_N,
            1-1/cgwasenv$.IND_SNP_N,
            length.out = cgwasenv$.SNP_N)
  ct <- tpcor2(tq, md)
  ntp[order(ntp)] <- ntp[order(ntp)] * ct
  ntp[ntp>1] <- 1
  rm(owp); gc()

  logOutput("MinGWAS Simulation start .. (totally ",
            cgwasenv$.SIMUL_DEP * cgwasenv$.IND_SNP_N,
            " simulated independent SNPs)\n", cgwasenv = cgwasenv)
  gtm <- foreach(i=cgwasenv$.TRAIT_NAME, .combine="cbind", .inorder=T) %dopar%
           as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, paste0(i, ".stat")),
                                           header = F, stringsAsFactors = F, nThread = 1))[,1]

  simup <- foreach(i = 1: cgwasenv$.SIMUL_N,
                   .inorder=F,
                   .combine="c") %dopar%
             efvn(posmg[[1]], cgwasenv$.SIMUL_SNP_N)

  logOutput("MinGWAS GetCoef start ..\n", cgwasenv = cgwasenv)
  md <- nullcorrection(simup, "Minp", cgwasenv)
  logOutput("MinGWAS Qcoef obtained\n", cgwasenv = cgwasenv)
  logOutput("LOESS model samples written to Details/NullMinpcorrection.txt\n", cgwasenv = cgwasenv)
  logOutput("Simulated MinGWAS correction effect plotted to Details/NullMinpdistribution.jpg\n\n", cgwasenv = cgwasenv)
  logOutput("Equivalent independent number of test at study-wide significance : ",
            signif(tpcor2(sum(simup<cgwasenv$.P_THRD_STUDY)/length(simup), md),7), '\n',
            cgwasenv = cgwasenv)
  logOutput("Equivalent independent number of test at nominal significance : ",
            signif(tpcor2(0.05, md),7), '\n', cgwasenv = cgwasenv)

  logOutput("Applying the correction to MinGWAS raw p-value .. \n\n", cgwasenv = cgwasenv)
  ngtp <- gtp <- foreach(i = iter(gtm, by = "row", chunksize = chunkSize),
                         .inorder = T,
                         .combine="c") %dopar% minv(i)

  ct <- tpcor2(tq, md)
  gtp[order(gtp)] <- gtp[order(gtp)] * ct
  gtp[gtp>1] <- 1
  ntp[ntp <= 0] <- max(ntp, rm.na = T)
  gtp[gtp <= 0] <- max(gtp, rm.na = T)
  pm <- cbind(ntp,nntp,gtp,ngtp)
  colnames(pm) <- c("C-GWASadjustedP","C-GWASrawP","MinGWASadjustedP","MinGWASrawP")
  if (cgwasenv$.MAF_FILE_EXIST) {
    MAF <- signif(data.table::fread(file.path(cgwasenv$.CGWAS_iEbICoW_PATH, "MAF"),
                                    header=F, stringsAsFactors=F)[,1], 7)
    MAF <- as.data.frame(MAF)
    data.table::fwrite(cbind(Sind, as.data.frame(signif(pm,7)), MAF),
                       file.path(cgwasenv$.CGWAS_RESULT_PATH, "C-GWAS.p"),
                       sep=" ", na="NA", row.names=F, quote=F)
  } else {
    data.table::fwrite(cbind(Sind, as.data.frame(signif(pm,7))),
                       file.path(cgwasenv$.CGWAS_RESULT_PATH, "C-GWAS.p"),
                       sep=" ", na="NA", row.names=F, quote=F)
  }
  logOutput("Raw P and adjusted P of C-GWAS and MinGWAS written to Result/C-GWAS.p\n", cgwasenv = cgwasenv)

  if (cgwasenv$.SNP_N > 2e6) {
    ssid <- c(1:2e4, seq(2e4, 2e5, 10), seq(2e5, 2e6, 100), seq(2e6, cgwasenv$.SNP_N, 1000), (cgwasenv$.SNP_N-1000):cgwasenv$.SNP_N)
  } else if(cgwasenv$.SNP_N>2e5) {
    ssid <- c(1:2e4, seq(2e4, 2e5, 10), seq(2e5, cgwasenv$.SNP_N, 100), (cgwasenv$.SNP_N-1000):cgwasenv$.SNP_N)
  } else {
    ssid <- c(1:2e4,seq(2e4,cgwasenv$.SNP_N,10),(cgwasenv$.SNP_N-1000):cgwasenv$.SNP_N)
  }

  xx <- -log10(ppoints(cgwasenv$.SNP_N)[ssid])
  yyc <- -log10(ntp[order(ntp)][ssid])
  yyg <- -log10(gtp[order(gtp)][ssid])
  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWASminpQQ.jpg"),
       width=3000, height=3000, res=600)
  par(mar=c(3.5, 3.5, 1, 1))
  plot(NA, xlim=c(0, max(xx)), ylim=c(0, max(c(yyc[1], yyg[1]))),
       xlab=expression(paste("Expected   ", -Log[10](italic(p)))),
       ylab=expression(paste("Observed   ", -Log[10](italic(p)))),
       mgp=c(2, 0.7, 0), las=1, cex.axis=0.85, tck=-0.015)
  points(xx, yyg, pch=20, col="#1161C9", cex=0.6)
  points(xx, yyc, pch=20, col="#FD9001", cex=0.6)
  abline(c(0, 1), lwd=1)
  dev.off()
  logOutput("QQ plots of C-GWAS and MinGWAS adjusted P plotted to Result/CGWASminpQQ.jpg\n", cgwasenv = cgwasenv)

  if (!cgwasenv$.KEEP_EbICoW) {
    unlink(cgwasenv$.CGWAS_iEbICoW_PATH, recursive=TRUE)
    logOutput("All i-EbICoW statistics dropped\n", cgwasenv = cgwasenv)
  } else{
    logOutput("All i-EbICoW statistics kept in Details/i-EbICoW/\n", cgwasenv = cgwasenv)
  }

  logOutput("\nC-GWAS step 5 completed\n", cgwasenv = cgwasenv)
}

# Summary
step6 <- function(cgwasenv) {
  logOutput("========== C-GWAS step 6 : Summary of C-GWAS ==========\n\n", cgwasenv = cgwasenv)

  cl <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(cl)
  # globalVariables(c('i'))
  i <- 1 # assign parallel control variants

  cgwasm <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_RESULT_PATH, "C-GWAS.p"),
                                            header=T))
  Sind <- cgwasm[,1:3]
  newtick <- manpos(Sind)
  Sind <- cbind(Sind, newtick[[1]])
  newtick <- newtick[[2]]
  if(cgwasenv$.MAF_FILE_EXIST) {
    mafv <- cgwasm[,8]
  }
  pm <- pm <- as.matrix(cgwasm[,4:6])

  fdrv <- rep(cgwasenv$.P_THRD_SUGST, 3)
  gv <- rep(cgwasenv$.P_THRD_STUDY,3)

  ns <- which(pm[,1]<=fdrv[1])
  os <- which(pm[,3]<=fdrv[3])

  lcidn <- locisearch(pm[ns, 1], Sind[ns, 1:2], fdrv[1], cgwasenv)
  lcido <- locisearch(pm[os, 3], Sind[os, 1:2], fdrv[3], cgwasenv)

  n1 <- length(os)-length(intersect(os,ns))
  n2 <- length(intersect(os,ns))
  n3 <- length(ns)-length(intersect(os,ns))

  tableid <- tableid <- union(ns,os)
  tableid <- tableid[order(tableid)]
  temnp <- pm[ns,1]
  temop <- pm[os,3]
  rkn <- tlcidn <- rep(NA, length(ns))
  rko <- tlcido <- rep(NA, length(os))
  if(length(lcidn)!=0) {
    for(j in 1:length(lcidn)) {
      tlcidn[lcidn[[j]]] <- j
      rkn[lcidn[[j]]] <- rank(temnp[lcidn[[j]]], ties.method="random")
    }
  }
  if(length(lcido)!=0) {
    for(j in 1:length(lcido)) {
      tlcido[lcido[[j]]] <- j
      rko[lcido[[j]]] <- rank(temop[lcido[[j]]], ties.method="random")
    }
  }
  if(length(tableid)!=1) {
    nm <- rbind(cbind(tlcidn, rkn, rep(1, length(rkn))),
                matrix(NA, length(tableid)-length(ns), 3))[order(c(ns, setdiff(tableid, ns))),]
    om <- rbind(cbind(tlcido, rko, rep(1, length(rko))),
                matrix(NA, length(tableid)-length(os), 3))[order(c(os, setdiff(tableid, os))),]
  } else{
    nm <- matrix(rbind(cbind(tlcidn, rkn, rep(1, length(rkn))),
                       matrix(NA, length(tableid)-length(ns), 3))
                 [order(c(ns, setdiff(tableid, ns))),], nrow=length(tableid))
    om <- matrix(rbind(cbind(tlcido, rko, rep(1, length(rko))),
                       matrix(NA, length(tableid)-length(os), 3))
                 [order(c(os, setdiff(tableid, os))),], nrow=length(tableid))
  }
  nm <- cbind(signif(pm[tableid, 1], 7), nm)
  om <- cbind(signif(pm[tableid, 3], 7), om)
  nm[nm[,1]<=5e-8,4] <- 2
  om[om[,1]<=5e-8,4] <- 2

  if (cgwasenv$.MAF_FILE_EXIST) {
    sumtable <- cbind(Sind[tableid,1:3], mafv[tableid], nm[,1:2], om[,1:2])
    colnames(sumtable) <- c("CHR","BP","SNP","MAF","C-GWAS-P","C-GWAS-Loci","MinGWAS-P","MinGWAS-Loci")
  } else {
    sumtable <- cbind(Sind[tableid,1:3], nm, om)
    colnames(sumtable) <- c("CHR","BP","SNP","C-GWAS-P","C-GWAS-Loci","MinGWAS-P","MinGWAS-Loci")
  }
  write.csv(sumtable,
            file.path(cgwasenv$.CGWAS_RESULT_PATH, "SummarySugSigSNP.csv"),
            row.names=F,quote=F)
  logOutput("Suggestive significant SNP summary table written to Result/SummarySugSigSNP.csv\n", cgwasenv = cgwasenv)

  ns <- which(pm[,1]<=5e-8)
  os <- which(pm[,3]<=5e-8)
  lcidn <- locisearch(pm[ns, 1], Sind[ns, 1:2], 5e-8, cgwasenv)
  lcido <- locisearch(pm[os, 3], Sind[os, 1:2], 5e-8, cgwasenv)
  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWAS-GWAS-S.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm3 <- manhattan(pm[,3], pm[,1], os, ns, lcido, lcidn, fdrv[3], fdrv[1], gv[3], gv[1],
                     Sind, newtick)
  dev.off()
  logOutput("Manhattan plots of C-GWAS and MinGWAS plotted to Result/CGWAS-GWAS.jpg\n\n", cgwasenv = cgwasenv)

  n4 <- length(os)-length(intersect(os,ns))
  n5 <- length(intersect(os,ns))
  n6 <- length(ns)-length(intersect(os,ns))

  logOutput("At suggestive significance ", signif(cgwasenv$.P_THRD_SUGST),
            " : (SNP)\n", cgwasenv = cgwasenv)
  logOutput("MinGWAS unique : ", n1,
            ",  MinGWAS and C-GWAS shared : ", n2,
            ",  C-GWAS unique : ", n3, '\n\n', cgwasenv = cgwasenv)
  logOutput("At study-wide significance ", signif(cgwasenv$.P_THRD_STUDY,7),
            " : (Loci/SNP)\n", cgwasenv = cgwasenv)
  logOutput("MinGWAS unique : ", lcsm3[1,1], "/", n4,
            ",  MinGWAS and C-GWAS shared : ", min(lcsm3[2,1],lcsm3[1,2]), "/", n5,
            ",  C-GWAS unique : ", lcsm3[2,2], "/", n6, '\n', cgwasenv = cgwasenv)

  stopCluster(cl)

  logOutput("\nC-GWAS step 6 completed\n", cgwasenv = cgwasenv)
}
