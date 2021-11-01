# Statistic Summary Extraction
step1 <- function(cgwasenv) {
  if (cgwasenv$.TRAIT_NUM <= 1)  {
    stop("Error: Traits number must larger than 1.")
  } # check equality among all element in snp.N

  Sepblock <- min(cgwasenv$.PARAL_NUM, cgwasenv$.TRAIT_NUM)
  cl <- makeCluster(Sepblock)
  registerDoParallel(cl)
  # globalVariables(c("i"))
  i <- 1 # assign parallel control variants

  naid <- unique(foreach(i = seq_len(cgwasenv$.TRAIT_NUM), .combine="c", .inorder=F) %dopar%
                   StatE1(i, cgwasenv))
  if (cgwasenv$.MAF_FILE_EXIST) {
    mafv <- as.data.frame(data.table::fread(cgwasenv$.MAF_FILE_PATH,
                                            header=F, stringsAsFactors=F,
                                            nThread = cgwasenv$.PARAL_NUM))[,1]
    if (length(naid) != 0) {
      mafv <- mafv[-naid]
    }
    write.table(mafv, file.path(cgwasenv$.CGWAS_COLDATA_PATH, "MAF"), row.names=F, col.names=F, quote=F)
  }
  snp.N <- foreach(i = seq_len(cgwasenv$.TRAIT_NUM), .combine="c", .inorder = F) %dopar%
    StatE2(i, naid, cgwasenv)

  if (cgwasenv$.MAF_FILE_EXIST) {snp.N <- c(snp.N, length(mafv))}
  if (length(unique(snp.N)) > 1)  {
    stop("Error: SNP number in each files are not equal.")
  } # check equality among all element in snp.N

  write.table(naid, file.path(cgwasenv$.CGWAS_COLDATA_PATH, 'NA_list'), row.names = F, col.names = F)

  stopCluster(cl)
}

# Inflation Correlation Estimation
step2 <- function(cgwasenv) {
  inseq <- ppoints(2000)
  ranidlist <- ridl.ICE(repn = 1, minsnpn = 1e5, cgwasenv)

  Sepblock <- min(cgwasenv$.PARAL_NUM, nrow(ranidlist))
  cl <- makeCluster(Sepblock)
  registerDoParallel(cl)
  # globalVariables(c('i', 'j'))
  i <- j <- 1 # assign parallel control variants

  resinfm <- matrix(NA, cgwasenv$.TRAIT_NUM , 5)

  # system.time({
  for(i in seq_len(cgwasenv$.TRAIT_NUM)) {
    tm <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[i], ".stat")),
                                      header = F,
                                      nThread = cgwasenv$.PARAL_NUM))[,1]
    tm2 <- tm^2
    resinfm[i,1] <- mean(tm2)
    resinfm[i,2] <- median(tm2)/qchisq(0.5, 1)
    mp <- 0.1
    rantm <- matrix(tm2[ranidlist], nrow(ranidlist))
    infm <- foreach(j=iter(rantm, by="row", chunksize=1),.combine="rbind",.inorder=T) %dopar%
              calinf(j, inseq, mp)
    infm[infm[,2] > mp, 2] <- mp
    while((sum(infm[,2] == mp) > (nrow(ranidlist)*0.1)) & (mp!=0.5)) {
      mp <- mp+0.1
      infm <- foreach(j=iter(rantm, by="row", chunksize=1),.combine="rbind",.inorder=T) %dopar%
                calinf(j, inseq, mp)
      infm[infm[,2] > mp, 2] <- mp
    }
    if(mp!=0.5) {
      resinfm[i, 3] <- mean(infm[infm[, 2]!=mp, 1])
      resinfm[i, 5] <- mean(infm[infm[, 2]!=mp, 2])
    } else{
      resinfm[i, 3] <- mean(infm[, 1])
      resinfm[i, 5] <- mean(infm[, 2])
    }
    if((resinfm[i, 1]-resinfm[i, 3]) < 1e-3) {
      resinfm[i, 4] <- 1.001
    } else{
      resinfm[i, 4] <- 1+resinfm[i, 1]-resinfm[i, 3]
    }
    data.table::fwrite(as.data.frame(signif(tm/sqrt(resinfm[i, 1]/resinfm[i, 4]), 7)),
                       file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(cgwasenv$.TRAIT_NAME[i], ".efstat")),
                       row.names = F, col.names = F, quote = F)
  }
  # })

  colnames(resinfm) <- c("MeanX2", "Lambda", "EstInf", "ReaEff", "Eprop")
  write.table(cbind(TraitName = cgwasenv$.TRAIT_NAME, signif(resinfm, 7)),
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "InflationStat.txt"),
              col.names = T, row.names = F, quote = F, sep = "\t")

  stopCluster(cl)

  pairma <- t(combn(seq_len(cgwasenv$.TRAIT_NUM), 2))
  n <- nrow(pairma)

  upnum <- ceiling(n / cgwasenv$.PARAL_NUM)
  upblock <- ceiling(n / upnum)

  Sepblock <- min(cgwasenv$.PARAL_NUM, upblock)

  cl <- makeCluster(Sepblock)
  registerDoParallel(cl)

  # system.time({
  tresm <- foreach(i = 1:upblock,
                   .combine = "rbind",
                   .inorder=T) %dopar%
             CorE.ICE(i, upnum, n, pairma, resinfm, ranidlist, cgwasenv)
  # })
  tresm <- signif(tresm, 7)

  corm <- cbind(cgwasenv$.TRAIT_NAME[pairma[,1]],
                cgwasenv$.TRAIT_NAME[pairma[,2]],
                signif(resinfm[pairma[,1], 4]-1, 7),
                signif(resinfm[pairma[,2], 4]-1, 7),
                tresm)

  colnames(corm) <- c("Trait1", "Trait2", "T1Eff", "T2Eff", "StatCor", "BgCor", "EffCov", "EffCor", "T1sEff", "T2sEff", "EffsCov", "EffsCor")
  write.table(corm,
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "BCCorrelationStat.txt"),
              row.names = F, quote = F, sep = "\t")

  stopCluster(cl)
}

# EbICo
step3 <- function(cgwasenv) {
  Sepblock <- min(cgwasenv$.PARAL_NUM, cgwasenv$.TRAIT_NUM)
  cl <- makeCluster(Sepblock)
  registerDoParallel(cl)
  # globalVariables(c('i'))
  i <- 1 # assign parallel control variants

  ranidlist <- ridl.ebico(repn = 1, minsnpn = 1e5, cgwasenv)

  thresc <- qchisq(cgwasenv$.P_THRD, 1, lower.tail=F)
  thresw <- qchisq(cgwasenv$.P_THRD, 2, lower.tail=F)

  if(cgwasenv$.MAF_FILE_EXIST) {
    mafv <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "MAF"),
                                            header=F, stringsAsFactors=F,
                                            nThread = cgwasenv$.PARAL_NUM))[,1]
  }

  corm.tmp <- read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "BCCorrelationStat.txt"),
                       header=T, stringsAsFactors=F)
  coridm <- as.matrix(corm.tmp[,1:2])
  cordatm <- as.matrix(corm.tmp[,3:12])

  localsep <- 1e5

  n <- cgwasenv$.TRAIT_NUM+1
  AllEff <- as.numeric(c(cordatm[-(1:(n-2))^2/2+(n-0.5)*(1:(n-2))-n+2, 1],
                         cordatm[-(n-2)^2/2+(n-0.5)*(n-2)-n+2, 2]))
  Ess <- unlist(foreach(i = cgwasenv$.TRAIT_NAME, .inorder=T) %dopar%
                  Essfun(i, mafv, cgwasenv))
  TNm <- cbind(data.frame(TraitName = cgwasenv$.TRAIT_NAME, stringsAsFactors=F), AllEff, Ess)
  idlist <- as.list(cgwasenv$.TRAIT_NAME)
  bflist <- as.list(rep(1, cgwasenv$.TRAIT_NUM))
  tflist <- as.list(rep(1, cgwasenv$.TRAIT_NUM))
  cdlist <- as.list(rep(1, cgwasenv$.TRAIT_NUM))
  rnlist <- as.list(rep(NA, cgwasenv$.TRAIT_NUM))

  bkcordatm <- c()
  bknum <- 1
  logfm <- c()
  rn <- 1

  while(!(all(abs(cordatm[,6] - cordatm[,4]) + abs(cordatm[,10] - cordatm[,4]) < cgwasenv$.MIN_CORR_DIFF) &
          all(abs(cordatm[,4]) < cgwasenv$.HIGH_CORR_RES))) {
    mtarid <- order(abs(cordatm[,6]-cordatm[,4]+cordatm[,10]-cordatm[,4]), decreasing=T)[1]
    selidv <- coridm[mtarid,]
    newname <- paste0(selidv[1], "_", selidv[2])
    t2m <- cbind(
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(selidv[1], ".efstat")),
                                         header=F, stringsAsFactors=F,
                                         nThread = cgwasenv$.PARAL_NUM))[,1],
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(selidv[2], ".efstat")),
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

    x2m <- apply(temsm, 2, mean)

    efin <- 2 -
            0.0725615289639716 * abs(cordatm[mtarid, 4]) -
            0.196105591438997 * cordatm[mtarid, 4] ^ 2 -
            0.192036017551772 * cordatm[mtarid, 4] ^ 4 -
            0.134862184571824 * cordatm[mtarid, 4] ^ 8 -
            0.0969596906274313 * cordatm[mtarid, 4] ^ 16 -
            0.260927307332186 * cordatm[mtarid, 4] ^ 64 +
            0.232255565305109 * cordatm[mtarid, 4] ^ 128 -
            0.27795687154795 * cordatm[mtarid, 4] ^ 256
    ge1 <- which((ttm[,1]^2) > qchisq(cgwasenv$.P_THRD / efin, 1, lower.tail=F))
    ge2 <- which((ttm[,2]^2) > qchisq(cgwasenv$.P_THRD / efin, 1, lower.tail=F))
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
    gen <- max(c(length(idg), length(ge1c), length(ge2c)))
    geon <- length(idog)
    een <- length(ide)
    esn <- length(ids)
    ean <- length(ida)
    ewn <- length(idw)
    jee <- length(intersect(idgc, ide))
    jes <- length(intersect(idgc, ids))
    jea <- length(intersect(idgc, ida))
    jew <- length(intersect(idgc, idw))
    jee1 <- length(intersect(ge1c, ide))
    jeeo <- length(intersect(idog, ide))
    jee2 <- length(intersect(ge2c, ide))
    jes1 <- length(intersect(ge1c, ids))
    jeso <- length(intersect(idog, ids))
    jes2 <- length(intersect(ge2c, ids))
    jea1 <- length(intersect(ge1c, ida))
    jeao <- length(intersect(idog, ida))
    jea2 <- length(intersect(ge2c, ida))
    jew1 <- length(intersect(ge1c, idw))
    jewo <- length(intersect(idog, idw))
    jew2 <- length(intersect(ge2c, idw))

    print(paste0(rn, " Round    ", sum((abs(cordatm[,6]-cordatm[,4])+abs(cordatm[,10]-cordatm[,4])>cgwasenv$.MIN_CORR_DIFF)|(abs(cordatm[,4])>cgwasenv$.HIGH_CORR_RES))-1, " Left"))
    print(paste0(selidv[1], "  ", signif(cordatm[mtarid, 1], 3), "  ", signif(cordatm[mtarid, 7], 3), "  ", gen1))
    print(paste0(selidv[2], "  ", signif(cordatm[mtarid, 2], 3), "  ", signif(cordatm[mtarid, 8], 3), "  ", gen2))
    print(paste0(newname, "  (", signif(cordatm[mtarid, 4], 3), ")  ", geon, "/", gen, "(", signif(geon/gen*100, 3), "%)"))
    print(paste0("using AEC  (", signif(cordatm[mtarid, 6], 3), ")  ", signif(cofm[1, 1], 3), "  ", signif(cofm[2, 1], 3), " -> ", signif(x2m[1], 3), "  ", een, "/", gen, "(", signif(een/gen*100, 3), "%)"))
    print(paste0("using SEC  (", signif(cordatm[mtarid, 10], 3), ")  ", signif(cofm[1, 2], 3), "  ", signif(cofm[2, 2], 3), " -> ", signif(x2m[2], 3), "  ", esn, "/", gen, "(", signif(esn/gen*100, 3), "%)"))
    print(paste0("using SSC  (", signif(fm[2, 9], 3), ")  ", signif(cofm[1, 3], 3), "  ", signif(cofm[2, 3], 3), " -> ", signif(x2m[3], 3), "  ", ean, "/", gen, "(", signif(ean/gen*100, 3), "%)"))
    print(paste0("using WDT -> ", ewn, "/", gen, "(", signif(ewn/gen*100, 3), "%)"))
    print(paste0("using AEC  ", jee, "/", gen, "(", signif(jee/gen*100, 3), "%)  ", jee1, "/", gen1, "(", signif(jee1/gen1*100, 3), "%)  ", jeeo, "/", geon, "(", signif(jeeo/geon*100, 3), "%)  ", jee2, "/", gen2, "(", signif(jee2/gen2*100, 3), "%)"))
    print(paste0("using SEC  ", jes, "/", gen, "(", signif(jes/gen*100, 3), "%)  ", jes1, "/", gen1, "(", signif(jes1/gen1*100, 3), "%)  ", jeso, "/", geon, "(", signif(jeso/geon*100, 3), "%)  ", jes2, "/", gen2, "(", signif(jes2/gen2*100, 3), "%)"))
    print(paste0("using SSC  ", jea, "/", gen, "(", signif(jea/gen*100, 3), "%)  ", jea1, "/", gen1, "(", signif(jea1/gen1*100, 3), "%)  ", jeao, "/", geon, "(", signif(jeao/geon*100, 3), "%)  ", jea2, "/", gen2, "(", signif(jea2/gen2*100, 3), "%)"))
    print(paste0("using WDT  ", jew, "/", gen, "(", signif(jew/gen*100, 3), "%)  ", jew1, "/", gen1, "(", signif(jew1/gen1*100, 3), "%)  ", jewo, "/", geon, "(", signif(jewo/geon*100, 3), "%)  ", jew2, "/", gen2, "(", signif(jew2/gen2*100, 3), "%)"))
    rn <- rn+1

    if((een>=esn)&(jeeo>=jeso)&(jee>=jes)) {
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

    if((c(een, esn, ean)[did]<max(c(gen*cgwasenv$.MIN_EbiCo_POWER_INC, ewn*cgwasenv$.MIN_EbiCo_POWER_INC)))&(abs(cordatm[mtarid, 4])<cgwasenv$.HIGH_CORR_RES)) {
      print(paste0("Decision -> do not combine"))
      print("")
      bkcordatm <- rbind(bkcordatm, cordatm[mtarid,])
      logfm <- rbind(logfm, c(rn-1, selidv, newname, NA, NA, signif(c(cordatm[mtarid, c(4, 1, 2, 6, 7, 8, 10)], gen1, gen2, gen, geon, x2m, cofm[,1], cofm[,2], cofm[,3], een, esn, ean, ewn, jee, jes, jea, jew, jee1, jeeo, jee2, jes1, jeso, jes2, jea1, jeao, jea2, jew1, jewo, jew2), 7), rep(NA, 5)))
      cordatm[mtarid,4] <- cordatm[mtarid,6] <- cordatm[mtarid,10] <- 0
      cordatm[mtarid,3] <- bknum+1
      bknum <- bknum+1
      next
    }

    newt <- t2m %*% cofm[,did]
    x2m[did] <- mean(newt^2)-1
    if(x2m[did]<1e-4) {
      newt <- newt*sqrt(1.001/(x2m[did]+1))
      x2m[did] <- mean(newt^2)-1
    }

    bofm <- fm[,4:6]
    dofm <- fm[,7:9]

    b2m <- cbind(
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(selidv[1], ".beta")),
                                         header=F, stringsAsFactors=F,
                                         nThread = cgwasenv$.PARAL_NUM))[,1],
                 as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(selidv[2], ".beta")),
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

    print(paste0("Equivalent sample size  ", round(mse[1]), " ", round(mse[2]), " -> ", round(newmse)))
    print(paste0("Mean beta weights  ", signif(bofv[1], 3), " ", signif(bofv[2], 3), "   total ", sum((bofv*dofm[,did]))))
    if((newmse<max(c(mse[1]+mse[2]*cgwasenv$.SAMPLE_SIZE_INC, mse[2]+mse[1]*cgwasenv$.SAMPLE_SIZE_INC))) &
       (abs(cordatm[mtarid, 4])<cgwasenv$.HIGH_CORR_RES)) {
      print(paste0("Decision -> do not combine"))
      print("")
      bkcordatm <- rbind(bkcordatm, cordatm[mtarid,])
      logfm <- rbind(logfm, c(rn-1, selidv, newname, c("AEC", "SEC", "SSC")[did], NA, signif(c(cordatm[mtarid, c(4, 1, 2, 6, 7, 8, 10)], gen1, gen2, gen, geon, x2m, cofm[,1], cofm[,2], cofm[,3], een, esn, ean, ewn, jee, jes, jea, jew, jee1, jeeo, jee2, jes1, jeso, jes2, jea1, jeao, jea2, jew1, jewo, jew2), 7), round(c(mse, newmse)), signif(bofv, 7)))
      cordatm[mtarid,4] <- cordatm[mtarid,6] <- cordatm[mtarid,10] <- 0
      cordatm[mtarid,3] <- bknum+1
      bknum <- bknum+1
      next
    } else{
      print(paste0("Decision -> use ", c("AEC", "SEC", "SSC")[did]))
      print("")
    }

    logfm <- rbind(logfm, c(rn-1, selidv, newname, c("AEC", "SEC", "SSC")[did], TRUE, signif(c(cordatm[mtarid, c(4, 1, 2, 6, 7, 8, 10)], gen1, gen2, gen, geon, x2m, cofm[,1], cofm[,2], cofm[,3], een, esn, ean, ewn, jee, jes, jea, jew, jee1, jeeo, jee2, jes1, jeso, jes2, jea1, jeao, jea2, jew1, jewo, jew2), 7), round(c(mse, newmse)), signif(bofv, 7)))

    if (nrow(coridm) < 2) break

    selid1 <- which(TNm[,1]==selidv[1])
    selid2 <- which(TNm[,1]==selidv[2])
    TNm[selid1,] <- NA
    TNm[selid2,] <- NA
    newcordatm <- foreach(i=which(!is.na(TNm[,1])), .combine="rbind", .inorder=T) %dopar%
                    CorE.ebico(i, newt, TNm, x2m, did, ranidlist, cgwasenv)
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
    TNm[nrow(TNm),2] <- x2m[did]
    TNm[nrow(TNm),3] <- newmse
    idlist[[n]] <- c(idlist[[selid1]], idlist[[selid2]])
    idlist[[selid1]] <- idlist[[selid2]] <- NA
    bflist[[n]] <- c(bflist[[selid1]]*bofv[1], bflist[[selid2]]*bofv[2])
    bflist[[selid1]] <- bflist[[selid2]] <- NA
    tflist[[n]] <- c(tflist[[selid1]]*cofm[1, did], tflist[[selid2]]*cofm[2, did])
    tflist[[selid1]] <- tflist[[selid2]] <- NA
    cdlist[[n]] <- c(cdlist[[selid1]]*dofm[1, did], cdlist[[selid2]]*dofm[2, did])
    cdlist[[selid1]] <- cdlist[[selid2]] <- NA
    rnlist[[n]] <- c(rnlist[[selid1]], rnlist[[selid2]], rn-1)
    rnlist[[n]] <- rnlist[[n]][!is.na(rnlist[[n]])]
    rnlist[[selid1]] <- rnlist[[selid2]] <- NA
    n <- n+1

    data.table::fwrite(as.data.frame(signif(newt, 7)),
                       file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(newname, ".efstat")),
                       row.names=F,col.names=F,quote=F)
    data.table::fwrite(as.data.frame(signif(newb, 7)),
                       file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(newname, ".beta")),
                       row.names=F,col.names=F,quote=F)
  }

  compm <- bfm <- tfm <- cdm <- rnm <- as.data.frame(rep(NA, sum(!is.na(TNm[,1]))))
  tn <- 1
  for(i in 1:(n-1)) {
    if(!is.na(idlist[[i]][1])) {
      if(length(idlist[[i]])>ncol(compm)) {
        compm <- cbind(compm, t(rep(NA, length(idlist[[i]])-ncol(compm))))
        bfm <- cbind(bfm, t(rep(NA, length(bflist[[i]])-ncol(bfm))))
        tfm <- cbind(tfm, t(rep(NA, length(tflist[[i]])-ncol(tfm))))
        cdm <- cbind(cdm, t(rep(NA, length(cdlist[[i]])-ncol(cdm))))
        rnm <- cbind(rnm, t(rep(NA, length(rnlist[[i]])-ncol(rnm))))
      }
      compm[tn,1:length(idlist[[i]])] <- idlist[[i]]
      bfm[tn,1:length(bflist[[i]])] <- bflist[[i]]
      tfm[tn,1:length(tflist[[i]])] <- tflist[[i]]
      cdm[tn,1:length(cdlist[[i]])] <- cdlist[[i]]
      rnm[tn,1:length(rnlist[[i]])] <- rnlist[[i]]
      tn <- tn+1
    }
  }

  colnames(compm) <- paste0("N", 1:ncol(compm))
  colnames(bfm) <- paste0("bf", 1:ncol(bfm))
  colnames(tfm) <- paste0("tf", 1:ncol(tfm))
  colnames(cdm) <- paste0("d", 1:ncol(cdm))
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

  colnames(logfm) <- c("Round", "Trait1", "Trait2", "Newtrait", "CovType", "IFcom",
                       colnames(cordatm)[c(4, 1, 2, 6, 7, 8, 10)],
                       "SSNPN1", "SSNPN2", "SSNPNall", "SSNPNcom",
                       "meanX21", "meanX22", "meanX23",
                       "AECtf1", "AECtf2", "SECtf1", "SECtf2", "SSCtf1", "SSCtf2",
                       "SSNPAEC", "SSNPSEC", "SSNPSSC", "SSNPWD",
                       "SSNPAECall", "SSNPSECall", "SSNPSSCall", "SSNPWDall",
                       "SSNPAECN1", "SSNPAECcom", "SSNPAECN2", "SSNPSECN1",
                       "SSNPSECcom", "SSNPSECN2", "SSNPSSCN1", "SSNPSSCcom",
                       "SSNPSSCN2", "SSNPWDN1", "SSNPWDcom", "SSNPWDN2",
                       "Ess1", "Ess2", "EssNew", "bf1", "bf2")

  write.table(cbind(TNm, compm, signif(bfm, 7), signif(tfm, 7), cdm, rnm),
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "ACEffect.txt"),
              row.names=F,quote=F,sep="\t")

  accorm <- cbind(matrix(coridm[orderid,], ncol=ncol(coridm)), matrix(signif(cordatm[orderid,], 7), ncol=ncol(cordatm)))
  colnames(accorm) <- c(colnames(coridm), colnames(cordatm))
  write.table(accorm,
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "ACCorrelationStat.txt"),
              row.names=F,quote=F,sep="\t")
  write.table(logfm,
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "EBICOlog.txt"),
              row.names=F,quote=F,sep="\t")

  stopCluster(cl)
}

# SWaT
step4 <- function(cgwasenv) {
  cl <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(cl)
  # globalVariables(c('i', 'j'))
  i <- j <- 1 # assign parallel control variants

  localsep <- 5e4

  Sind <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SnpIndex"),
                                          header=T,
                                          nThread = cgwasenv$.PARAL_NUM))
  ACstatm <- read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "ACEffect.txt"),
                        header=T,stringsAsFactors=F)
  maxcn <- (ncol(ACstatm)-2)/5

  statcorm <- gcrom(read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "ACCorrelationStat.txt"),
                               header=T,stringsAsFactors=F)[,5],length(ACstatm[,1]))
  bgcorm <- gcrom(read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "ACCorrelationStat.txt"),
                             header=T,stringsAsFactors=F)[,6],length(ACstatm[,1]))
  posm <- ptc(bgcorm, statcorm)

  corm.tmp <- read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "BCCorrelationStat.txt"),
                         header=T, stringsAsFactors=F)
  statcorm <- gcrom(corm.tmp[,5],cgwasenv$.TRAIT_NUM)
  bgcorm <- gcrom(corm.tmp[,6],cgwasenv$.TRAIT_NUM)

  posmg <- ptc(bgcorm, statcorm)

  cutoff.thv <- qchisq(cgwasenv$.SWaT_STRAT_CUT, 1, lower.tail = F)
  ts <- Sys.time() # test
  simumin <- foreach(i = 1: cgwasenv$.SIMUL_N,
                     .inorder = F,
                     .combine = "rbind") %dopar%
               swtrtsimu(cgwasenv$.SIMUL_SNP_N, posmg[[1]], posm[[1]], ACstatm, maxcn, cutoff.thv, cgwasenv)
  clusterEvalQ(cl, gc()) # collect garbage for each thread
  print("time consuming of swtrtsimu: "); print(difftime(Sys.time(), ts)) # test
  stopCluster(cl)
  gc()

  ppn <- ppoints(cgwasenv$.SIMUL_N * cgwasenv$.SIMUL_SNP_N)
  ts <- Sys.time() # test
  trtc <- unlist(apply(simumin, 2, function(x) calnna(x, ppn)))
  gc()
  print("time consuming of calnna: "); print(difftime(Sys.time(), ts)) # test

  md1 <- nullcorrection(simumin[,ncol(simumin)], "EbICoW", cgwasenv)
  pop1 <- sum(simumin[,ncol(simumin)]<5e-8)/nrow(simumin)

  cl <- makeCluster(4) # 4 threads is enough
  registerDoParallel(cl)
  ts <- Sys.time() # test
  simumin1 <- foreach(i = iter(simumin, by="row", chunksize = localsep),
                     .inorder = T) %dopar%
                tpcor.m(i, trtc)

  print("time consuming of tpcor: "); print(difftime(Sys.time(), ts)) # test
  stopCluster(cl)
  sww <- unlist(lapply(simumin1, appmin))
  rm(simumin); rm(simumin1); gc()

  ts <- Sys.time() # test
  md <- nullcorrection(sww, "CGWAS", cgwasenv)
  pop <- sum(sww<5e-8) / length(sww)
  pop6 <- sum(sww<5e-8 * trtc[length(trtc)]) / length(sww)
  pop5 <- sum(sww<0.05) / length(sww)
  MT <- matrix(c("Combination.study-wide", tpcor2(pop, md), tpcor2(pop, md)*5e-8,
                 "EbICoW.study-wide", tpcor2(pop1, md1), tpcor2(pop1, md1)*5e-8,
                 "C-GWAS.study-wide", tpcor2(pop6, md)*trtc[length(trtc)], tpcor2(pop6, md)*trtc[length(trtc)]*5e-8,
                 "C-GWAS.nominal", tpcor2(pop5, md)*trtc[length(trtc)], tpcor2(pop5, md)*trtc[length(trtc)]*5e-8),
               ncol=3, byrow=T)
  print("time consuming of calc MT: "); print(difftime(Sys.time(), ts)) # test

  cl <- makeCluster(cgwanenv$.PARAL_NUM)
  compn <- as.character(ACstatm[,1])
  tm <- foreach(i = compn, .combine="cbind", .inorder=T) %dopar%
          as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(i, ".efstat")),
                                                      header = F, stringsAsFactors = F, nThread = 1))[,1]
  ts <- Sys.time() # test
  owp <- foreach(i = iter(as.matrix(tm), by = "row", chunksize = localsep),
                 .combine="rbind",
                 .inorder=T) %dopar%
           swtrt(i, posm[[1]], cutoff.thv)
  clusterEvalQ(cl, gc()) # collect garbage for each thread
  print("time consuming of swtrt: "); print(difftime(Sys.time(), ts)) # test

  otp <- owp[, ncol(owp)]
  tq <- seq(1/cgwasenv$.IND_SNP_N, 1-1/cgwasenv$.IND_SNP_N, length.out=length(otp))
  ct <- tpcor2(tq, md1)
  otp[order(otp)] <- otp[order(otp)]*ct
  otp[otp>1] <- 1

  ts <- Sys.time() # test
  owp <- foreach(i = iter(owp, by="row", chunksize = localsep),
                 .inorder = T,
                 .combine = "rbind") %dopar%
           tpcor.m(i, trtc)
  clusterEvalQ(cl, gc()) # collect garbage for each thread
  print("time consuming of tpcor: "); print(difftime(Sys.time(), ts)) # test

  data.table::fwrite(as.data.frame(signif(owp, 7)),
                     file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SWaT.p"),
                     sep=" ",na="NA",row.names=F,col.names=F,quote=F)

  ntp <- appmin(owp)
  ct <- tpcor2(tq, md)
  ntp[order(ntp)] <- ntp[order(ntp)] * ct
  ntp[ntp>1] <- 1
  rm(owp); gc()

  gtm <- foreach(i=cgwasenv$.TRAIT_NAME, .combine="cbind", .inorder=T) %dopar%
           as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(i, ".efstat")),
                                           header = F, stringsAsFactors = F, nThread = 1))[,1]

  simup <- foreach(i = 1: cgwasenv$.SIMUL_N,
                   .inorder=F,
                   .combine="c") %dopar%
             efvn(posmg[[1]], cgwasenv$.SIMUL_SNP_N)

  gtp <- foreach(i = iter(gtm, by = "row", chunksize = localsep),
                 .inorder = T,
                 .combine="c") %dopar%
           minv(i)

  md <- nullcorrection(simup, "Minp", cgwasenv)
  pop <- sum(simup<5e-8)/length(simup)
  MT <- rbind(MT, matrix(c("Minp.study-wide", tpcor2(pop, md), tpcor2(pop, md)*5e-8),
                         ncol=3, byrow=T))
  ct <- tpcor2(tq, md)
  gtp[order(gtp)] <- gtp[order(gtp)] * ct
  gtp[gtp>1] <- 1

  ntp[ntp <= 0] <- max(ntp, rm.na = T)
  otp[otp <= 0] <- max(otp, rm.na = T)
  gtp[gtp <= 0] <- max(gtp, rm.na = T)
  pm <- cbind(ntp, otp, gtp)
  data.table::fwrite(as.data.frame(signif(pm, 7)),
                     file.path(cgwasenv$.CGWAS_COLDATA_PATH, "Min.p"),
                     sep=" ", na="NA", row.names=F, col.names=F, quote=F)

  fdrm <- cbind(cgwasenv$.FDR_SET,
                fdrf(pm[,1], Sind, cgwasenv),
                fdrf(pm[,2], Sind, cgwasenv),
                fdrf(pm[,3], Sind, cgwasenv))
  colnames(fdrm) <- c("FDR", "cCutoff", "cN", "cLocin", "cFPLocin", "eCutoff", "eN", "eLocin", "eFPLocin", "gCutoff", "gN", "gLocin", "gFPLocin")
  write.table(signif(fdrm, 7),
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "FDRsummary.txt"),
              row.names=F,quote=F,sep="\t")
  write.table(as.data.frame(signif(trtc, 7)),
              file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "MTcorrection.txt"),
              row.names=F,col.names=F,quote=F)

  stopCluster(cl)
  gc()

  MT[,2:3] <- signif(as.numeric(MT[,2:3]), 7)
  write.table(as.data.frame(MT),
              file.path(cgwasenv$.CGWAS_RESULT_PATH, "MutipleTestingCorrection.txt")
              ,row.names=F, col.names=c("", "Ind.Tes.Num", "Threshold"), quote=F)

  if(length(ntp)>2e6) {
    ssid <- c(1:2e4, seq(2e4, 2e5, 10), seq(2e5, 2e6, 100), seq(2e6, length(ntp), 1000), (length(ntp)-1000):length(ntp))
  } else{
    ssid <- c(1:2e4, seq(2e4, 2e5, 10), seq(2e5, length(ntp), 100), (length(ntp)-1000):length(ntp))
  }

  xx <- -log10(ppoints(length(ntp))[ssid])
  yyc <- -log10(ntp[order(ntp)][ssid])
  yye <- -log10(otp[order(otp)][ssid])
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

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWASminpQQEbICoW.jpg"),
       width=3000, height=3000, res=600)
  par(mar=c(3.5, 3.5, 1, 1))
  plot(NA, xlim=c(0, max(xx)), ylim=c(0, max(c(yyc[1], yyg[1]))),
       xlab=expression(paste("Expected   ", -Log[10](italic(p)))),
       ylab=expression(paste("Observed   ", -Log[10](italic(p)))),
       mgp=c(2, 0.7, 0), las=1, cex.axis=0.85, tck=-0.015)
  points(xx, yyg, pch=20, col="#1161C9", cex=0.6)
  points(xx, yye, pch=20, col="#C0CD28", cex=0.6)
  points(xx, yyc, pch=20, col="#FD9001", cex=0.6)
  abline(c(0, 1), lwd=1)
  dev.off()
}

# Summary
step5 <- function(cgwasenv) {
  cl <- makeCluster(cgwasenv$.PARAL_NUM)
  registerDoParallel(cl)
  # globalVariables(c('i'))
  i <- 1 # assign parallel control variants

  Sind <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SnpIndex"),
                                          header = T,
                                          nThread = cgwasenv$.PARAL_NUM))
  newtick <- manpos(Sind)
  Sind <- cbind(Sind, newtick[[1]])
  newtick <- newtick[[2]]

  if(cgwasenv$.MAF_FILE_EXIST) {
    mafv <- as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "MAF"),
                                            header = F, stringsAsFactors = F,
                                            nThread = cgwasenv$.PARAL_NUM))[,1]
  }

  pm <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "Min.p"),
                                    header = F,
                                    nThread = cgwasenv$.PARAL_NUM))

  # fdrv <- apply(pm, 2, ffdrv, fv=fdr, Sind)

  fdrv <- rep(1e-6,3)

  ns <- which(pm[,1]<=fdrv[1])
  ss <- which(pm[,2]<=fdrv[2])
  os <- which(pm[,3]<=fdrv[3])

  lcidn <- locisearch(pm[ns, 1], Sind[ns, 1:2], fdrv[1], cgwasenv)
  lcids <- locisearch(pm[ss, 2], Sind[ss, 1:2], fdrv[2], cgwasenv)
  lcido <- locisearch(pm[os, 3], Sind[os, 1:2], fdrv[3], cgwasenv)

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "EBICO-GWAS.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm1 <- manhattan(pm[,3], pm[,2], os, ss, lcido, lcids, fdrv[3], fdrv[2], Sind, newtick)
  dev.off()

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWAS-EBICO.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm2 <- manhattan(pm[,2], pm[,1], ss, ns, lcids, lcidn, fdrv[2], fdrv[1], Sind, newtick)
  dev.off()

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWAS-GWAS.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm3 <- manhattan(pm[,3], pm[,1], os, ns, lcido, lcidn, fdrv[3], fdrv[1], Sind, newtick)
  dev.off()

  ressum <- rbind(c(length(os), length(os)-length(intersect(os, ss)),
                    rep(length(intersect(os, ss)), 2),
                    length(ss)-length(intersect(os, ss)),
                    length(ss)), c(sum(lcsm1[,1]), as.vector(lcsm1), sum(lcsm1[,2])),
                  c(length(ss), length(ss)-length(intersect(ss, ns)),
                    rep(length(intersect(ss, ns)), 2),
                    length(ns)-length(intersect(ss, ns)), length(ns)),
                  c(sum(lcsm2[,1]), as.vector(lcsm2), sum(lcsm2[,2])),
                  c(length(os), length(os)-length(intersect(os, ns)),
                    rep(length(intersect(os, ns)), 2),
                    length(ns)-length(intersect(os, ns)), length(ns)),
                  c(sum(lcsm3[,1]), as.vector(lcsm3), sum(lcsm3[,2])))

  tableid <- union(union(ns, ss), os)
  tableid <- tableid[order(tableid)]
  temnp <- pm[ns,1]
  temsp <- pm[ss,2]
  temop <- pm[os,3]
  rkn <- tlcidn <- rep(NA, length(ns))
  rks <- tlcids <- rep(NA, length(ss))
  rko <- tlcido <- rep(NA, length(os))
  if(length(lcidn)!=0) {
    for(j in 1:length(lcidn)) {
      tlcidn[lcidn[[j]]] <- j
      rkn[lcidn[[j]]] <- rank(temnp[lcidn[[j]]], ties.method="random")
    }
  }
  if(length(lcids)!=0) {
    for(j in 1:length(lcids)) {
      tlcids[lcids[[j]]] <- j
      rks[lcids[[j]]] <- rank(temsp[lcids[[j]]], ties.method="random")
    }
  }
  if(length(lcido)!=0) {
    for(j in 1:length(lcido)) {
      tlcido[lcido[[j]]] <- j
      rko[lcido[[j]]] <- rank(temop[lcido[[j]]], ties.method="random")
    }
  }
  if(length(tableid)!=1) {
    nm <- rbind(cbind(tlcidn, rkn, rep(1, length(rkn))), matrix(NA, length(tableid)-length(ns), 3))[order(c(ns, setdiff(tableid, ns))),]
    sm <- rbind(cbind(tlcids, rks, rep(1, length(rks))), matrix(NA, length(tableid)-length(ss), 3))[order(c(ss, setdiff(tableid, ss))),]
    om <- rbind(cbind(tlcido, rko, rep(1, length(rko))), matrix(NA, length(tableid)-length(os), 3))[order(c(os, setdiff(tableid, os))),]
  } else{
    nm <- matrix(rbind(cbind(tlcidn, rkn, rep(1, length(rkn))), matrix(NA, length(tableid)-length(ns), 3))[order(c(ns, setdiff(tableid, ns))),], nrow=length(tableid))
    sm <- matrix(rbind(cbind(tlcids, rks, rep(1, length(rks))), matrix(NA, length(tableid)-length(ss), 3))[order(c(ss, setdiff(tableid, ss))),], nrow=length(tableid))
    om <- matrix(rbind(cbind(tlcido, rko, rep(1, length(rko))), matrix(NA, length(tableid)-length(os), 3))[order(c(os, setdiff(tableid, os))),], nrow=length(tableid))
  }
  nm <- cbind(signif(pm[tableid, 1], 7), nm)
  sm <- cbind(signif(pm[tableid, 2], 7), sm)
  om <- cbind(signif(pm[tableid, 3], 7), om)
  nm[nm[,1]<=5e-8,4] <- 2
  sm[sm[,1]<=5e-8,4] <- 2
  om[om[,1]<=5e-8,4] <- 2
  if (cgwasenv$.MAF_FILE_EXIST) {
    sumtable <- cbind(Sind[tableid, 1:3], tableid, mafv[tableid], nm, sm, om)
    colnames(sumtable) <- c("CHR", "BP", "SNP", "ID", "MAF", "CP", "Clc", "Crk", "Csf", "EP", "Elc", "Erk", "Esf", "GP", "Glc", "Grk", "Gsf")
  } else {
    sumtable <- cbind(Sind[tableid, 1:3], tableid, nm, sm, om)
    colnames(sumtable) <- c("CHR", "BP", "SNP", "ID", "CP", "Clc", "Crk", "Csf", "EP", "Elc", "Erk", "Esf", "GP", "Glc", "Grk", "Gsf")
  }
  write.csv(sumtable,
            file.path(cgwasenv$.CGWAS_RESULT_PATH, "Summaryhits.csv"),
            row.names=F,quote=F)

  ns <- which(pm[,1]<=5e-8)
  ss <- which(pm[,2]<=5e-8)
  os <- which(pm[,3]<=5e-8)

  lcidn <- locisearch(pm[ns, 1], Sind[ns, 1:2], 5e-8, cgwasenv)
  lcids <- locisearch(pm[ss, 2], Sind[ss, 1:2], 5e-8, cgwasenv)
  lcido <- locisearch(pm[os, 3], Sind[os, 1:2], 5e-8, cgwasenv)

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "EBICO-GWAS-S.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm1 <- manhattan(pm[,3], pm[,2], os, ss, lcido, lcids, fdrv[3], fdrv[2], Sind, newtick)
  dev.off()

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWAS-EBICO-S.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm2 <- manhattan(pm[,2], pm[,1], ss, ns, lcids, lcidn, fdrv[2], fdrv[1], Sind, newtick)
  dev.off()

  jpeg(file.path(cgwasenv$.CGWAS_RESULT_PATH, "CGWAS-GWAS-S.jpg"),
       width=6000,height=4500,res=600)
  par(mar=c(1, 4, 1, 2))
  lcsm3 <- manhattan(pm[,3], pm[,1], os, ns, lcido, lcidn, fdrv[3], fdrv[1], Sind, newtick)
  dev.off()

  ressum <- rbind(ressum, c(length(os), length(os)-length(intersect(os, ss)), rep(length(intersect(os, ss)), 2), length(ss)-length(intersect(os, ss)), length(ss)), c(sum(lcsm1[,1]), as.vector(lcsm1), sum(lcsm1[,2])), c(length(ss), length(ss)-length(intersect(ss, ns)), rep(length(intersect(ss, ns)), 2), length(ns)-length(intersect(ss, ns)), length(ns)), c(sum(lcsm2[,1]), as.vector(lcsm2), sum(lcsm2[,2])), c(length(os), length(os)-length(intersect(os, ns)), rep(length(intersect(os, ns)), 2), length(ns)-length(intersect(os, ns)), length(ns)), c(sum(lcsm3[,1]), as.vector(lcsm3), sum(lcsm3[,2])))

  colnames(ressum) <- c("T1all", "T1unique", "T1T2overlap", "T2T1overlap", "T2unique", "T2all")
  rownames(ressum) <- c("G-E-FDR-SNP", "G-E-FDR-LN", "E-C-FDR-SNP", "E-C-FDR-LN", "G-C-FDR-SNP", "G-C-FDR-LN", "G-E-SBC-SNP", "G-E-SBC-LN", "E-C-SBC-SNP", "E-C-SBC-LN", "G-C-SBC-SNP", "G-C-SBC-LN")
  write.csv(ressum,
            file.path(cgwasenv$.CGWAS_RESULT_PATH, "Summary.csv"),
            quote=F)


  localsep <- 1e5

  ACstatm <- read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "ACEffect.txt"),
                        header=T,stringsAsFactors=F)
  corm.tmp <- read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "BCCorrelationStat.txt"),
                         header=T, stringsAsFactors=F)
  statcorm <- gcrom(corm.tmp[,5],cgwasenv$.TRAIT_NUM)
  bgcorm <- gcrom(corm.tmp[,6],cgwasenv$.TRAIT_NUM)
  logm <- read.table(file.path(cgwasenv$.CGWAS_TEMPDATA_PATH, "EBICOlog.txt"),
                     header=T,stringsAsFactors=F)

  swpm <- as.matrix(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, "SWaT.p"),
                                      header=F,
                                      nThread = cgwasenv$.PARAL_NUM))
  swtable <- cbind(Sind[tableid, 1:3], tableid, swpm[tableid,])
  colnames(swtable) <- c("CHR", "BP", "SNP", "ID", paste0("SWp", 1:ncol(swpm)))
  write.csv(swtable,
            file.path(cgwasenv$.CGWAS_RESULT_PATH, "SWaThits.csv"),
            row.names=F,quote=F)

  compn <- as.character(ACstatm[,1])
  bm <- foreach(i=compn, .combine="cbind", .inorder=T) %dopar%
          as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(i, ".beta")),
                                          header = F, stringsAsFactors = F, nThread = 1))[,1]
  tm <- foreach(i=compn, .combine="cbind", .inorder=T) %dopar%
          as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(i, ".efstat")),
                                          header = F, stringsAsFactors = F, nThread = 1))[,1]
  ebtable <- cbind(Sind[tableid, 1:3], tableid, bm[tableid,], tm[tableid,])
  colnames(ebtable) <- c("CHR", "BP", "SNP", "ID", paste0(compn, "_beta"), paste0(compn, "_stat"))
  write.csv(ebtable, file.path(cgwasenv$.CGWAS_RESULT_PATH, "EbICohits.csv"),
            row.names=F, quote=F)

  bm <- foreach(i=cgwasenv$.TRAIT_NAME, .combine="cbind", .inorder=T) %dopar%
          as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(i, ".beta")),
                                          header = F, stringsAsFactors = F, nThread = 1))[,1]
  tm <- foreach(i=cgwasenv$.TRAIT_NAME, .combine="cbind", .inorder=T) %dopar%
          as.data.frame(data.table::fread(file.path(cgwasenv$.CGWAS_COLDATA_PATH, paste0(i, ".efstat")),
                                        header = F, stringsAsFactors = F, nThread = 1))[,1]
  gtable <- cbind(Sind[tableid, 1:3], tableid, bm[tableid,], tm[tableid,])
  colnames(gtable) <- c("CHR", "BP", "SNP", "ID", paste0(cgwasenv$.TRAIT_NAME, "_beta"), paste0(cgwasenv$.TRAIT_NAME, "_stat"))
  write.csv(gtable, file.path(cgwasenv$.CGWAS_RESULT_PATH, "GWAShits.csv"), row.names=F, quote=F)
  stopCluster(cl)
}
