step1 <- function(ST_Datapath, ST_Filename, Clv, Bdd, Fdd){
  for(i in 1:cgwasenv$.pheNum){
    bpm <- as.data.frame(data.table::fread(paste0(ST_Datapath,ST_Filename[i]),header=T))
    if (i == 1){
      RsnpnumTmp <- nrow(bpm)
      tm <- matrix(0,RsnpnumTmp,cgwasenv$.pheNum)
      snpm <- bpm[, 1:3]
    }
    if(Fdd){
      sig <- 2*as.numeric(nchar(bpm[,4])==(Bdd+2))-1
      tm[,i] <- -qnorm(as.numeric(bpm[,5])/2)*sig
    } else{
      tm[,i] <- -qnorm(as.numeric(bpm[,5])/2)*abs(as.numeric(bpm[,4]))/as.numeric(bpm[,4])
    }
    print(paste0(i," Trait Statistics Generated"))
  }

  snpInx = apply(tm, 1,function(x) sum(is.na(x)) == 0)
  tm <- tm[snpInx, ]
  snpm <- snpm[snpInx, ]

  data.table::fwrite(snpm, file.path(cgwasenv$.cgwas_result, "SnpIndex"), quote = F, col.names = T, row.names = F)
  print("SnpIndex Completed")

  cgwasenv$.Rsnpnum <- nrow(snpm)
  if (is.na(cgwasenv$.Gt)) cgwasenv$.Gt <- max(5e-8, 0.05/cgwasenv$.Rsnpnum)
  if (is.na(cgwasenv$.St)) cgwasenv$.St <- max(1.6e-6, 1/cgwasenv$.Rsnpnum)
  cgwasenv$.Esnpnum <- round(0.05/cgwasenv$.Gt)

  if (is.na(Clv)){
    Correctlambda <- T
  } else {
    Correctlambda <- F
  }


  if(Correctlambda){
    Sind <- as.numeric(snpm[,2])
    clv <- c()
    cxv <- c()
    gif <- c()
  }
  lv <- c()
  xv <- c()
  LDr <- 500000
  preset <- 5000
  print("")
  for(i in 1:cgwasenv$.pheNum){
    abst <- abs(tm[,i])
    lv <- c(lv,quantile(abst,0.5)/qnorm(0.75))
    xv <- c(xv,mean(abst^2))
    print(paste0(cgwasenv$.traitName[i],"   lambda: ",round(lv[i],4),"   meanX2: ",round(xv[i],4)))
    if(Correctlambda){
      tarid <- which(abst>abs(qnorm(1e-4/2)))
      tarid <- tarid[which(tarid>preset)]
      tarid <- tarid[which(tarid<(cgwasenv$.Rsnpnum-preset))]
      id <- c()
      tid <- c()
      bid <- 1
      doub <- 1
      if (length(tarid) >= 2){
        for (j in 2:length(tarid)){
          if (tarid[j]>tail(bid,1)){
            if (doub==1){
              tbid <- Sind[(tarid[j-1]-preset):(tarid[j-1]+preset)]
              bid <- which((tbid<(Sind[tarid[j-1]]+LDr))&(tbid>(Sind[tarid[j-1]]-LDr)))+tarid[j-1]-preset-1
              id <- unique(c(id,bid))
              doub <- 0
            }
            tbid <- Sind[(tarid[j]-preset):(tarid[j]+preset)]
            bid <- which((tbid<(Sind[tarid[j]]+LDr))&(tbid>(Sind[tarid[j]]-LDr)))+tarid[j]-preset-1
            id <- unique(c(id,bid))
            if(length(id)>LDr){
              tid <- unique(c(tid,id))
              id <- c()
            }
          } else{
            doub <- 1
          }
        }
      }
      tid <- unique(c(tid,id))
      gif <- c(gif,quantile(abst[setdiff(1:cgwasenv$.Rsnpnum,tid)],0.5)/qnorm(0.75))
      print(paste0("Inflation factor: ",round(gif[i],4)))
      if(gif[i]>=1){
        tm[,i] <- tm[,i]/gif[i]
        clv <- c(clv,quantile((abst/gif[i]),0.5)/qnorm(0.75))
        cxv <- c(cxv,mean((abst/gif[i])^2))
      } else{
        gif[i] <- 1
        clv <- c(clv,lv[i])
        cxv <- c(cxv,xv[i])
      }
      print(paste0(cgwasenv$.traitName[i],"   corrected lambda: ",round(clv[i],4),"   meanX2: ",round(cxv[i],4)))
    }
    print("")
  }
  if(Correctlambda){
    data.table::fwrite(cbind(cgwasenv$.traitName,lv,xv,gif,clv,cxv),file.path(cgwasenv$.cgwas_result, "GwasDataCheck.txt"),row.names = F,col.names = F,quote = F)
    data.table::fwrite(as.data.frame(rep(1,cgwasenv$.pheNum)), file.path(cgwasenv$.cgwas_simuldata, "lambda.txt"), row.names = F, col.names = F, quote = F)
  } else{
    data.table::fwrite(cbind(cgwasenv$.traitName,lv,xv),file.path(cgwasenv$.cgwas_result, "GwasDataCheck.txt"),row.names = F,col.names = F,quote = F)
    data.table::fwrite(as.data.frame(Clv), file.path(cgwasenv$.cgwas_simuldata, "lambda.txt"), row.names = F, col.names = F, quote = F)
  }

  data.table::fwrite(as.data.frame(tm),file.path(cgwasenv$.cgwas_rowdata, "All.s"),col.names=F,row.names=F,sep=" ")
  print(paste0("All Data Block Saved"))

  print(paste0("Calculate Correlation Matrix.."))
  orgtemp <- qnorm(ppoints(1:cgwasenv$.Rsnpnum))
  for(i in 1:cgwasenv$.pheNum){
    tm[order(tm[,i]),i] <- orgtemp
    print(paste0(i," Z score Transformed"))
  }
  ctm <- cor(tm)
  colnames(ctm) <- cgwasenv$.traitName
  write.csv(ctm,file.path(cgwasenv$.cgwas_result, paste0(cgwasenv$.pheNum, "StatsCor.csv")),row.names=F,quote=F)
  print(paste0("Correlation Matrix Completed"))

  print(paste0("Step1 Completed"))
}

step2_parallel <- function(){
    lambda <- as.numeric(read.table(file.path(cgwasenv$.cgwas_simuldata, "lambda.txt"), header = F)[,1])
    ocorm <- as.matrix(read.csv(file.path(cgwasenv$.cgwas_result, paste0(cgwasenv$.pheNum, "StatsCor.csv")),header=T))
    preCorrection()
    qvLen <- length(cgwasenv$.qv)
    fminm <- matrix(0,ncol=cgwasenv$.pheNum,nrow=cgwasenv$.simulNum*qvLen)
    ocorm2 <- ocorm^2
    ocorm3 <- ocorm^3
    ocorm4 <- ocorm^4

    step2Time1 <- Sys.time()

    snowfall::sfSetMaxCPUs(parallel::detectCores() - 1)
    if (cgwasenv$.paral == 1){
      snowfall::sfInit(parallel = F)
    } else {
      snowfall::sfInit(parallel = T, cpus = cgwasenv$.paral)
    }
    snowfall::sfLibrary(MASS)
    pheNum <- cgwasenv$.pheNum
    snowfall::sfExport("ocorm", "ocorm2", "ocorm3", "ocorm4", "pheNum", "lambda")
    snowfall::sfClusterSetupRNGstream(seed=2)

    for(i in 1:cgwasenv$.simulNum){
      smltTime1 <- Sys.time()

      fillind <- (qvLen*(i-1)+1):(qvLen*i)

      c3aTime1 <- Sys.time()
      resm = matrix(unlist(snowfall::sfLapply(1:cgwasenv$.Esnpnum, mkdf_corterm3all_parallel)), byrow = T, ncol = cgwasenv$.pheNum)
      c3aTime2 <- Sys.time()
      print(paste0(i, "th c3a use time:"))
      print(c3aTime2 - c3aTime1)

      fminm[fillind,] <- apply(resm,2,function(a){return(sort(a))})[cgwasenv$.qv,]

      smltTime2 <- Sys.time()
      print(paste0(i, "th simulation use time:"))
      print(smltTime2 - smltTime1)
    }

    snowfall::sfStop()
    step2Time2 <- Sys.time()
    print(paste0("step2 use time:"))
    print(step2Time2 - step2Time1)

    for(i in 1:cgwasenv$.pheNum){
      data.table::fwrite(matrix(fminm[,i],ncol=qvLen,byrow = T),file.path(cgwasenv$.cgwas_allp,paste0("qm",cgwasenv$.simulNum,"_",i)),row.names=F,col.names=F,quote=F)
    }

    print(paste0("Step2 Completed"))
}

step3_parallel <- function(SimuReg){
  data.table::fwrite(as.data.frame((qbeta(0.05,1:cgwasenv$.Esnpnum,cgwasenv$.Esnpnum:1)/(((qbeta(0.05,1,cgwasenv$.Esnpnum)/cgwasenv$.Gt-1)*((cgwasenv$.Esnpnum-1):0)+cgwasenv$.Esnpnum-1)/(cgwasenv$.Esnpnum-1)))[cgwasenv$.qv]),file.path(cgwasenv$.cgwas_minp, cgwasenv$.pheNum+1),row.names=F,col.names=F,quote=F)
  basecf <- as.numeric(read.table(file.path(cgwasenv$.cgwas_minp, cgwasenv$.pheNum+1),header=F,stringsAsFactors=F)[,1])

  for(i in 1:cgwasenv$.pheNum){
    qm <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_allp,paste0("qm",cgwasenv$.simulNum,"_",i)),header=F))
    curcf <- apply(qm,2,quantile,probs=0.05)
    data.table::fwrite(as.data.frame(curcf),file.path(cgwasenv$.cgwas_minp, i),row.names=F,col.names=F,quote=F)

    jpeg(file.path(cgwasenv$.cgwas_correction, paste0(i, ".jpg")), width=3600, height=1350, res=200)
    par(mfrow=c(2,2))
    simulQQplot(qm,"Original",cgwasenv$.simulNum,cgwasenv$.Esnpnum)
    simulQQplot(t(t(qm)*basecf/curcf),"Corrected",cgwasenv$.simulNum,cgwasenv$.Esnpnum)
    simultime(t(t(qm)*basecf/curcf),cgwasenv$.simulNum,cgwasenv$.Esnpnum,SimuReg,"Simulate FDR")
    inflationplot(basecf/curcf,cgwasenv$.Esnpnum,"Inflation Coefficient")
    dev.off()
  }

  qm <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_allp,paste0("qm",cgwasenv$.simulNum,"_",1)),header=F))
  basecf <- as.numeric(read.table(file.path(cgwasenv$.cgwas_minp, cgwasenv$.pheNum+1),header=F,stringsAsFactors=F)[,1])
  curcf <- as.numeric(read.table(file.path(cgwasenv$.cgwas_minp, 1),header=F,stringsAsFactors=F)[,1])
  aa <- t(t(qm)*basecf/curcf)
  print("1 Quantile Data Extracted")
  for(i in 2:cgwasenv$.pheNum){
    qm <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_allp,paste0("qm",cgwasenv$.simulNum,"_",i)),header=F))
    curcf <- as.numeric(read.table(file.path(cgwasenv$.cgwas_minp, i),header=F,stringsAsFactors=F)[,1])
    bb <- t(t(qm)*basecf/curcf)
    aa[bb<aa] <- bb[bb<aa]
    print(paste0(i," Quantile Data Extracted"))
  }

  qm <- aa
  curcf <- apply(qm,2,quantile,probs=0.05)
  data.table::fwrite(as.data.frame(curcf),file.path(cgwasenv$.cgwas_minp, cgwasenv$.pheNum+2),row.names=F,col.names=F,quote=F)
  jpeg(file.path(cgwasenv$.cgwas_correction, paste0(cgwasenv$.pheNum+1, ".jpg")), width=2400, height=2600, res=200)
  par(mfrow=c(2,2))
  simulQQplot(qm,"Original",cgwasenv$.simulNum,cgwasenv$.Esnpnum)
  simulQQplot(t(t(qm)*basecf/curcf),"Corrected",cgwasenv$.simulNum,cgwasenv$.Esnpnum)
  simultime(t(t(qm)*basecf/curcf),cgwasenv$.simulNum,cgwasenv$.Esnpnum,SimuReg,"Simulate FDR")
  inflationplot(basecf/curcf,cgwasenv$.Esnpnum,"Inflation Coefficient")
  dev.off()

  print(paste0("Step3 Completed"))
}

step4_parallel <- function(){
  step4Time1 <- Sys.time()
  ocorm <- as.matrix(read.csv(file.path(cgwasenv$.cgwas_result, paste0(cgwasenv$.pheNum, "StatsCor.csv")),header=T))
  dm <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_rowdata, "All.s"),header=F))
  pm <- pnorm(-abs(dm))*2
  ordm <- matrix(0,ncol=cgwasenv$.pheNum,cgwasenv$.Rsnpnum)
  ocorm2 <- ocorm^2
  ocorm3 <- ocorm^3
  ocorm4 <- ocorm^4

  snowfall::sfSetMaxCPUs(parallel::detectCores() - 1)
  if (cgwasenv$.paral == 1){
    snowfall::sfInit(parallel = F)
  } else {
    snowfall::sfInit(parallel = T, cpus = cgwasenv$.paral)
  }
  pheNum <- cgwasenv$.pheNum
  rowdataPath <- cgwasenv$.cgwas_rowdata
  snowfall::sfExport("ocorm", "ocorm2", "ocorm3", "ocorm4", "pheNum", "rowdataPath")
  paralInx <- rep(1:cgwasenv$.paral, each = ceiling(cgwasenv$.Rsnpnum / cgwasenv$.paral))[1:cgwasenv$.Rsnpnum]
  for (i in 1: cgwasenv$.paral){
    data.table::fwrite(as.data.frame(dm[paralInx == i,]), file.path(cgwasenv$.cgwas_rowdata, paste0("All.s.", i)),row.names=F,col.names=F,quote=F)
    data.table::fwrite(as.data.frame(pm[paralInx == i,]), file.path(cgwasenv$.cgwas_rowdata, paste0("All.sp.", i)),row.names=F,col.names=F,quote=F)
  }

  resm <- matrix(unlist(snowfall::sfLapply(1:cgwasenv$.paral, corterm3all_parallel)), byrow = T, ncol = cgwasenv$.pheNum)
  data.table::fwrite(as.data.frame(signif(resm,8)), file.path(cgwasenv$.cgwas_result, "All.res"),row.names=F,col.names=F,quote=F)
  rm(resm)
  gc()

  snowfall::sfStop()

  stepw <- matrix(unlist(lapply(1:cgwasenv$.Rsnpnum, function(i) order(pm[i, ]))), byrow = T, ncol = cgwasenv$.pheNum)
  data.table::fwrite(as.data.frame(stepw), file.path(cgwasenv$.cgwas_result, "All.labspInx"),row.names=F,col.names=F,quote=F)
  rm(stepw)
  gc()
  ordm <- matrix(unlist(lapply(1:cgwasenv$.Rsnpnum, function(i) sort(pm[i, ]))), byrow = T, ncol = cgwasenv$.pheNum)
  data.table::fwrite(as.data.frame(signif(ordm,8)), file.path(cgwasenv$.cgwas_result, "All.sortsp"),row.names=F,col.names=F,quote=F)
  rm(ordm)
  data.table::fwrite(as.data.frame(signif(pm,4)), file.path(cgwasenv$.cgwas_result, "All.sp"),row.names=F,col.names=F,quote=F)
  rm(pm)
  gc()

  step4Time2 <- Sys.time()
  print(paste0("step4 use time:"))
  print(step4Time2 - step4Time1)
  print(paste0("Step4 Completed"))
}

step5 <- function(Ht){
  step5Time1 <- Sys.time()
  pm <- data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "All.res"),header=F))
  print("CGWAS Data Extracted")

  insertm1 <- insertm2 <- matrix(0,nrow=cgwasenv$.Rsnpnum,2)
  simuper <- (cgwasenv$.qv-1)*(cgwasenv$.Rsnpnum-1)/(round(0.05/cgwasenv$.Gt)-1)+1
  simuperf <- floor(simuper)
  simuperc <- ceiling(simuper)
  simuperc[which(simuper==simuperc)] <- simuperc[which(simuper==simuperc)]+1
  simuperc[1] <- 1
  repv <- simuperf[-1]-simuperc[-length(simuperc)]+1
  for(i in 1:(length(simuperc)-1)){
    insertm1[simuperc[i]:(simuperc[i+1]-1),2] <- (1:repv[i])-1
  }
  temre <- simuper[-1]-simuper[-length(simuper)]
  insertm1[,2] <- (insertm1[,2]+rep((simuperc-simuper)[-length(simuperc)],repv))/rep(temre,repv)
  insertm1[,1] <- 1-insertm1[,2]

  basecf <- as.numeric(read.table(file.path(cgwasenv$.cgwas_minp, cgwasenv$.pheNum+1),header=F,stringsAsFactors=F)[,1])
  tind <- scalev <- rep(0,cgwasenv$.Rsnpnum)
  for(i in c(1:cgwasenv$.pheNum,cgwasenv$.pheNum+2)){
    curcf <- as.numeric(read.table(file.path(cgwasenv$.cgwas_minp, i),header=F,stringsAsFactors=F)[,1])
    baseper <- basecf/curcf
    insertm2[,1] <- rep(baseper[-length(baseper)],repv)
    insertm2[,2] <- rep(baseper[-1],repv)
    insertm2 <- insertm1*insertm2
    scalev <- insertm2[,1]+insertm2[,2]
    if(i==(cgwasenv$.pheNum+2)){
      cgsv <- scalev
    } else{
      tind <- order(pm[,i])
      pm[tind,i] <- pm[tind,i]*scalev
      print(paste0(i," CGWAS Data Corrected"))
    }
    if(i==1){
      sgsv <- scalev
    }
  }

  pm[pm>1] <- 1
  mincp <- apply(pm,1,min)
  data.table::fwrite(as.data.frame(signif(mincp,4)),file.path(cgwasenv$.cgwas_result, "Min.cp"),row.names=F,col.names=F,quote=F)
  print(paste0("Derived CGWAS Minp"))
  tind <- order(mincp)
  mincp[tind] <- mincp[tind]*cgsv
  mincp[mincp>1] <- 1
  data.table::fwrite(as.data.frame(signif(mincp,4)),file.path(cgwasenv$.cgwas_result, "Min.corrcp"),row.names=F,col.names=F,quote=F)
  print(paste0("Derived CGWAS corrMinp"))

  pm <- signif(pm,4)

  sugghitind2 <- which(mincp<cgwasenv$.St)
  HJind2 <- which(mincp<Ht)

  spo <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_result, "All.sortsp"),header=F))
  file.remove(file.path(cgwasenv$.cgwas_result, "All.sortsp"))

  minsp <- spo[,1]
  tind <- order(minsp)
  minsp[tind] <- minsp[tind]*sgsv
  minsp[minsp>1] <- 1
  data.table::fwrite(as.data.frame(signif(minsp,4)),file.path(cgwasenv$.cgwas_result, "Min.corrsp"),row.names=F,col.names=F,quote=F)
  print(paste0("Derived corrGWAS Minp"))

  spo <- signif(spo,4)
  data.table::fwrite(as.data.frame(spo),file.path(cgwasenv$.cgwas_result, "All.sortsp"),row.names=F,col.names=F,quote=F,sep=" ")
  print(paste0("sortGWAS Result Saved"))
  data.table::fwrite(as.data.frame(spo[,1]),file.path(cgwasenv$.cgwas_result, "Min.sp"),row.names=F,col.names=F,quote=F)
  print(paste0("Derived GWAS Minp"))

  sugghitind <- which(as.numeric(minsp)<cgwasenv$.St)
  HJind <- which(minsp<Ht)

  sugghit3 <- union(sugghitind,sugghitind2)
  HJ3 <- union(HJind,HJind2)
  sugghit <- matrix(0,length(sugghit3)*cgwasenv$.pheNum,3)
  HJM <- matrix(0,length(HJ3),2)
  sugghit[,1] <- rep(sugghit3,each=cgwasenv$.pheNum)
  sugghit[,2] <- as.vector(t(pm[sugghit3,]))
  sugghit[,3] <- rep(c(1:cgwasenv$.pheNum),length(sugghit3))
  HJM[,1] <- as.numeric(mincp)[HJ3]
  HJM[,2] <- as.numeric(minsp)[HJ3]
  HJpm <- pm[HJ3,]
  data.table::fwrite(as.data.frame(pm),file.path(cgwasenv$.cgwas_result, "All.cp"),row.names=F,col.names=F,quote=F,sep=" ")
  print(paste0("CGWAS Result Saved"))
  rm(pm)
  gc()

  Sind <- as.data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "SnpIndex"),header=T))
  temres <- cbind(Sind[as.numeric(sugghit[,1]),],sugghit)
  temHJ <- cbind(Sind[HJ3,],HJM)
  temres <- temres[order(temres[,5]),]
  temres <- temres[order(temres[,2]),]
  temres <- temres[order(temres[,1]),]
  temHJ <- temHJ[order(temHJ[,2]),]
  temHJ <- temHJ[order(temHJ[,1]),]
  colnames(temres)[4:6] <- c("SNPid","CorrCP","CombNum")
  colnames(temHJ)[4:5] <- c("CP","BP")
  write.csv(temres,file.path(cgwasenv$.cgwas_result, "CgwasRawSuggHits.csv"),row.names=F,quote=F)
  data.table::fwrite(as.data.frame(temHJ),file.path(cgwasenv$.cgwas_result, "HJ_man.txt"),row.names=F,col.names=F,quote=F)
  print(paste0("CGWAS Suggestive Hits Saved"))

  sspo <- spo[sugghit3,]
  HJspo <- spo[HJ3,]
  rm(spo)
  gc()
  spsignum <- c()
  cpsignum <- c()
  for(i in 1:nrow(HJspo)){
    spsignum <- c(spsignum,sum(HJspo[i,]<=cgwasenv$.St))
    cpsignum <- c(cpsignum,which(HJpm[i,]==min(HJpm[i,]))[1])
  }
  spl <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_result, "All.labspInx"),header=F))
  sspl <- matrix(unlist(lapply(sugghit3,function(i) {cgwasenv$.traitName[spl[i,]]})), byrow = T, ncol = ncol(spl))
  HJspl <- matrix(unlist(lapply(HJ3,function(i) {cgwasenv$.traitName[spl[i,]]})), byrow = T, ncol = ncol(spl))
  rm(spl)
  gc()
  sugghit <- matrix(0,length(sugghit3)*cgwasenv$.pheNum,3)
  sugghit[,1] <- rep(sugghit3,each=cgwasenv$.pheNum)
  sugghit[,2] <- as.vector(t(sspo[,]))
  sugghit[,3] <- as.vector(t(sspl[,]))

  temres <- cbind(Sind[as.numeric(sugghit[,1]),],sugghit)
  temHJall <- cbind(Sind[HJ3,],cpsignum,spsignum,HJpm,HJspl,HJspo)
  temres <- temres[order(temres[,2]),]
  temres <- temres[order(temres[,1]),]
  temHJall <- temHJall[order(temHJall[,2]),]
  temHJall <- temHJall[order(temHJall[,1]),]
  temHJall <- temHJall[,-c(1:2)]
  colnames(temres)[4:6] <- c("SNPid","SP","Trait")
  write.csv(as.matrix(temres),file.path(cgwasenv$.cgwas_result, "GwasRawSuggHits.csv"),row.names=F,quote=F)
  data.table::fwrite(as.data.frame(temHJall),file.path(cgwasenv$.cgwas_result, "HJ_bar.txt"),row.names=F,col.names=F,quote=F)
  print(paste0("GWAS Suggestive Hits Saved"))

  print(object.size(sspl), units = 'auto')
  print(object.size(HJspl), units = 'auto')
  print(object.size(sugghit), units = 'auto')
  print(object.size(Sind), units = 'auto')
  print(object.size(temres), units = 'auto')
  print(object.size(temHJall), units = 'auto')

  step5Time2 <- Sys.time()
  print(paste0("step5 use time:"))
  print(step5Time2 - step5Time1)
  print(paste0("Step5 Completed"))
}

step6 <- function(m_wid=3500,m_hei=3500,sing_col=c('chartreuse2','chartreuse4'),comb_col=c('darkgoldenrod2','darkgoldenrod4'),speeduppoints=T,q_wid=3500,q_hei=2400){
  step6Time1 <- Sys.time()
  print(paste0("Read Data.."))
  Sind <- as.data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "SnpIndex"),header=T))
  singpm <- as.data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "Min.corrsp"),header=F))
  cpm <- as.data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "Min.corrcp"),header=F))
  cpfm <- cbind(Sind,-log10(singpm),-log10(cpm))

  if(speeduppoints){
    cpfm <- cpfm[-which((cpfm[,4]<1)&(cpfm[,5]<1)),]
    tem <- which((cpfm[,4]<2)&(cpfm[,5]<2))
    cpfm <- cpfm[-sample(tem,round(length(tem)*0.75)),]
    tem <- which((cpfm[,4]<3)&(cpfm[,5]<3))
    cpfm <- cpfm[-sample(tem,round(length(tem)*0.5)),]
    tem <- which((cpfm[,4]<4)&(cpfm[,5]<4))
    cpfm <- cpfm[-sample(tem,round(length(tem)*0.25)),]
  }
  cpfm[which(cpfm[,4]<1),4] <- 0
  cpfm[which(cpfm[,5]<1),5] <- 0

  print(paste0("Produce Manhattan Plots.."))
  jpeg(file.path(cgwasenv$.cgwas_result, "man.jpg"),width=m_wid,height=m_hei,res=300)
  par(mar=c(2,5,2,2))
  manhattanPlot(cpfm,cex.axis=1,spt.cex=1,cpt.cex=1,spt.col=sing_col,cpt.col=comb_col,pch=20)
  dev.off()
  print(paste0("Manhattan Plots Completed"))
  rm(cpfm)
  gc()

  print(paste0("Read Data For QQplots.."))
  tsingm <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_result, "All.sp"),header=F))
  tcombm <- as.matrix(data.table::fread(file.path(cgwasenv$.cgwas_result, "All.cp"),header=F))
  minqqm <- cbind(singpm,cpm)

  tsingm <- apply(tsingm,2,function(a){return(a[order(a)])})
  print(paste0("QQplots: GWAS Data Extracted"))
  tcombm <- apply(tcombm,2,function(a){return(a[order(a)])})
  print(paste0("QQplots: CGWAS Data Extracted"))

  sgif <- as.numeric(qnorm(apply(tsingm/2,2,function(a){return(mean(a[c(floor(cgwasenv$.Rsnpnum/2),ceiling(cgwasenv$.Rsnpnum/2))]))}))/qnorm(0.25))
  cgif <- as.numeric(qnorm(apply(tcombm/2,2,function(a){return(mean(a[c(floor(cgwasenv$.Rsnpnum/2),ceiling(cgwasenv$.Rsnpnum/2))]))}))/qnorm(0.25))
  print(paste0("Lambda Calculated"))

  c1 <- round(cgwasenv$.Rsnpnum*0.0001)
  cc1 <- round(cgwasenv$.Rsnpnum*0.000001)
  c2 <- round(cgwasenv$.Rsnpnum*0.001)
  cc2 <- round(cgwasenv$.Rsnpnum*0.00001)
  c3 <- round(cgwasenv$.Rsnpnum*0.01)
  c4 <- round(cgwasenv$.Rsnpnum*0.1)
  if(cgwasenv$.Rsnpnum>=1e6){
    cgwasenv$.qv <- c(1:(c1-1),seq(c1,(c2-1),cc1),seq(c2,(c3-1),cc2),seq(c3,(c4-1),c1),seq(c4,cgwasenv$.Rsnpnum-100,c2),(cgwasenv$.Rsnpnum-99):cgwasenv$.Rsnpnum)
  } else if(cgwasenv$.Rsnpnum>=1e5){
    cgwasenv$.qv <- c(1:(c2-1),seq(c2,(c3-1),cc2),seq(c3,(c4-1),c1),seq(c4,cgwasenv$.Rsnpnum-100,c2),(cgwasenv$.Rsnpnum-99):cgwasenv$.Rsnpnum)
  } else{
    cgwasenv$.qv <- 1:cgwasenv$.Rsnpnum
  }

  minqqm <- apply(minqqm,2,function(a){return(a[order(a)])})
  tsingm <- -log10(tsingm[cgwasenv$.qv,])
  tcombm <- -log10(tcombm[cgwasenv$.qv,])
  minqqm <- -log10(minqqm[cgwasenv$.qv,])
  ee <- -log10(ppoints(cgwasenv$.Rsnpnum)[cgwasenv$.qv])

  print(paste0("Produce QQ Plots.."))
  jpeg(file.path(cgwasenv$.cgwas_result, "qq.jpg"),width=q_wid,height=q_hei,res=300)
  par(mfrow=c(2,3))
  par(mar=c(5,5,2,2))
  qqPlot(tsingm,tcombm,minqqm,ee,sgif,cgif,c(sing_col[1],comb_col[1]))
  dev.off()
  print(paste0("QQ Plots Completed"))

  step6Time2 <- Sys.time()
  print(paste0("step6 use time:"))
  print(step6Time2 - step6Time1)
  print(paste0("Step6 Completed"))
}

step7 <- function(){
  step7Time1 <- Sys.time()
  basegm <- read.csv(file.path(cgwasenv$.cgwas_result, "GwasRawSuggHits.csv"),header=T,stringsAsFactors=F)
  basecm <- read.csv(file.path(cgwasenv$.cgwas_result, "CgwasRawSuggHits.csv"),header=T,stringsAsFactors=F)
  singpm <- as.data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "Min.corrsp"),header=F))
  cpm <- as.data.frame(data.table::fread(file.path(cgwasenv$.cgwas_result, "Min.corrcp"),header=F))
  print(paste0("Read Data Completed"))
  SNPn <- unique(basegm$SNP)
  tem <- c()
  for(i in 1:length(SNPn)){
    temgm <- basegm[(((i-1)*cgwasenv$.pheNum+1):(i*cgwasenv$.pheNum)),]
    temcm <- basecm[(((i-1)*cgwasenv$.pheNum+1):(i*cgwasenv$.pheNum)),]
    fdrsp <- as.numeric(temgm[,5])
    fdrcp <- as.numeric(temcm[,5])
    sgt <- 0
    cgt <- 0
    temv <- c(temgm[1,1:5],temcm[1,5],temgm[1,6],sgt,cgt,temcm[1,6],c(temgm[,6],fdrsp,fdrcp[order(temcm[,6])])[rep(seq(1,3*cgwasenv$.pheNum,cgwasenv$.pheNum),cgwasenv$.pheNum)+rep(0:(cgwasenv$.pheNum-1),each=3)])
    tem <- rbind(tem,temv)
  }
  print(paste0("Format Data Completed"))

  locim <- matrix(unlist(tem[,1:2]),ncol=2)
  locic <- c()
  startchr <- locim[1,1]
  startbp <- locim[1,2]
  startloci <- 1
  locinum <- 1
  for(i in 1:nrow(locim)){
    if(locim[i,1]==startchr){
      if((locim[i,2]-startbp)<500000){
        locic[i] <- paste(startchr,"_",startloci,sep="")
        startbp <- locim[i,2]
      } else{
        startloci <- startloci + 1
        locinum <- locinum + 1
        locic[i] <- paste(startchr,"_",startloci,sep="")
        startbp <- locim[i,2]
      }
    } else{
      startchr <- locim[i,1]
      startloci <- 1
      locinum <- locinum + 1
      locic[i] <- paste(startchr,"_",startloci,sep="")
      startbp <- locim[i,2]
    }
  }
  print(paste0("Loci Zoning Completed"))

  tem <- cbind(tem[,1:4],locic,singpm[as.numeric(tem[,4]),1],cpm[as.numeric(tem[,4]),1],tem[,5:ncol(tem)])
  if(length(which((tem[,6]>cgwasenv$.St)&(tem[,7]>cgwasenv$.St)))!=0){
    tem <- tem[-which((tem[,6]>cgwasenv$.St)&(tem[,7]>cgwasenv$.St)),]
  }

  tem[which(tem[,6]<=cgwasenv$.St),11] <- 1
  tem[which(tem[,6]<=cgwasenv$.Gt),11] <- 2
  tem[which(tem[,7]<=cgwasenv$.St),12] <- 1
  tem[which(tem[,7]<=cgwasenv$.Gt),12] <- 2

  colnames(tem) <- c(colnames(tem)[1:4],"Loci","MinCorrSp","MinCorrCp","MinSp","MinCp","MinSpTrait","SpSigType","CpSigType","MinCpNum",paste0(c("Rank","Sp","Cp"),rep(1:cgwasenv$.pheNum,each=3)))
  write.csv(tem,file.path(cgwasenv$.cgwas_result, "AllHits.csv"),row.names=F,quote=F)

  topsnpi <- c()
  Sp_or_Cp <- c()
  for(i in unique(locic)){
    gt <- intersect(which(locic==i),as.numeric(which(unlist(tem[,6])==min(as.numeric(tem[which(locic==i),6])))))
    ct <- intersect(which(locic==i),as.numeric(which(unlist(tem[,7])==min(as.numeric(tem[which(locic==i),7])))))
    temv <- union(ct,gt)
    indv <- order(temv)
    topsnpi <- c(topsnpi,temv[indv])
    gv <- rep(0,length(temv))
    for(j in 1:length(gt)){
      gv[which(temv==gt[j])] <- 1
    }
    cv <- rep(0,length(temv))
    for(j in 1:length(ct)){
      cv[which(temv==ct[j])] <- 2
    }
    mv <- c(gv+cv)[indv]
    mv[mv==1] <- "S"
    mv[mv==2] <- "C"
    mv[mv==3] <- "SC"
    Sp_or_Cp <- c(Sp_or_Cp,mv)
  }
  print(paste0("Pick Top Hits Completed"))
  write.csv(cbind(tem[topsnpi,1:5],Sp_or_Cp,tem[topsnpi,6:ncol(tem)]),file.path(cgwasenv$.cgwas_result, "TopSNP.csv"),row.names=F,quote=F)

  step7Time2 <- Sys.time()
  print(paste0("step7 use time:"))
  print(step7Time2 - step7Time1)
  print(paste0("Step7 Completed"))
}
