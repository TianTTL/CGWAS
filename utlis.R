corterm3all_parallel <- function(i){
  dm <- as.matrix(data.table::fread(file.path(rowdataPath, paste0("All.s.", i)),header=F))
  pm <- as.matrix(data.table::fread(file.path(rowdataPath, paste0("All.sp.", i)),header=F))

  cp <- matrix(0, nrow = nrow(dm), ncol = ncol(dm))
  for (j in 1:nrow(dm)){
    b <- dm[j,]
    p <- pm[j,]
    temv <- which(b<0)
    ocorms <- ocorm
    ocorms3 <- ocorm3
    ocorms[temv,] <- -ocorm[temv,]
    ocorms[,temv] <- -ocorm[,temv]
    ocorms3[temv,] <- -ocorm3[temv,]
    ocorms3[,temv] <- -ocorm3[,temv]
    corm <- 3.2630398097*ocorms + 0.7095678755*ocorm2 + 0.0268257772*ocorms3 + 0.0005732151*ocorm4
    stepw <- order(p)
    p <- p[stepw]
    corm <- corm[stepw,stepw]
    stat <- -2*log(p)
    curcov <- rep(0,pheNum)
    for(k in 2:pheNum){
      stat[k] <- stat[k]+stat[k-1]
      curcov[k] <- curcov[k-1]+sum(corm[1:(k-1),k])
    }
    E <- 2*(1:pheNum)
    V <- 2*curcov + 2*E
    cp[j,] <- pchisq(2*E*stat/V,2*E^2/V,lower.tail=F)
  }
  return(cp)
}


mkdf_corterm3all_parallel <- function(i){
  b <- mvrnorm(1,mu=rep(0,pheNum),Sigma=ocorm)
  b <- b*lambda
  p <- pnorm(-abs(b))*2
  temv <- which(b<0)
  ocorm[temv,] <- -ocorm[temv,]
  ocorm[,temv] <- -ocorm[,temv]
  ocorm3[temv,] <- -ocorm3[temv,]
  ocorm3[,temv] <- -ocorm3[,temv]
  corm <- 3.2630398097*ocorm + 0.7095678755*ocorm2 + 0.0268257772*ocorm3 + 0.0005732151*ocorm4
  stepw <- order(p)
  p <- p[stepw]
  corm <- corm[stepw,stepw]
  stat <- -2*log(p)
  curcov <- rep(0,pheNum)
  for(j in 2:pheNum){
    stat[j] <- stat[j]+stat[j-1]
    curcov[j] <- curcov[j-1]+sum(corm[1:(j-1),j])
  }
  E <- 2*(1:pheNum)
  V <- 2*curcov + 2*E
  cp <- pchisq(2*E*stat/V,2*E^2/V,lower.tail=F)
  return(cp)
}

inflationplot <- function(infcoef,Esp,title){
  e <- -log10(ppoints(Esp))[cgwasenv$.qv]
  plot(c(0,max(e)), c(min(infcoef),max(infcoef)), type="n" ,xlab = "Expected -LogP", ylab = "Inflation Factor", main = title)
  points(e, infcoef, pch=20)
  abline(c(1,0),col="red")
  abline(c(2,0),col="red")
  abline(c(5,0),col="red")
  abline(c(10,0),col="red")
  abline(c(20,0),col="red")
  abline(c(50,0),col="red")
  abline(c(100,0),col="red")
  ifpos <- order(abs(e+log10(0.5)))[1:2]
  points(mean(e[ifpos]), mean(infcoef[ifpos]),cex=2,pch=21,col="black",bg="white")
}

manhattanPlot <- function(indscp,spt.col=c('gray10','gray50'),cpt.col=c('gray10','gray50'),spt.bg=c('gray10','gray50'),cpt.bg=c('gray10','gray50'),spt.cex=0.7,cpt.cex=0.7,pch=21,cex.axis=1.1,suggestiveline.col='blue',suggestiveline.lwd=1.5,suggestiveline.lty=1,genomewideline.col='red',genomewideline.lwd=1.5,genomewideline.lty=1){
  d=indscp
  colnames(d)=c("CHR","BP","SNP","SP","CP")
  d$pos=NA

  ytop=max(ceiling(max(d$CP)),ceiling(max(d$SP)))
  if(ytop>10){
    d$CP[d$CP>10]=d$CP[d$CP>10]/5+8
    d$SP[d$SP>10]=d$SP[d$SP>10]/5+8
    if(ytop>50){
      d$CP[d$CP>18]=d$CP[d$CP>18]/5+14.4
      d$SP[d$SP>18]=d$SP[d$SP>18]/5+14.4
      if(ytop>200){
        d$CP[d$CP>24]=24.4
        d$SP[d$SP>24]=24.4
      }
    }
  }
  ymax=ceiling(max(d$CP))
  ymin=-ceiling(max(d$SP))

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
      ticks[i+1]=max(d[d$index==i,]$pos)
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

  l=c(-seq(2,24,2)[floor(-ymin/2):1],0,seq(2,24,2)[1:floor(ymax/2)])
  ll=c(c(2,4,6,8,10,20,30,40,50,100,150,200)[floor(-ymin/2):1],0,c(2,4,6,8,10,20,30,40,50,100,150,200)[1:floor(ymax/2)])
  axis(2,las=1,at=l,labels=as.character(ll),lwd=0,lwd.ticks=1,cex.axis=cex.axis)

  d$SP[which(d$SP==0)] <- NA
  d$CP[which(d$CP==0)] <- NA
  pt.col = rep(cpt.col,max(d$CHR))[1:max(d$CHR)]
  pt.bg = rep(cpt.bg,max(d$CHR))[1:max(d$CHR)]

  icol=1
  for (i in unique(d$CHR)) {
    with(d[d$CHR==i, ],points(pos, CP, col=pt.col[icol],bg=pt.bg[icol],cex=cpt.cex,pch=pch))
    icol=icol+1
  }

  pt.col = rep(spt.col,max(d$CHR))[1:max(d$CHR)]
  pt.bg = rep(spt.bg,max(d$CHR))[1:max(d$CHR)]

  icol=1
  for (i in unique(d$CHR)) {
    with(d[d$CHR==i, ],points(pos, -SP, col=pt.col[icol],bg=pt.bg[icol],cex=spt.cex,pch=pch))
    icol=icol+1
  }

  text(newtick,-0.25,lowerlabs,cex=1.1)
  text(newtick,0.25,upperlabs,cex=1.1)

  abline(h=-log10(cgwasenv$.St), col=suggestiveline.col[1],lwd=suggestiveline.lwd,lty=suggestiveline.lty)
  abline(h=log10(cgwasenv$.St), col=suggestiveline.col[1],lwd=suggestiveline.lwd,lty=suggestiveline.lty)
  abline(h=-log10(cgwasenv$.Gt), col=genomewideline.col[1],lwd=genomewideline.lwd,lty=genomewideline.lty)
  abline(h=log10(cgwasenv$.Gt), col=genomewideline.col[1],lwd=genomewideline.lwd,lty=genomewideline.lty)

  box()
}

preCorrection <- function(){
  c1 <- ceiling(cgwasenv$.Esnpnum*0.001)
  cc1 <- ceiling(cgwasenv$.Esnpnum*0.00001)
  c2 <- ceiling(cgwasenv$.Esnpnum*0.01)
  cc2 <- ceiling(cgwasenv$.Esnpnum*0.0001)
  c3 <- ceiling(cgwasenv$.Esnpnum*0.1)
  if(cgwasenv$.Esnpnum>=1e6){
    cgwasenv$.qv <- c(1:(c1-1),seq(c1,(c2-1),cc1),seq(c2,(c3-1),cc2),seq(c3,cgwasenv$.Esnpnum-100,c1),(cgwasenv$.Esnpnum-99):cgwasenv$.Esnpnum)
  } else if(cgwasenv$.Esnpnum>=1e5){
    cgwasenv$.qv <- c(1:(c2-1),seq(c2,(c3-1),cc2),seq(c3,cgwasenv$.Esnpnum-100,c1),(cgwasenv$.Esnpnum-99):cgwasenv$.Esnpnum)
  } else if(cgwasenv$.Esnpnum>=1e4){
    cgwasenv$.qv <- c(1:(c3-1),seq(c3,cgwasenv$.Esnpnum-100,c1),(cgwasenv$.Esnpnum-99):cgwasenv$.Esnpnum)
  } else{
    cgwasenv$.qv <- 1:cgwasenv$.Esnpnum
  }
}

qqPlot <- function(datasm,datacm,minm,e,sgif,cgif,twoc=c('gray10','gray50')){
  myc <- c(colorRampPalette(c("navy","blue","cyan","chartreuse","yellow","orange","Red","darkred"))(1000))

  axis2 <- function(ym){
    l=c(0,seq(2,24,2)[1:floor(ym/2)])
    ll=c(0,c(2,4,6,8,10,20,30,40,50,100,150,200)[1:floor(ym/2)])
    axis(2,las=1,at=l,labels=as.character(ll),lwd=0,lwd.ticks=1,cex.axis=0.95)
  }

  ytop <- max(datasm[1,])
  if(ytop>10){
    datasm[datasm>10]=datasm[datasm>10]/5+8
    if(ytop>50){
      datasm[datasm>18]=datasm[datasm>18]/5+14.4
      if(ytop>200){
        datasm[datasm>24]=24.4
      }
    }
  }

  xmax <- max(e)*1.03
  xmin <- max(e)*-0.03
  ymax <- max(datasm[1,])*1.03
  ymin <- -ymax*0.03
  plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),col=F,las=1,xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=0.95)
  axis(side=1,labels=seq(0,xmax,1),at=seq(0,xmax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  axis2(ymax)

  for(i in 1:pheNum){
    points(e,datasm[,i],pch=20,cex=0.9,col=myc[round(i*1000/pheNum+1)])
    print(paste0(i," GWAS QQ Plots Completed"))
  }
  abline(0,1,col="black",lwd=1.5,lty=1)
  box()
  print(paste0("Overlap GWAS QQ Plots Completed"))


  xmax <- 3*1.03
  xmin <- 3*-0.03
  ymax <- 3*1.03
  ymin <- -ymax*0.03
  plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),col=F,las=1,xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=0.95)
  axis(side=1,labels=seq(0,xmax,1),at=seq(0,xmax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  axis(side=2,las=1,labels=seq(0,ymax,1),at=seq(0,ymax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  for(i in 1:pheNum){
    cpind <- datasm[,i]<ymax
    points(e[cpind],datasm[cpind,i],pch=20,cex=0.9,col=myc[round(i*1000/(pheNum+1))])
    print(paste0(i," Local GWAS QQ Plots Completed"))
  }
  abline(0,1,col="black",lwd=1.5,lty=1)
  box()
  print(paste0("Local GWAS QQ Plots Completed"))


  ytop <- max(minm)
  if(ytop>10){
    minm[minm>10]=minm[minm>10]/5+8
    if(ytop>50){
      minm[minm>18]=minm[minm>18]/5+14.4
      if(ytop>200){
        minm[minm>24]=24.4
      }
    }
  }

  xmax <- max(e)*1.03
  xmin <- max(e)*-0.03
  ymax <- max(minm)*1.03
  ymin <- -ymax*0.03
  plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),col=F,las=1,xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=0.95)
  axis(side=1,labels=seq(0,xmax,1),at=seq(0,xmax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  axis2(ymax)

  points(e,minm[,1],pch=20,cex=0.9,col=twoc[1])
  points(e,minm[,2],pch=20,cex=0.9,col=twoc[2])
  abline(0,1,col="black",lwd=1.5,lty=1)
  box()
  print(paste0("Minp QQ Plots Completed"))


  ytop <- max(datacm[1,])
  if(ytop>10){
    datacm[datacm>10]=datacm[datacm>10]/5+8
    if(ytop>50){
      datacm[datacm>18]=datacm[datacm>18]/5+14.4
      if(ytop>200){
        datacm[datacm>24]=24.4
      }
    }
  }

  xmax <- max(e)*1.03
  xmin <- max(e)*-0.03
  ymax <- max(datacm[1,])*1.03
  ymin <- -ymax*0.03
  plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),col=F,las=1,xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=0.95)
  axis(side=1,labels=seq(0,xmax,1),at=seq(0,xmax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  axis2(ymax)

  for(i in 1:pheNum){
    points(e,datacm[,i],pch=20,cex=0.9,col=myc[round(i*1000/(pheNum+1))])
    print(paste0(i," CGWAS QQ Plots Completed"))
  }
  abline(0,1,col="black",lwd=1.5,lty=1)
  box()
  print(paste0("Overlap CGWAS QQ Plots Completed"))


  xmax <- 3*1.03
  xmin <- 3*-0.03
  ymax <- 3*1.03
  ymin <- -ymax*0.03
  plot(0,xlab=expression(Expected~~-log[10](italic(p))),ylab=expression(Observed~~-log[10](italic(p))),col=F,las=1,xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=0.95)
  axis(side=1,labels=seq(0,xmax,1),at=seq(0,xmax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  axis(side=2,las=1,labels=seq(0,ymax,1),at=seq(0,ymax,1),cex.axis=0.95,lwd=0,lwd.ticks=1)
  for(i in 1:pheNum){
    cpind <- datacm[,i]<ymax
    points(e[cpind],datacm[cpind,i],pch=20,cex=0.9,col=myc[round(i*1000/(pheNum+1))])
    print(paste0(i," Local CGWAS QQ Plots Completed"))
  }
  abline(0,1,col="black",lwd=1.5,lty=1)
  box()
  print(paste0("Local CGWAS QQ Plots Completed"))


  xmax <- pheNum*1.03
  xmin <- pheNum*-0.03+1
  ymax <- max(max(sgif),max(cgif),1)+(max(max(sgif),max(cgif),1)-min(min(sgif),min(cgif),1))*0.03
  ymin <- min(min(sgif),min(cgif),1)-(max(max(sgif),max(cgif),1)-min(min(sgif),min(cgif),1))*0.03
  plot(0,xlab="Combined Trait Number",ylab="Inflation Factor",col=F,las=1,xaxt='n',yaxt='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty='n',xaxs='i',yaxs='i',cex.axis=0.95)
  axis(side=1,cex.axis=0.95,lwd=0,lwd.ticks=1)
  axis(side=2,las=1,cex.axis=0.95,lwd=0,lwd.ticks=1)

  for(ii in seq(0.05,0.45,1/150)){
    rect(xleft=xmin,ybottom=quantile(sgif,1-ii),xright=xmax,ytop=quantile(sgif,ii),col=paste0("grey",97.5-150*ii),border=paste0("grey",97.5-150*ii))
  }
  abline(c(quantile(sgif,1),0),lty=3,lwd=1.5)
  abline(c(quantile(sgif,0.95),0),lty=2,lwd=1.5)
  abline(c(quantile(sgif,0.5),0),lwd=1.5)
  abline(c(quantile(sgif,0.05),0),lty=2,lwd=1.5)
  abline(c(quantile(sgif,0),0),lty=3,lwd=1.5)
  for(ii in 1:(pheNum-1)){
    lines(c(ii,ii+1),c(cgif[ii],cgif[ii+1]),col=myc[round(ii*1000/pheNum)],lwd=1.5)
  }
  points(1:pheNum,cgif,col=myc[round((1:pheNum)*1000/(pheNum+1))],bg=myc[round((1:pheNum)*1000/(pheNum+1))],pch=23,cex=1.5)
  box()

  print(paste0("Lambda Plots Completed"))
}

simulQQplot <- function(pm,title,simulNum,Esp){
  e <- -log10(ppoints(Esp))[cgwasenv$.qv]
  hpv <- apply(pm,2,function(a){return(quantile(a,0.05))})
  eh <- -log10(hpv)
  em <- -log10(apply(pm,2,function(a){return(quantile(a,0.5))}))
  el <- -log10(apply(pm,2,function(a){return(quantile(a,0.95))}))
  ttp <- -log10(hpv[1])
  plot(c(0,max(e)), c(0,ttp), type="n" ,xlab = "Expected -LogP", ylab = "Observed -LogP", main = title)
  polygon(c(e,e[length(e):1]),c(el,eh[length(e):1]),border="black",col="grey80")
  points(e, em, pch=20)
  text(0.5,max(ttp),expression(italic(P)[Sim]),adj=1)
  text(0.5,max(ttp)*0.95,expression(italic(P)[Set]),adj=1)
  text(0.5,max(ttp)*0.9,expression(italic(FDR)[Set]),adj=1)
  text(0.5,max(ttp),paste0(" = ",signif(hpv[1],3)),adj=0)
  text(0.5,max(ttp)*0.95,paste0(" = ",signif(0.05/Esp,3)),adj=0)
  text(0.5,max(ttp)*0.9,paste0(" = ",signif(mean(order(abs(pm[order(pm[,1]),1]-0.05/Esp),decreasing=F)[1:2])/simulNum,3)),adj=0)
  abline(c(0,1),col="red")
}

simultime <- function(pm,simulNum,Esp,SimuReg,title){
  tseq <- seq(40,round(simulNum*SimuReg),30)
  pv95 <- matrix(0,1000,length(tseq))
  for(j in 1:1000){
    temord <- sample(1:simulNum,round(simulNum*SimuReg))
    for(i in 1:length(tseq)){
      pv95[j,i] <- mean(order(abs(pm[temord[1:tseq[i]][order(pm[temord[1:tseq[i]],1])],1]-0.05/Esp),decreasing=F)[1:2])/tseq[i]
    }
  }
  eh <- apply(pv95,2,function(a){return(quantile(a,0.05))})
  em <- apply(pv95,2,function(a){return(quantile(a,0.5))})
  el <- apply(pv95,2,function(a){return(quantile(a,0.95))})
  plot(c(40,round(simulNum*SimuReg)), c(min(min(c(eh,el)),0.05),max(max(c(eh,el)),0.05)), type="n" ,xlab = "Simu Time", ylab = "FDR", main = title)
  polygon(c(tseq,tseq[length(tseq):1]),c(el,eh[length(em):1]),border="black",col="grey80")
  lines(tseq,em,lwd=2)
  abline(c(0.05,0),col="red")
  abline(c(el[length(em)],0),col="red",lty=2)
  abline(c(eh[length(em)],0),col="red",lty=2)
}
