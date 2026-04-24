library(rethinking)
source("00_functions.R")
source("01_process_data.R")
source("02_define_traditions.R")
traditions<-keep_first_pioneer(traditions) #Keep just the most important pioneer

#Load the posterior
load("post_reduced_step_back7.Rdata")

precis(post.r)

#nicer visualization of the posterior
basecols<-c("#FFE000","#FF2255","#22FF55","#0060FF")

postPlot<-function(post,K,sc=0.03,cols=rep("#FF2255",length(post)),gcol="#A0A0A0",
                   labs=names(post),
                   ofs=rep(0,length(post))){
  
  h<-length(post)+3
  y<-c(1:length(labs))+cumsum(ofs)
  
  tiff("Posterior.tiff",width=7.8,height=max(y)+1,units="cm",res=300,compression="lzw")
  #layout(matrix(c(1,2),ncol=1),heights=c(h,h2))
  par(mar=c(8,5,1,1),mgp=c(2,0.65,0))
  plot(NULL,xlim=c(0,1),ylim=c(max(y),1)+c(1,-1)*0.2,
       xlab="Parameter value", ylab="", yaxt="n",xaxs="i")
  
  rect(-1000,-1000,1000,1000,col="#FBFBFA",border=NA)
  axis(2,at=y,labels=labs,las=1)
  
  #White background  
  for(i in 1:length(post)) {
    abline(h=y[i],col=gcol,lty=1,lwd=1)
    toplot<- post[[i]]
    dens <- density(toplot)
    take<-which(dens$x >= min(toplot)  & dens$x <= max(toplot))
    polygon(c(dens$x[take],rev(dens$x[take])),y[i]+c(dens$y[take],rev(-dens$y)[take])*sc,col="white", border="white")
  }
  
  #Grid
  abline(v=seq(0.1,0.9,0.1),col=gcol,lty=2,lwd=0.6)
  abline(v=seq(0.2,0.8,0.2),col=gcol,lty=1,lwd=0.8)
  
  #color overlay
  for(i in 1:length(post)) {
    abline(h=y[i],col=gcol,lty=1,lwd=1)
    toplot<- post[[i]]
    dens <- density(toplot)
    take<-which(dens$x >= min(toplot)  & dens$x <= max(toplot))
    polygon(c(dens$x[take],rev(dens$x[take])),y[i]+c(dens$y[take],rev(-dens$y)[take])*sc,col=paste0(cols[i],"AA"), border="#808080A0",xpd=NA)
    lines(HPDI(toplot,prob=0.95),rep(y[i],2),lwd=2,col=1)
    points(mean(toplot),y[i],pch=21,bg="white",lwd=1.5,cex=1,col=1)
  }
  
  #title("Posterior distribution",adj=0,cex=0.9)
  
  legend("bottomleft",legend=c("attraction","error term","exploration term 1 (own style)","exploration term 2 (within career)"),
         fill=paste0(basecols,"AA"),bty="n",
         inset=c(-0.55, -.39),xpd=NA,ncol=1)
  
  box()
  dev.off()
}


postPlot(post.r,cols=rep(basecols,c(1,4,4,4)),
         ofs=c(0,0.3,0,0,0,0.3,0,0,0,0.3,0,0,0),
         labs=c(expression(bar(alpha)),rep(c(expression(eta["baseline"]),
                          expression(nu["dimension"]),
                          expression(nu["pioneer (avg.)"]),
                          expression(nu["dim., pioneer"])),times=3)))
