library(rethinking)
source("00_functions.R")
source("01_process_data.R")
source("02_define_traditions.R")
traditions<-keep_first_pioneer(traditions) #Keep just the most important pioneer

#Reconstruct the data and Load the posterior
data <- prepare_stan_data(D, d$namlong, traditions)
load("post_reduced_step_back7.Rdata")

#You must have object data from the steps that precede model fitting boiled dows
str(post.r)
post <- post.r

w<-10

tiff("Prediction.tiff",width=w*1.3,height=w*1.3,units="cm",res=300,compression="lzw")
layout(matrix(c(1,3,2,3),ncol=2),heights=c(2.8,1))

cols<-c("#FF2255","#22FF55","#0060FF")
pcol<-"#000000"
subbg<-"#FFFFFFB5" #color to tone under compatibility intervals

gcol<-"#C0C0C0"

#Which pioneers have largest and smallest average dim variance?
sd_rm<-data$sd_rm_avg
sds_rm<-data$sds_rm #Here each tradition is a row

which.max(sd_rm)
which.min(sd_rm)

#Take the minimum variance
rID<-which.min(sd_rm)
names(traditions)[rID]

#extract sds per dimension
sds<-data$sds

#How dimensions will be ordered
ord<-order(sds_rm[rID,],decreasing=T)

#Start adding posterior predictions
sdd<- sds[ord]
sdr<- rep(sd_rm[rID],300)
sdp<- sds_rm[rID,][ord]

#create new contrafactual values for sdp to fill the whole image with predictions
sdp_new<-seq(0,1.5,by=0.01)
sdr_one<-sd_rm[rID] #Sometimes we will use role model average sd as a single number and not a repeated vector

#how many samples and dimensions
s<-dim(post[[1]])[1]
dims<-300

#Error variance per epigone
lsig1<-with(post,sapply(1:s,function(s){sqrt(lsb[s,]^2 + lnur[s,]^2 * sdr_one^2 + lnup[s,]^2 * sdp_new^2)}))
lmu1<-apply(lsig1,1,mean)
lPI1<-apply(lsig1,1,PI,prob=0.9)

#Error variance per epigone variance total
lsig<-with(post,sapply(1:s,function(s){sqrt(lsb[s,]^2 + lnud[s,]^2 * sdd^2 + lnur[s,]^2 * sdr^2 + lnup[s,]^2 * sdp^2)}))
lmu<-apply(lsig,1,mean)
lPI<-apply(lsig,1,PI,prob=0.9)

#Exploratory variance per epigone
dsig1<-with(post,sapply(1:s,function(s){sqrt(dsb[s,]^2 + dnur[s,]^2 * sdr_one^2 + dnup[s,]^2 * sdp_new^2)}))
dmu1<-apply(dsig1,1,mean)
dPI1<-apply(dsig1,1,PI,prob=0.9)

#Exploratory variance per epigone total
dsig<-with(post,sapply(1:s,function(s){sqrt(dsb[s,]^2 + dnud[s,]^2 * sdd^2 + dnur[s,]^2 * sdr^2 + dnup[s,]^2 * sdp^2)}))
dmu<-apply(dsig,1,mean)
dPI<-apply(dsig,1,PI,prob=0.9)

#Exploratory variance, no general dimension sd
sig1<-with(post,sapply(1:s,function(s){sqrt(sb[s,]^2 + nur[s,]^2 * sdr_one^2 + nup[s,]^2 * sdp_new^2)}))
mu1<-apply(sig1,1,mean)
PI1<-apply(sig1,1,PI,prob=0.9)

#Exploratory variance total
sig<-with(post,sapply(1:s,function(s){sqrt(sb[s,]^2 + nud[s,]^2 * sdd^2 + nur[s,]^2 * sdr^2 + nup[s,]^2 * sdp^2)}))
mu<-apply(sig,1,mean)
PI<-apply(sig,1,PI,prob=0.9)

#Linear sd predition plot
par(mar=c(3.5,3.5,2,1),mgp=c(2,0.65,0))
plot(NULL,xlim=c(0.4,1.2),ylim=c(0,2),xlab="Pioneer SD",ylab="SD",xaxs="i",yaxs="i", asp=1)
title("Crime fiction",adj=0)

hwr<-0.005 # half width of rectangle
rect(sdp-hwr,rep(0,dims),sdp+hwr,sdd,col=paste0("#008080",60),border=NA)

points(sdp,sdp,col="white",pch=19)

shade(lPI1,sdp_new,col=subbg)
shade(dPI1,sdp_new,col=subbg)
shade(PI1,sdp_new,col=subbg)

abline(h=seq(0,2,0.1),col=gcol,lty=2)
abline(v=seq(0,2,0.1),col=gcol,lty=2)
abline(0,1,lty=2,col=1)

lines(sdp,lmu,col=cols[1])
lines(sdp,dmu,col=cols[2])
lines(sdp,mu,col=cols[3])

abline(h=sdr_one,col="#FF00FF") #Parental average

points(sdp,sdp,col=paste0(pcol,80),pch=19)

shade(lPI1,sdp_new,col=paste0(cols[1],20))
lines(sdp_new,lmu1,col=cols[1])

shade(dPI1,sdp_new,col=paste0(cols[2],20))
lines(sdp_new,dmu1,col=cols[2])

shade(PI1,sdp_new,col=paste0(cols[3],20))
lines(sdp_new,mu1,col=cols[3])


#Take the maximum variance group
rID<-which.max(sd_rm)
names(traditions)[rID]

#extract sds per dimension
sds<-data$sds

#How dimensions will be ordered
ord<-order(sds_rm[rID,],decreasing=T)

#Start adding posterior predictions
sdd<- sds[ord]
sdr<- rep(sd_rm[rID],300)
sdp<- sds_rm[rID,][ord]

#create new contrafactual values for sdp to fill the whole image with predictions
sdp_new<-seq(0,1.5,by=0.01)
sdr_one<-sd_rm[rID] #Sometimes we will use role model average sd as a single number and not a repeated vector

#how many samples and dimensions
s<-dim(post[[1]])[1]
dims<-300

#Error variance per epigone
lsig1<-with(post,sapply(1:s,function(s){sqrt(lsb[s,]^2 + lnur[s,]^2 * sdr_one^2 + lnup[s,]^2 * sdp_new^2)}))
lmu1<-apply(lsig1,1,mean)
lPI1<-apply(lsig1,1,PI,prob=0.9)

#Error variance per epigone variance total
lsig<-with(post,sapply(1:s,function(s){sqrt(lsb[s,]^2 + lnud[s,]^2 * sdd^2 + lnur[s,]^2 * sdr^2 + lnup[s,]^2 * sdp^2)}))
lmu<-apply(lsig,1,mean)
lPI<-apply(lsig,1,PI,prob=0.9)

#Exploratory variance per epigone
dsig1<-with(post,sapply(1:s,function(s){sqrt(dsb[s,]^2 + dnur[s,]^2 * sdr_one^2 + dnup[s,]^2 * sdp_new^2)}))
dmu1<-apply(dsig1,1,mean)
dPI1<-apply(dsig1,1,PI,prob=0.9)

#Exploratory variance per epigone total
dsig<-with(post,sapply(1:s,function(s){sqrt(dsb[s,]^2 + dnud[s,]^2 * sdd^2 + dnur[s,]^2 * sdr^2 + dnup[s,]^2 * sdp^2)}))
dmu<-apply(dsig,1,mean)
dPI<-apply(dsig,1,PI,prob=0.9)

#Exploratory variance, no general dimension sd
sig1<-with(post,sapply(1:s,function(s){sqrt(sb[s,]^2 + nur[s,]^2 * sdr_one^2 + nup[s,]^2 * sdp_new^2)}))
mu1<-apply(sig1,1,mean)
PI1<-apply(sig1,1,PI,prob=0.9)

#Exploratory variance total
sig<-with(post,sapply(1:s,function(s){sqrt(sb[s,]^2 + nud[s,]^2 * sdd^2 + nur[s,]^2 * sdr^2 + nup[s,]^2 * sdp^2)}))
mu<-apply(sig,1,mean)
PI<-apply(sig,1,PI,prob=0.9)


#Second panel
plot(NULL,xlim=c(0.4,1.2),ylim=c(0,2),xlab="Pioneer SD",ylab="SD",xaxs="i",yaxs="i", asp=1)
title("Sci-fi",adj=0)

hwr<-0.005 # half width of rectangle
rect(sdp-hwr,rep(0,dims),sdp+hwr,sdd,col=paste0("#008080",60),border=NA)

points(sdp,sdp,col="white",pch=19)

shade(lPI1,sdp_new,col=subbg)
shade(dPI1,sdp_new,col=subbg)
shade(PI1,sdp_new,col=subbg)

abline(h=seq(0,2,0.1),col=gcol,lty=2)
abline(v=seq(0,2,0.1),col=gcol,lty=2)
abline(0,1,lty=2,col=1)

lines(sdp,lmu,col=cols[1])
lines(sdp,dmu,col=cols[2])
lines(sdp,mu,col=cols[3])

abline(h=sdr_one,col="#FF00FF") #Parental average

points(sdp,sdp,col=paste0(pcol,80),pch=19)

shade(lPI1,sdp_new,col=paste0(cols[1],20))
lines(sdp_new,lmu1,col=cols[1])

shade(dPI1,sdp_new,col=paste0(cols[2],20))
lines(sdp_new,dmu1,col=cols[2])

shade(PI1,sdp_new,col=paste0(cols[3],20))
lines(sdp_new,mu1,col=cols[3])


#Legend
par(mar=c(0,0,0,0))
plot(NULL,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
legend=c("Mean pioneer SD",
         "Referential line, x = y",
         "Pioneer SD per dimension",
         "(bars) Dimension SD",
         "Predicted epigone SD, without Dimension SD ± 90% CI",
         "Predicted epigone SD, total")

y<-seq(0.9,0.1,length.out=length(legend))
x1<-0.05
of<-0.05
of2<-0.07

x2<-x1+of
for(i in 1:length(legend)){
  text(x1+of2,y[i],legend[i],col="black",pos=4)
}
lines(c(x1,x2),rep(y[1],2),col="#FF00FF",lty=1)
lines(c(x1,x2),rep(y[2],2),col=1,lty=2)
points(mean(c(x1,x2)),y[3],col=paste0(pcol,80),pch=19)
bars<-c(1,1.5,0.7,1.8,0.5,0.3,0.8,1.5,1,0.8)
wr<-0.005
for(i in 1:length(bars)){
  rect(x1+(i-1)*wr,y[4]-0.03,x1+i*wr,y[4]-0.03+bars[i]*0.05,col=paste0("#008080",60),border=NA)
}
rcol<-"#808080"
rect(x1,y[5]-0.03,x2,y[5]+0.03,col=paste0(rcol,40),border=NA)
lines(c(x1,x2),rep(y[5],2),col=rcol,lwd=1)
lines(x1+((1:(length(bars)*2))-1)*wr/2,y[6]-0.03+sample(rep(bars,2))*0.05,col=1,lwd=1)

#Second columns of legend
x1<-0.5
x2<-x1+of

legend=c("error term",
         "exploration term 1 (own style)",
         "exploration term 2 (within career)")

for(i in 1:length(legend)){
  text(x1+of2,y[i],legend[i],col="black",pos=4)
}
for(i in 1:3){
  rect(x1,y[i]-0.03,x2,y[i]+0.03,col=paste0(cols[i],40),border=NA)
  lines(c(x1,x2),rep(y[i],2),col=cols[i],lty=1)  
}
dev.off()





