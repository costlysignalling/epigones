library(rethinking)
source("00_functions.R")
source("01_process_data.R")
source("02_define_traditions.R")
traditions<-keep_first_pioneer(traditions)

#Run something that established traditions first
traditions_mu <- lapply(traditions, function(tr) {
  # pooled SD for all pioneers
  D_p       <- D[ d$namlong %in% tr$pioneers, , drop = FALSE ]
  mu_pioneers <- apply(D_p, 2, mean)
  # individual SD for each epigone
  mu_epigones <- setNames(
    lapply(tr$epigones, function(ep) {
      D_e      <- D[ d$namlong == ep, , drop = FALSE ]
      apply(D_e, 2, mean)
    }),
    tr$epigones
  )
  list(
    pioneers = mu_pioneers,
    epigones = mu_epigones
  )
})

#SD of epigone mu
traditions_mu_sd <- lapply(traditions_mu, function(tr) {
  mat <- do.call(cbind, tr$epigones)    # 300 rows (dims) × K cols (epigones)
  apply(mat, 1, sd)
})

# compute SDs
traditions_sd <- lapply(traditions, function(tr) {
  # pooled SD for all pioneers
  D_p       <- D[ d$namlong %in% tr$pioneers, , drop = FALSE ]
  sd_pioneers <- apply(D_p, 2, sd2)
  # individual SD for each epigone
  sd_epigones <- setNames(
    lapply(tr$epigones, function(ep) {
      D_e      <- D[ d$namlong == ep, , drop = FALSE ]
      apply(D_e, 2, sd2)
    }),
    tr$epigones
  )
  list(
    pioneers = sd_pioneers,
    epigones = sd_epigones
  )
})

str(traditions_mu_sd)
str(traditions_sd)

#Average SD per pioneer
sd_rm_avg<-sapply(traditions_sd, function(tr) { sqrt(mean(tr$pioneers^2)) }) 

#convert to long format
# Between epigone mean sd
mu_mat      <- do.call(cbind, traditions_mu_sd)
p_mat       <- do.call(cbind, lapply(traditions_sd, `[[`, "pioneers"))

# number of traditions (T) and dimensions (D)
T           <- length(traditions_mu_sd)
dim           <- nrow(mu_mat)

# build the long data.frame
dlb  <- data.frame(
  pID           = rep(seq_len(T), each = dim),
  dim           = rep(seq_len(dim), times = T),
  mu_epigone_sd = as.vector(mu_mat),
  sdp    = as.vector(p_mat))

#add general dim variance
# we can also limit D to persons limited
Dsub<-D[d$namlong %in% unlist(traditions),]
Dpi<-D[d$namlong %in% unlist(lapply(traditions,function(tr){tr$pioneers})),] #Or only pioneers

sds  <- apply(D, 2, sd2)
#sds <- apply(Dsub, 2, sd2)
#sds  <- apply(Dpi, 2, sd2)

dlb$sdr <- sd_rm_avg[dlb$pID]
dlb$sdd <- sds[dlb$dim]

head(dlb, 6) # quick check

# Per epigone sd
# count total number of epigones
total_epig <- sum(sapply(traditions_sd, function(x) length(x$epigones)))

# prepare an empty list to collect each epigone’s block
res_list <- vector("list", total_epig)

ep_counter <- 1
for(p in names(traditions_sd)) {
  # integer code for this tradition
  pID_val <- which(names(traditions_sd) == p)
  
  # pooled pioneers SD (named length-300 vector)
  sdp <- traditions_sd[[p]]$pioneers
  dims       <- names(sdp)    # "D1","D2",…,"D300"
  n_dims     <- length(dims)          # should be 300
  
  # loop through each epigone in this tradition
  for(e in names(traditions_sd[[p]]$epigones)) {
    ep_sd   <- traditions_sd[[p]]$epigones[[e]]
    
    # build a small data.frame for this epigone × all dims
    res_list[[ep_counter]] <- data.frame(
      eID         = rep(ep_counter, n_dims),          # global epigone index
      pID         = rep(pID_val,   n_dims),           # tradition index
      dim         = seq_len(n_dims),                  # 1..300
      epigone_sd  = as.numeric(ep_sd),                # that author’s SD₂
      sdp  = as.numeric(sdp),          # pooled pioneers SD₂
      stringsAsFactors = FALSE
    )
    ep_counter <- ep_counter + 1
  }
}

# bind all blocks together
dl <- do.call(rbind, res_list)
rownames(dl) <- NULL

# inspect
head(dl, 10)

#Add general dim variance
dl$sdd <- sds[dl$dim]
dl$sdr <- sd_rm_avg[dl$pID]

pcol<-"#106090"

#Plots and models
h<-10
tiff("Basic_sd_visualization.tiff", width=2*h, height=h, units="cm", res=300, compression="lzw")
par(mfrow=c(1,2), mar=c(3.5,3.5,1,1), mgp=c(2,0.65,0))
plot(dl$sdp,dl$epigone_sd,col=paste0(pcol,20),pch=19, xlim=c(0,1.8),ylim=c(0,1.8),
     xlab="Within-Pioneer SD", ylab="Within-Epigone SD",xaxs="i", yaxs="i")
mW<-lm(epigone_sd ~ sdp, data = dl)
newx<-data.frame(sdp = seq(0,2,l= 100))
pred <- predict(mW, newdata = newx, interval = "confidence")
polygon(c(newx$sdp,rev(newx$sdp)),
         c(pred[, "lwr"], rev(pred[, "upr"])),
         col = "#80808050", border = NA)
abline(mW,lwd=1.5)
text(rep(0.02,2),c(1.7,1.55),paste(c("a =","b ="),format(round(summary(mW)$coefficients[,1],2),nsmall=2)),pos=4)

#Epigone sd vs pioneer sd
plot(dlb$sdp,dlb$mu_epigone_sd,col=paste0(pcol,20),pch=19,xlim=c(0,1.8),ylim=c(0,1.8),
     xlab="Within-Pioneer SD", ylab="Between-Epigone SD",xaxs="i", yaxs="i")
mB<-lm(mu_epigone_sd ~ sdp, data = dlb)
pred <- predict(mB, newdata = newx, interval = "confidence")
polygon(c(newx$sdp,rev(newx$sdp)),
        c(pred[, "lwr"], rev(pred[, "upr"])),
        col = "#80808050", border = NA)
abline(mB,lwd=1.5)
text(rep(0.02,2),c(1.7,1.55),paste(c("a =","b ="),format(round(summary(mB)$coefficients[,1],2),nsmall=2)),pos=4)
dev.off()

