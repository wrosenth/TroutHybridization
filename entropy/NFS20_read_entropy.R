#Processing entropy output
#WRosenthal 3/10/21
library(RColorBrewer)
library(rhdf5)
library(abind)

meta <- read.csv("./NFS20_meta.csv",header=T)


#start with k=2
qrep1 <- h5read("./NFS20_0.5_0.03_filt_2rep1.hdf5","q")
qrep2 <- h5read("./NFS20_0.5_0.03_filt_2rep2.hdf5","q")
qrep3 <- h5read("./NFS20_0.5_0.03_filt_2rep3.hdf5","q")

Q_rep1 <- h5read("./NFS20_0.5_0.03_filt_2rep1.hdf5","Q")
Q_rep2 <- h5read("./NFS20_0.5_0.03_filt_2rep2.hdf5","Q")
Q_rep3 <- h5read("./NFS20_0.5_0.03_filt_2rep3.hdf5","Q")

all_q <- abind(qrep1,qrep2,qrep3,along=1)
all_Q <- abind(Q_rep1,Q_rep2,Q_rep3,along=1)

q <- t(apply(all_q,2:3,mean))
Q <- t(apply(all_Q,2:3,mean))

q.ci <- apply(all_q,2:3,quantile,probs=c(0.025,0.975))
Q.ci <- apply(all_Q,2:3,quantile,probs=c(0.025,0.975))

q.ci_Width <- q.ci[2,2,] - q.ci[1,2,]
Q.ci_width <- Q.ci[2,2,] - Q.ci[1,2,]

j <- 1
pdf(file=paste0("convergence_plots_Q",j,".pdf"),onefile = T)
par(mfrow=c(2,2))
for (i in 1:nrow(meta)) {
  plot(Q_rep1[,2,i],type="l",ylim = c(0,1),xlim=c(0,5900),ylab=meta$Gen_ID[i])
  lines(Q_rep2[,2,i],col="blue")
  lines(Q_rep3[,2,i],col="red")
}
dev.off()

pdf(file=paste0("convergence_plots_lilq",j,".pdf"),onefile = T)
par(mfrow=c(2,2))
for (i in 1:nrow(meta)) {
  plot(qrep1[,1,i],type="l",ylim = c(0,1),xlim=c(0,5900),ylab=meta$Gen_ID[i])
  lines(qrep2[,1,i],col="blue")
  lines(qrep3[,1,i],col="red")
}
dev.off()

#dev_rep1 <- h5read("./NFS20_0.5_0.03_filt_2rep1.hdf5","deviance")
#dev_rep2 <- h5read("./NFS20_0.5_0.03_filt_2rep2.hdf5","deviance")
#dev_rep3 <- h5read("./NFS20_0.5_0.03_filt_2rep3.hdf5","deviance")

meta$q <- q[,1]
meta$big_Q <- Q[,2]
meta$q_ci_width <- q.ci_Width
meta$Q_ci_width <- Q.ci_width
write.csv(meta,file="NFS20_results_meta.csv",row.names = F,quote=F)

#-----------------------------------------------------------------------
results_meta <- read.csv("./NFS20_results_meta.csv",header=T)

###categorize fish into hybrid classes using q and Q

hybrid_status <- seq(1:nrow(results_meta))
for (i in 1:length(hybrid_status)){
  if (results_meta$q[i] < 0.1){hybrid_status[i] <- "RBT"}
  if (results_meta$q[i] > 0.9){hybrid_status[i] <- "YSC"}
  if (results_meta$q[i] < 0.6 & results_meta$q[i] > 0.4){
    if (results_meta$big_Q[i] > 0.8){hybrid_status[i] <- "F1"}
    if (results_meta$big_Q[i] > 0.4 & results_meta$big_Q[i] < 0.6){hybrid_status[i] <- "F2"}
  }
  if (hybrid_status[i] == i){
    calc <- results_meta$big_Q[i] - 2*results_meta$q[i]
    if (calc > -0.1 & calc < 0.1){hybrid_status[i] <- "BC RBT"}
    if (hybrid_status[i] == i){
      calc <- results_meta$big_Q[i] + 2*results_meta$q[i]
      if (calc > 1.9 & calc < 2.1){hybrid_status[i] <- "BC YSC"}
    }
  }
  if (hybrid_status[i] == i){hybrid_status[i] <- "Other"}
}
hybrid_status <- factor(hybrid_status,c("RBT","BC RBT","F1","F2","Other","BC YSC","YSC"))

#Middle Creek 2020 adults
results_meta_middle2020 <- results_meta[which(results_meta$Gen_ID %in% results_meta$Gen_ID[grep("WCR2020",results_meta$Gen_ID)] & !is.na(results_meta$Sex)),]
hybrid_status_middle2020 <- hybrid_status[which(results_meta$Gen_ID %in% results_meta$Gen_ID[grep("WCR2020",results_meta$Gen_ID)] & !is.na(results_meta$Sex))]
{par(pty="s",mfrow=c(1,2),mar=c(2,3,3,1))
  plot(results_meta_middle2020$q,results_meta_middle2020$big_Q, xlim=c(0,1), ylim=c(0,1),xlab="Proportion of YSC ancestry (q)",ylab="Interspecific ancestry (Q)",type="n",axes=F,main="Middle Creek 2020 adults")
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  #arrows(results_meta_middle2020$q.ci.width[1,2,], Q[,2], q.ci[2,2,], Q[,2], length=0, col="gray45")
  #arrows(q[,2], Q.ci[1,2,], q[,2],  Q.ci[2,2,], length=0, col="gray45")
  abline(v=0.1)
  abline(v=0.9)
  pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
  points(results_meta_middle2020$q,results_meta_middle2020$big_Q,col=pal[as.factor(hybrid_status_middle2020)],pch=21,cex=1.2,bg=adjustcolor(pal[as.factor(hybrid_status_middle2020)],alpha.f = 0.6))
  plot(1:10,1:10, type="n", xlab="", ylab="", axes=F)
  legend("center",legend=unique(hybrid_status_middle2020),text.width=2.3,x.intersp = 0.3,col=pal[as.factor(unique(hybrid_status_middle2020))],pt.bg = adjustcolor(pal[as.factor(unique(hybrid_status_middle2020))],alpha.f = 0.8),pch=21,bty="n",ncol=2)
}

#Middle Creek 2020 juveniles
results_meta_middle2020_juv <- results_meta[which(results_meta$Gen_ID %in% results_meta$Gen_ID[grep("WCR2020",results_meta$Gen_ID)] & is.na(results_meta$Sex)),]
hybrid_status_middle2020_juv <- hybrid_status[which(results_meta$Gen_ID %in% results_meta$Gen_ID[grep("WCR2020",results_meta$Gen_ID)] & is.na(results_meta$Sex))]
{par(pty="s",mfrow=c(1,2),mar=c(2,3,3,1))
  plot(results_meta_middle2020_juv$q,results_meta_middle2020_juv$big_Q, xlim=c(0,1), ylim=c(0,1),xlab="Proportion of YSC ancestry (q)",ylab="Interspecific ancestry (Q)",type="n",axes=F, main="Middle Creek 2020 juveniles")
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  #arrows(results_meta_middle2020$q.ci.width[1,2,], Q[,2], q.ci[2,2,], Q[,2], length=0, col="gray45")
  #arrows(q[,2], Q.ci[1,2,], q[,2],  Q.ci[2,2,], length=0, col="gray45")
  abline(v=0.1)
  abline(v=0.9)
  pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
  points(results_meta_middle2020_juv$q,results_meta_middle2020_juv$big_Q,col=pal[as.factor(hybrid_status_middle2020_juv)],pch=21,cex=1.2,bg=adjustcolor(pal[as.factor(hybrid_status_middle2020_juv)],alpha.f = 0.6))
  plot(1:10,1:10, type="n", xlab="", ylab="", axes=F)
  legend("center",legend=unique(hybrid_status_middle2020_juv),text.width=2.3,x.intersp = 0.3,col=pal[as.factor(unique(hybrid_status_middle2020_juv))],pt.bg = adjustcolor(pal[as.factor(unique(hybrid_status_middle2020_juv))],alpha.f = 0.8),pch=21,bty="n",ncol=2)
}

#Both
{par(pty="s",mfrow=c(1,3),mar=c(2,3,3,1))
  plot(results_meta_middle2020$q,results_meta_middle2020$big_Q, xlim=c(0,1), ylim=c(0,1),xlab="Proportion of YSC ancestry (q)",ylab="Interspecific ancestry (Q)",type="n",axes=F,main="Middle Creek 2020 adults")
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  #arrows(results_meta_middle2020$q.ci.width[1,2,], Q[,2], q.ci[2,2,], Q[,2], length=0, col="gray45")
  #arrows(q[,2], Q.ci[1,2,], q[,2],  Q.ci[2,2,], length=0, col="gray45")
  abline(v=0.1)
  abline(v=0.9)
  pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
  points(results_meta_middle2020$q,results_meta_middle2020$big_Q,col=pal[as.factor(hybrid_status_middle2020)],pch=21,cex=1.8,bg=adjustcolor(pal[as.factor(hybrid_status_middle2020)],alpha.f = 0.6))
  
  plot(results_meta_middle2020_juv$q,results_meta_middle2020_juv$big_Q, xlim=c(0,1), ylim=c(0,1),xlab="Proportion of YSC ancestry (q)",ylab="Interspecific ancestry (Q)",type="n",axes=F, main="Middle Creek 2020 juveniles")
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  #arrows(results_meta_middle2020$q.ci.width[1,2,], Q[,2], q.ci[2,2,], Q[,2], length=0, col="gray45")
  #arrows(q[,2], Q.ci[1,2,], q[,2],  Q.ci[2,2,], length=0, col="gray45")
  abline(v=0.1)
  abline(v=0.9)
  pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
  points(results_meta_middle2020_juv$q,results_meta_middle2020_juv$big_Q,col=pal[as.factor(hybrid_status_middle2020_juv)],pch=21,cex=1.8,bg=adjustcolor(pal[as.factor(hybrid_status_middle2020_juv)],alpha.f = 0.6))
  plot(1:10,1:10, type="n", xlab="", ylab="", axes=F)
  legend("center",legend=unique(hybrid_status_middle2020_juv),text.width=2.3,x.intersp = 0.3,col=pal[as.factor(unique(hybrid_status_middle2020_juv))],pt.bg = adjustcolor(pal[as.factor(unique(hybrid_status_middle2020_juv))],alpha.f = 0.8),cex=1.5,pch=21,bty="n",ncol=2)
}


#plot of 2019 Trout Creek & 2020 Middle Creek, adults & juvs
#triangle plots first
library(ggplot2)
meta20 <- read.csv("./NFS20_results_meta.csv")

hybrid_status <- seq(1:nrow(meta20))
for (i in 1:length(hybrid_status)){
  if (meta20$q[i] < 0.1){hybrid_status[i] <- "RBT"}
  if (meta20$q[i] > 0.9){hybrid_status[i] <- "YCT"}
  if (meta20$q[i] < 0.6 & meta20$q[i] > 0.4){
    if (meta20$big_Q[i] > 0.8){hybrid_status[i] <- "F1"}
    if (meta20$big_Q[i] > 0.4 & meta20$big_Q[i] < 0.6){hybrid_status[i] <- "F2"}
  }
  if (hybrid_status[i] == i){
    calc <- meta20$big_Q[i] - 2*meta20$q[i]
    if (calc > -0.1 & calc < 0.1){hybrid_status[i] <- "BC RBT"}
    if (hybrid_status[i] == i){
      calc <- meta20$big_Q[i] + 2*meta20$q[i]
      if (calc > 1.9 & calc < 2.1){hybrid_status[i] <- "BC YCT"}
    }
  }
  if (hybrid_status[i] == i){hybrid_status[i] <- "Other"}
}
hybrid_status <- factor(hybrid_status,c("RBT","BC RBT","F1","F2","Other","BC YCT","YCT"))
meta20$hyb_stat <- hybrid_status
adults20 <- meta20[which(!is.na(meta20$Sex) & c(1:nrow(meta20)) %in% grep("WCR2020",meta20$Gen_ID)),]
juv20 <- meta20[which(is.na(meta20$Sex) & c(1:nrow(meta20)) %in% grep("WCR2020", meta20$Gen_ID)),]

meta19 <- read.csv("./../../summer_19/parentage/results_meta_offspring.csv")
hybrid_status <- seq(1:nrow(meta19))
for (i in 1:length(hybrid_status)){
  if (meta19$q[i] < 0.1){hybrid_status[i] <- "RBT"}
  if (meta19$q[i] > 0.9){hybrid_status[i] <- "YCT"}
  if (meta19$q[i] < 0.6 & meta19$q[i] > 0.4){
    if (meta19$Q[i] > 0.8){hybrid_status[i] <- "F1"}
    if (meta19$Q[i] > 0.4 & meta19$Q[i] < 0.6){hybrid_status[i] <- "F2"}
  }
  if (hybrid_status[i] == i){
    calc <- meta19$Q[i] - 2*meta19$q[i]
    if (calc > -0.1 & calc < 0.1){hybrid_status[i] <- "BC RBT"}
    if (hybrid_status[i] == i){
      calc <- meta19$Q[i] + 2*meta19$q[i]
      if (calc > 1.9 & calc < 2.1){hybrid_status[i] <- "BC YCT"}
    }
  }
  if (hybrid_status[i] == i){hybrid_status[i] <- "Other"}
}
hybrid_status <- factor(hybrid_status,c("RBT","BC RBT","F1","F2","Other","BC YCT","YCT"))
meta19$hyb_stat <- hybrid_status
adults19 <- meta19[which(!is.na(meta19$sex)),]
juv19 <- meta19[which(is.na(meta19$sex)),]


{par(pty="s",mfcol=c(2,2),mar=c(3,2,2,0),oma=c(2,4.5,1,2),xpd=TRUE)
  pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
  plot(adults19$q,adults19$Q, xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",type="n",axes=F)
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  arrows(adults19$q, adults19$Q, adults19$q - (adults19$q.ci.width/2), adults19$Q, length=0, col="gray45") #x - left arrow
  arrows(adults19$q, adults19$Q, adults19$q + (adults19$q.ci.width/2), adults19$Q, length=0, col="gray45") #x - right arrow
  arrows(adults19$q, adults19$Q, adults19$q, adults19$Q - (adults19$Q.ci.width/2), length=0, col="gray45") #y - up arrow
  arrows(adults19$q, adults19$Q, adults19$q, adults19$Q + (adults19$Q.ci.width/2), length=0, col="gray45") #y - down arrow
  points(adults19$q,adults19$Q,col=pal[as.factor(adults19$hyb_stat)],pch=21,cex=1.8,bg=adjustcolor(pal[as.factor(adults19$hyb_stat)],alpha.f = 0.6))
  text(-0.38,0.5,labels="Adults",cex = 2,xpd=NA)
  mtext("Trout Creek 2019",side=3,line=0.5,cex=1.8)
  
  plot(juv19$q,juv19$Q, xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",type="n",axes=F)
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  arrows(juv19$q, juv19$Q, juv19$q - (juv19$q.ci.width/2), juv19$Q, length=0, col="gray45") #x - left arrow
  arrows(juv19$q, juv19$Q, juv19$q + (juv19$q.ci.width/2), juv19$Q, length=0, col="gray45") #x - right arrow
  arrows(juv19$q, juv19$Q, juv19$q, juv19$Q - (juv19$Q.ci.width/2), length=0, col="gray45") #y - up arrow
  arrows(juv19$q, juv19$Q, juv19$q, juv19$Q + (juv19$Q.ci.width/2), length=0, col="gray45") #y - down arrow
  points(juv19$q,juv19$Q,col=pal[as.factor(juv19$hyb_stat)],pch=21,cex=1.8,bg=adjustcolor(pal[as.factor(juv19$hyb_stat)],alpha.f = 0.6))
  text(-0.38,0.5,labels="Juveniles",cex = 2,xpd=NA)
  
  
  plot(adults20$q,adults20$big_Q, xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",type="n",axes=F)
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  arrows(adults20$q, adults20$big_Q, adults20$q - (adults20$q_ci_width/2), adults20$big_Q, length=0, col="gray45") #x - left arrow
  arrows(adults20$q, adults20$big_Q, adults20$q + (adults20$q_ci_width/2), adults20$big_Q, length=0, col="gray45") #x - right arrow
  arrows(adults20$q, adults20$big_Q, adults20$q, adults20$big_Q - (adults20$Q_ci_width/2), length=0, col="gray45") #y - up arrow
  arrows(adults20$q, adults20$big_Q, adults20$q, adults20$big_Q + (adults20$Q_ci_width/2), length=0, col="gray45") #y - down arrow
  points(adults20$q,adults20$big_Q,col=pal[as.factor(adults20$hyb_stat)],pch=21,cex=1.8,bg=adjustcolor(pal[as.factor(adults20$hyb_stat)],alpha.f = 0.6))
  mtext("Middle Creek 2020",side=3,line=0.5,cex=1.8)
  legend(1.08,1,legend=levels(juv20$hyb_stat),text.width=2.3,x.intersp = 0.4,pt.cex = 1.8,col=pal,pt.bg = adjustcolor(pal,alpha.f = 0.8),pch=21,bty="n",ncol=1)

  
  plot(juv20$q,juv20$big_Q, xlim=c(0,1), ylim=c(0,1),xlab="",ylab="",type="n",axes=F)
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  arrows(juv20$q, juv20$big_Q, juv20$q - (juv20$q_ci_width/2), juv20$big_Q, length=0, col="gray45") #x - left arrow
  arrows(juv20$q, juv20$big_Q, juv20$q + (juv20$q_ci_width/2), juv20$big_Q, length=0, col="gray45") #x - right arrow
  arrows(juv20$q, juv20$big_Q, juv20$q, juv20$big_Q - (juv20$Q_ci_width/2), length=0, col="gray45") #y - up arrow
  arrows(juv20$q, juv20$big_Q, juv20$q, juv20$big_Q + (juv20$Q_ci_width/2), length=0, col="gray45") #y - down arrow
  points(juv20$q,juv20$big_Q,col=pal[as.factor(juv20$hyb_stat)],pch=21,cex=1.8,bg=adjustcolor(pal[as.factor(juv20$hyb_stat)],alpha.f = 0.6))
  
  mtext("Proportion interspecific ancestry (Q)",side=2,line=-6,cex=1.2,outer=T)
  mtext("Proportion YCT ancestry (q)",side=1,line=-0.5,cex=1.2,outer=T)
}


#try for ggplot with marginal histograms
library(gridExtra)
library(ggplot2)
library(grid)
library(ggExtra)
library(RGraphics)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

{pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
col_scale <- scale_color_manual(name="Hybrid status", values=pal)
fill_scale <- scale_color_manual(name="Hybrid status",values=alpha(pal,0.8))

adult19_p <- ggplot(adults19,aes(x=q,y=Q)) + geom_segment(aes(x=0,xend=0.5,y=0,yend=1),size=1.1,color="gray60") +
  geom_segment(aes(x=0.5,xend=1,y=1,yend=0),size=1.1,color="gray60") + guides(fill="none")
adult19_p <- adult19_p + geom_point(aes(color=hyb_stat,fill=hyb_stat),size=3.5)  +
  col_scale + fill_scale + theme_minimal() + coord_fixed() + theme(panel.grid.minor = element_blank())
legend <- get_legend(adult19_p)
adult19_p <- adult19_p + guides(fill="none",color="none") + annotate("text",x=0.65,y=0.1,label=paste0("n = ",nrow(adults19)))
adult19_p_hist <- ggMarginal(adult19_p,margins="x",color="gray60",type="histogram",binwidth=0.05)

juv19_p <- ggplot(juv19,aes(x=q,y=Q)) + geom_segment(aes(x=0,xend=0.5,y=0,yend=1),size=1.1,color="gray60") +
  geom_segment(aes(x=0.5,xend=1,y=1,yend=0),size=1.1,color="gray60") + guides(fill="none",color="none")
juv19_p <- juv19_p + geom_point(aes(color=hyb_stat,fill=hyb_stat),size=3.5) + labs(x="q",y="Q") +
  col_scale + fill_scale + theme_minimal() + coord_fixed() + theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.65,y=0.1,label=paste0("n = ",nrow(juv19)))
juv19_p_hist <- ggMarginal(juv19_p,margins="x",color="gray60",type="histogram",binwidth=0.05)

adult20_p <- ggplot(adults20,aes(x=q,y=big_Q)) + geom_segment(aes(x=0,xend=0.5,y=0,yend=1),size=1.1,color="gray60") +
  geom_segment(aes(x=0.5,xend=1,y=1,yend=0),size=1.1,color="gray60") + guides(fill="none",color="none")
adult20_p <- adult20_p + geom_point(aes(color=hyb_stat,fill=hyb_stat),size=3.5) + labs(x="q",y="Q") +
  col_scale + fill_scale + theme_minimal() + coord_fixed() + theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.65,y=0.1,label=paste0("n = ",nrow(adults20)))
adult20_p_hist <- ggMarginal(adult20_p,margins="x",color="gray60",type="histogram",binwidth=0.05)

juv20_p <- ggplot(juv20,aes(x=q,y=big_Q)) + geom_segment(aes(x=0,xend=0.5,y=0,yend=1),size=1.1,color="gray60") +
  geom_segment(aes(x=0.5,xend=1,y=1,yend=0),size=1.1,color="gray60") + guides(fill="none",color="none")
juv20_p <- juv20_p + geom_point(aes(color=hyb_stat,fill=hyb_stat),size=3.5) + labs(x="q",y="Q") +
  col_scale + fill_scale + theme_minimal() + coord_fixed() + theme(panel.grid.minor = element_blank()) +
  annotate("text",x=0.65,y=0.1,label=paste0("n = ",nrow(juv20)))
juv20_p_hist <- ggMarginal(juv20_p,margins="x",color="gray60",type="histogram",binwidth=0.05)

text1 <- "Trout Creek 2019"
text2 <- "Middle Creek 2020"
text3 <- "Adults"
text4 <- "Juveniles"
}

{grid.arrange(adult19_p_hist,juv19_p_hist,adult20_p_hist,juv20_p_hist,legend,textGrob(""),
             layout_matrix=cbind(c(6,1,2),c(6,3,4),c(6,5,5)),widths=c(3.4,2.3,0.5),heights=c(0.35,2.3,2.3))
grid.text(text1,x=unit(0.3,"npc"),y=unit(0.97,"npc"),gp=gpar(fontsize=24))
grid.text(text2,x=unit(0.76,"npc"),y=unit(0.97,"npc"),gp=gpar(fontsize=24))
grid.text(text3,x=unit(0.052,"npc"),y=unit(0.69,"npc"),gp=gpar(fontsize=24))
grid.text(text4,x=unit(0.058,"npc"),y=unit(0.23,"npc"),gp=gpar(fontsize=24))
}


#plot juv q through time
library(ggplot2)
library(lubridate)
library(gridExtra)

juv19$date <- as.Date(juv19$date,format="%m/%d/%Y")
sort_juv19 <- juv19[order(juv19$date),]
juv20$Date <- as.Date(juv20$Date,format="%m/%d/%y")
sort_juv20 <- juv20[order(juv20$Date),]

d1 <- ggplot(data=sort_juv19,aes(x=date,y=q)) + geom_boxplot(aes(group=date)) +
  theme_minimal() + labs(x="Date sampled",y="q",title="Trout Creek 2019") + 
  scale_x_date(breaks=as.Date(c("8/15/2019","9/1/2019","9/15/2019","10/1/2019","10/15/2019"),format="%m/%d/%Y"),
               date_labels = "%B-%d",limits=as.Date(c("8/14/2019","10/15/2019"),"%m/%d/%Y"))
d2 <- ggplot(data=sort_juv20,aes(x=Date,y=q)) + geom_boxplot(aes(group=Date)) +
  theme_minimal() + labs(x="Date sampled",y="q",title="Middle Creek 2020") +
  scale_x_date(breaks=as.Date(c("8/15/2020","9/1/2020","9/15/2020","10/1/2020","10/15/2020"),format="%m/%d/%Y"),
               date_labels = "%B-%d",limits=as.Date(c("8/14/2020","10/15/2020"),"%m/%d/%Y"))

grid.arrange(d1,d2,ncol=1)


#Fisher's exact test for more F2 indv in the juvenile generation
#Trout Creek 2019
adults19$f2 <- sapply(adults19$hyb_stat,function(x) ifelse(x == "F2","F2","other"))
juv19$f2 <- sapply(juv19$hyb_stat,function(x) ifelse(x=="F2","F2","other"))
table19 <- cbind(table(adults19$f2),table(juv19$f2))
colnames(table19) <- c("adults","juveniles")
fisher19 <- fisher.test(table19,alternative="less") #significant!

#Middle Creek 2020
adults20$f2 <- sapply(adults20$hyb_stat,function(x) ifelse(x == "F2","F2","other"))
juv20$f2 <- sapply(juv20$hyb_stat,function(x) ifelse(x=="F2","F2","other"))
table20 <- cbind(table(adults20$f2),table(juv20$f2))
colnames(table20) <- c("adults","juveniles")
fisher20 <- fisher.test(table20,alternative="less") #not significant


#try KS test for adult generations
ks.test(adults19$q,adults20$q)


#compare density across q between adults and juveniles
library(gridExtra)
library(ggplot2)
library(grid)
library(ggExtra)
library(RGraphics)
adult19_d <- density(adults19$q,from=0,to=1,adjust=0.5)
juv19_d <- density(juv19$q,from=0,to=1,adjust=0.5)
diff_d19 <- juv19_d$y - adult19_d$y
diff_d19 <- data.frame(x=juv19_d$x,y=diff_d19)

adult20_d <- density(adults20$q,from=0,to=1,adjust=0.5)
juv20_d <- density(juv20$q,from=0,to=1,adjust=0.5)
diff_d20 <- juv20_d$y - adult20_d$y
diff_d20 <- data.frame(x=juv20_d$x,y=diff_d20)


p1 <- ggplot(diff_d19,aes(x=x,y=y)) + geom_line(color="blue") + theme_minimal() +
  geom_segment(aes(x=0,xend=1,y=0,yend=0),lwd=2) + 
  xlab(element_blank()) + ylab(element_blank()) + ggtitle("Trout Creek 2019")
p2 <- ggplot(diff_d20,aes(x=x,y=y)) + geom_line(color="blue") + theme_minimal() +
  geom_segment(aes(x=0,xend=1,y=0,yend=0),lwd=2) +
  xlab(element_blank()) + ylab(element_blank()) + ggtitle("Middle Creek 2020")
tx <- "Proportion YCT ancestry (q)"
ty <- "Difference between\njuvenile and adult\npoint densities"

{grid.arrange(p1,p2,textGrob(""),ncol=2,layout_matrix=cbind(c(3,3,3),c(1,2,3)),widths=c(0.6,3.5),heights=c(2.4,2.4,0.3))
grid.text(tx,x=unit(0.6,"npc"),y=unit(0.05,"npc"),gp=gpar(fontsize=14))
grid.text(ty,x=unit(0.08,"npc"),y=unit(0.55,"npc"),gp=gpar(fontsize=14))
}
          
          
          
          
          
          
#make plot of Trout Creek adult female sizes
fem19 <- adults19[which(adults19$sex == "Female"),]
ggplot(data=fem19,aes(x=length)) + geom_histogram() + theme_minimal() +
  labs(x="Total length (mm)",y="Count")




















