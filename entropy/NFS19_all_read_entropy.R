#############################################################################
#subset entropy results
library(RColorBrewer)
library(rhdf5)
library(abind)

#entropy runs were split up to decrease run time
#each split had the same reference Yellowstone cutthroat trout indv.
#read in all metadata files first
meta1 <- read.csv("./NFS_split_entropy_all_meta1.csv",header=T)
meta2 <- read.csv("./NFS_split_entropy_all_meta2.csv",header=T)
meta3 <- read.csv("./NFS_split_entropy_all_meta3.csv",header=T)
meta4 <- read.csv("./NFS_split_entropy_all_meta4.csv",header=T)
meta5 <- read.csv("./NFS_split_entropy_all_meta5.csv",header=T)
meta6 <- read.csv("./NFS_split_entropy_all_meta6.csv",header=T)
meta7 <- read.csv("./NFS_split_entropy_all_meta7.csv",header=T)
meta8 <- read.csv("./NFS_split_entropy_all_meta8.csv",header=T)


for (j in 1:8){
  if (j == 1){
    meta <- meta1
  }
  if (j == 2){
    meta <- meta2
  }
  if (j == 3){
    meta <- meta3
  }
  if (j == 4){
    meta <- meta4
  }
  if (j == 5){
    meta <- meta5
  }
  if (j == 6){
    meta <- meta6
  }
  if (j == 7){
    meta <- meta7
  }
  if (j == 8){
    meta <- meta8
  }
  
  #read in q from HDF5 files
  qrep1 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep1.hdf5"),"q")
  qrep2 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep2.hdf5"),"q")
  qrep3 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep3.hdf5"),"q")
  
  #read in Q from HDF5 files
  Q_rep1 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep1.hdf5"),"Q")
  Q_rep2 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep2.hdf5"),"Q")
  Q_rep3 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep3.hdf5"),"Q")
  
  #make plots to assess convergence
  #should be quick because we provided entropy with starting values
  pdf(file=paste0("convergence_plots_Q",j,".pdf"),onefile = T)
  par(mfrow=c(2,2))
  for (i in 1:nrow(meta)) {
    plot(Q_rep1[,2,i],type="l",ylim = c(0,1),xlim=c(0,5900),ylab=meta$id[i])
    lines(Q_rep2[,2,i],col="blue")
    lines(Q_rep3[,2,i],col="red")
  }
  dev.off()
      
  pdf(file=paste0("convergence_plots_lilq",j,".pdf"),onefile = T)
  par(mfrow=c(2,2))
  for (i in 1:nrow(meta)) {
    plot(qrep1[,2,i],type="l",ylim = c(0,1),xlim=c(0,5900),ylab=meta$id[i])
    lines(qrep2[,2,i],col="blue")
    lines(qrep3[,2,i],col="red")
  }
  dev.off()
}

num_ref <- 47 #number of reference individuals
ref_ind <- as.vector(meta1$id[1:47]) #get names
clusters <- c(1,1,2,1,2,1,2,1) #to counteract cluster switching -- not necessary in newer entropy versions
ref_ind <- data.frame("id"=ref_ind,"q1"=rep(NA,num_ref),"q2"=rep(NA,num_ref),"q3"=rep(NA,num_ref),"q4"=rep(NA,num_ref),"q5"=rep(NA,num_ref),"q5"=rep(NA,num_ref),"q6"=rep(NA,num_ref),"q7"=rep(NA,num_ref),"q8"=rep(NA,num_ref))
#look at how reference q changes across runs
for (j in 1:8){
  if (j == 1){
    meta <- meta1
  }
  if (j == 2){
    meta <- meta2
  }
  if (j == 3){
    meta <- meta3
  }
  if (j == 4){
    meta <- meta4
  }
  if (j == 5){
    meta <- meta5
  }
  if (j == 6){
    meta <- meta6
  }
  if (j == 7){
    meta <- meta7
  }
  if (j == 8){
    meta <- meta8
  }
  
  rm <- which(is.na(meta$id))

  
  qrep1 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep1.hdf5"),"q")
  qrep2 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep2.hdf5"),"q")
  qrep3 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep3.hdf5"),"q")
  
  convergence <- 1000
  qrep1 <- qrep1[c(convergence:dim(qrep1)[1]),,1:num_ref]
  qrep2 <- qrep2[c(convergence:dim(qrep2)[1]),,1:num_ref]
  qrep3 <- qrep3[c(convergence:dim(qrep3)[1]),,1:num_ref]
  
  
  all_q <- abind(qrep1,qrep2,qrep3,along=1)
  
  q <- t(apply(all_q,2:3,mean))
  
  ref_ind[,j+1] <- q[,clusters[j]]
  
}
#make plot
{pdf(file="ref_across_dataset.pdf",onefile=T)
par(mfrow=c(2,2))
for (i in 1:nrow(ref_ind)){
  plot(c(1:4),ref_ind[i,2:5],ylim=c(0.85,1),main=ref_ind$id[i])
}
dev.off()
}
#consistent across runs! Yay!


#-----------------------------------------------------------------
#convergence at 1k steps (10k steps, actually, as we recorded every 10th MCMC step)
library(RColorBrewer)
library(rhdf5)
library(abind)

meta1 <- read.csv("./NFS_split_entropy_all_meta1.csv",header=T)
meta2 <- read.csv("./NFS_split_entropy_all_meta2.csv",header=T)
meta3 <- read.csv("./NFS_split_entropy_all_meta3.csv",header=T)
meta4 <- read.csv("./NFS_split_entropy_all_meta4.csv",header=T)
meta5 <- read.csv("./NFS_split_entropy_all_meta5.csv",header=T)
meta6 <- read.csv("./NFS_split_entropy_all_meta6.csv",header=T)
meta7 <- read.csv("./NFS_split_entropy_all_meta7.csv",header=T)
meta8 <- read.csv("./NFS_split_entropy_all_meta8.csv",header=T)

clusters <- c(1,1,2,1,2,1,2,1)

entropy_results <- as.data.frame(matrix(NA,nrow=1,ncol=5)) #initialize matrix
colnames(entropy_results) <- c("id","q","Q","q.ci.width","Q.ci.width")
for (j in 1:8){ #fill it in! Same as above
  if (j == 1){
    meta <- meta1
  }
  if (j == 2){
    meta <- meta2
  }
  if (j == 3){
    meta <- meta3
  }
  if (j == 4){
    meta <- meta4
  }
  if (j == 5){
    meta <- meta5
  }
  if (j == 6){
    meta <- meta6
  }
  if (j == 7){
    meta <- meta7
  }
  if (j == 8){
    meta <- meta8
  }
  
  qrep1 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep1.hdf5"),"q")
  qrep2 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep2.hdf5"),"q")
  qrep3 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep3.hdf5"),"q")
  
  Q_rep1 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep1.hdf5"),"Q")
  Q_rep2 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep2.hdf5"),"Q")
  Q_rep3 <- h5read(paste0("./NFS19_all_0.5_0.03_filt",j,"_2rep3.hdf5"),"Q")
  

  convergence <- 1000
  qrep1 <- qrep1[c(convergence:dim(qrep1)[1]),,]
  qrep2 <- qrep2[c(convergence:dim(qrep2)[1]),,]
  qrep3 <- qrep3[c(convergence:dim(qrep3)[1]),,]

  Q_rep1 <- Q_rep1[c(convergence:dim(Q_rep1)[1]),,]
  Q_rep2 <- Q_rep2[c(convergence:dim(Q_rep2)[1]),,]
  Q_rep3 <- Q_rep3[c(convergence:dim(Q_rep3)[1]),,]

  all_q <- abind(qrep1,qrep2,qrep3,along=1)
  all_Q <- abind(Q_rep1,Q_rep2,Q_rep3,along=1)

  q <- t(apply(all_q,2:3,mean))
  Q <- t(apply(all_Q,2:3,mean))

  q.ci <- apply(all_q,2:3,quantile,probs=c(0.025,0.975))
  Q.ci <- apply(all_Q,2:3,quantile,probs=c(0.025,0.975))

  q.ci_Width <- q.ci[2,2,] - q.ci[1,2,]
  Q.ci_width <- Q.ci[2,2,] - Q.ci[1,2,]
  
  temp <- as.data.frame(matrix(NA,nrow=nrow(meta)-47,ncol=5))
  colnames(temp) <- c("id","q","Q","q.ci.width","Q.ci.width")
  temp$id <- meta$id[48:nrow(meta)]
  
  temp$q <- q[c(48:nrow(meta)),clusters[j]]
  temp$Q <- Q[c(48:nrow(meta)),2]
  temp$q.ci.width <- q.ci_Width[48:nrow(meta)]
  temp$Q.ci.width <- Q.ci_width[48:nrow(meta)]
  
  entropy_results <- rbind(entropy_results,temp)
}
entropy_results <- entropy_results[-1,]

#combine the split metadata
all_meta <- rbind(meta1[48:nrow(meta1),],meta2[48:nrow(meta2),],meta3[48:nrow(meta3),],meta4[48:nrow(meta4),],meta5[48:nrow(meta5),],meta6[48:nrow(meta6),],meta7[48:nrow(meta7),],meta8[48:nrow(meta8),])
all_meta$id == entropy_results$id

#combine metadata and entropy results
results_meta <- cbind(all_meta,entropy_results)


#clean it up
num <- as.numeric(gsub("WCR2019_([[:digit:]]+)","\\1",results_meta$id))
results_meta <- results_meta[-which(is.na(num)),]
num <- num[-which(is.na(num))]


results_meta <- results_meta[order(num),]
results_meta <- results_meta[,-c(9:11)]
results_meta$id <- as.vector(results_meta$id)

#add in phenotype scores for individual traits -- not used
phenotype_meta <- read.csv("./../../metadata/YSCFieldData_2019_phenotype.csv",stringsAsFactors = F)
phenotype_meta$Genetic.ID <- gsub("WR","WCR",phenotype_meta$Genetic.ID)
phenotype_df <- phenotype_meta[match(results_meta$id,phenotype_meta$Genetic.ID),]
phenotype_df$Genetic.ID == results_meta$id

results_meta <- data.frame(results_meta[,c(1:8)],phenotype_df,results_meta[,c(9:12)])
results_meta$tributary <- gsub("Trout Creek","Trout",results_meta$tributary)
results_meta$tributary <- gsub("Middle Creek","Middle",results_meta$tributary)
write.csv(results_meta,file="./NFS19_all_0.5_0.03_results_meta.csv",row.names = F,quote = F)

YSC_meta <- results_meta[which(results_meta$q > 0.9),]
write.csv(YSC_meta,"./NFS19_ysc_results_meta.csv",row.names = F,quote=F)
#---------------------------------------------------------------------
results_meta <- read.csv("./NFS19_all_0.5_0.03_results_meta.csv",stringsAsFactors = F)

#categorize individuals into putative hybrid classes via q and Q
hybrid_status <- seq(1:nrow(results_meta))
for (i in 1:length(hybrid_status)){
  if (results_meta$q[i] < 0.1){hybrid_status[i] <- "RBT"}
  if (results_meta$q[i] > 0.9){hybrid_status[i] <- "YSC"}
  if (results_meta$q[i] < 0.6 & results_meta$q[i] > 0.4){
    if (results_meta$Q[i] > 0.8){hybrid_status[i] <- "F1"}
    if (results_meta$Q[i] > 0.4 & results_meta$Q[i] < 0.6){hybrid_status[i] <- "F2"}
  }
  if (hybrid_status[i] == i){
    calc <- results_meta$Q[i] - 2*results_meta$q[i]
    if (calc > -0.1 & calc < 0.1){hybrid_status[i] <- "BC RBT"}
    if (hybrid_status[i] == i){
      calc <- results_meta$Q[i] + 2*results_meta$q[i]
      if (calc > 1.9 & calc < 2.1){hybrid_status[i] <- "BC YSC"}
    }
  }
  if (hybrid_status[i] == i){hybrid_status[i] <- "Other"}
}
hybrid_status <- factor(hybrid_status,c("RBT","BC RBT","F1","F2","Other","BC YSC","YSC"))

#fix entry dates
library(lubridate)
dates <- as.Date(results_meta$date,"%m/%d/%y")
dates <- yday(dates)
dates <- dates - 111
dates_old <- as.Date(results_meta$date,"%m/%d/%Y")

{par(pty="s",mfrow=c(1,2),mar=c(2,3,3,1))
plot(results_meta$q,results_meta$Q, xlim=c(0,1), ylim=c(0,1),xlab="Proportion of YSC ancestry (q)",ylab="Interspecific ancestry (Q)",type="n",axes=F)
  axis(1, at=c(0,0.5,1), labels = c(0,0.5,1))
  axis(2,at=c(0,0.5,1), labels = c(0,0.5,1))
  arrows(0,0,0.5,1, length=0, col="gray")
  arrows(0.5,1,1,0, length=0, col="gray")
  #arrows(results_meta$q.ci.width[1,2,], Q[,2], q.ci[2,2,], Q[,2], length=0, col="gray45")
  #arrows(q[,2], Q.ci[1,2,], q[,2],  Q.ci[2,2,], length=0, col="gray45")
  #abline(v=0.1)
  #abline(v=0.9)
  #pal <- wesanderson::wes_palette("Zissou1",286-112,type="continuous")
  pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
  points(results_meta$q,results_meta$Q,col=pal[hybrid_status],pch=21,cex=1.2,bg=adjustcolor(pal[hybrid_status],alpha.f = 0.5))
  plot(1:10,1:10, type="n", xlab="", ylab="", axes=F)
  legend("center",legend=levels(hybrid_status),text.width=2.3,x.intersp = 0.3,col=pal,pt.bg = adjustcolor(pal,alpha.f = 0.8),pch=21,bty="n",ncol=2)
}

library(ggplot2)
res_met2 <- results_meta
res_met2$dates <- dates
res_met2 <- res_met2[-c(1:599),]
pal <- wesanderson::wes_palette("Zissou1",54,type="continuous")
ggplot(data=res_met2,aes(x=q,y=Q,alpha=0.5,size=1.2)) + geom_point(aes(x=q,y=Q,color=dates)) + scale_color_gradientn(colors = pal) +
  coord_equal() + theme_minimal() + guides(color="colorbar")

#plot of adult q and juv q
library(ggpubr)
library(gridExtra)
hybrid_loc <- "Trout"
trout_df <- results_meta[which(results_meta$tributary == hybrid_loc & !is.na(results_meta$sex)),]
pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
p <- ggplot(trout_df,aes(x=q,fill=..x..)) + geom_histogram(binwidth = 0.05) + scale_fill_gradientn(colours = pal) + theme_classic() + guides()
p2 <- p + labs(fill="",x="Prop. YSC ancestry",y="Count",title="Adults") + theme(legend.position = "none",axis.title = element_text(size=18),axis.text=element_text(size=14),plot.title=element_text(size=19))
juv_df <- results_meta[which(results_meta$tributary == hybrid_loc & is.na(results_meta$sex)),]
pal <- rev(c("#f16d22","#ffc629","#becb70","#7dcfb6","#00b2ca","#1988af","#1d4e89"))
p3 <- ggplot(juv_df,aes(x=q,fill=..x..)) + geom_histogram(binwidth = 0.05) + scale_fill_gradientn(colours = pal) + theme_classic() + guides() + ylim(0,150)
p4 <- p3 + labs(fill="",x="Prop. YSC ancestry",y="Count",title="Juveniles") + theme(legend.position = "none",axis.title = element_text(size=18),axis.text=element_text(size=14),plot.title = element_text(size=19))
grid.arrange(p2,p4,padding=unit(0.75,"line"))

#plot of juv q through time
juv_df <- results_meta[which(results_meta$tributary == hybrid_loc & is.na(results_meta$sex)),]
library(lubridate)
juv_df$yday <- yday(as.Date(juv_df$date, "%m/%d/%Y"))
juv_df$date <- as.Date(juv_df$date, "%m/%d/%Y")
ggplot(juv_df,aes(x=date,y=q,group=date)) + geom_boxplot(fill=pal[3]) + theme_classic() + guides() + labs(x="Date",y="Prop. YSC Ancestry")

summary(lm(q ~ scale(yday),data=juv_df))
