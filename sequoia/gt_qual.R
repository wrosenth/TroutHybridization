#filter vcf by genotype call quality

library(vcfR)
setwd("D:/Masters/trout/analysis/summer_19/parentage/")

vcf <- read.vcfR("./NFS19_parentage_0.9_0.15_filt.vcf")

chrom_pos <- vcf@fix[,c(1:2)]

metadata <- read.csv("./results_meta_offspring.csv",stringsAsFactors = F)

colnames(vcf@gt) <- gsub(".*aln_(WCR2019.*).sorted.bam","\\1",colnames(vcf@gt))

vcf_met <- metadata[which(metadata$id %in% colnames(vcf@gt)),]

trout_adults <- which(vcf_met$tributary == "Trout")

gt_info <- vcf@gt[,-1]
gt_info <- gt_info[,trout_adults]
#get max genotype phred-scale likelihood
pl <- apply(gt_info,c(1,2),function(x) sub("^./.:(.*):.*:.*:.*$","\\1",x))
pl <- unname(pl)

pl_to_prob <- function(phred){
  phred1 <- sub("(.*),.*,.*","\\1",phred)
  phred2 <- sub(".*,(.*),.*","\\1",phred)
  phred3 <- sub(".*,.*,(.*)","\\1",phred)
  phred <- c(as.integer(phred1),as.integer(phred2),as.integer(phred3))
  prob <- 10^(-0.1*phred)/sum(10^(-0.1*phred))
  return(prob)
}

gl_prob <- apply(pl,c(1,2),function(x) pl_to_prob(x))

#now get max genotype likelihood
gl_max <- matrix(NA,nrow=dim(gl_prob)[2],ncol=dim(gl_prob)[3])
for (i in 1:dim(gl_prob)[2]){
  for (j in 1:dim(gl_prob)[3]){
    gl_max[i,j] <- gl_prob[which.max(gl_prob[,i,j]),i,j]
  }
}

gt_filt <- gt_info
gt_filt[which(gl_max < 0.99)] <- NA
gt_filt <- gsub("[.]/[.]",NA,gt_filt)

vcf@gt <- vcf@gt[,c(1,trout_adults+1)]
if(ncol(gl_max) == (ncol(vcf@gt) - 1)){
  gl_max_mod <- cbind(rep(NA,nrow(gl_max)),gl_max) 
}

snp_miss <- apply(gt_filt,1,function(x) length(which(is.na(x)))/length(x))
ind_miss <- apply(gt_filt,2,function(x) length(which(is.na(x)))/length(x))

snps_keep <- which(snp_miss < 0.02)

#filter vcf
for (c in 2:ncol(vcf@gt)){
  vcf@gt[which(gl_max_mod[,c] < 0.99),c] <- gsub("^./.(:.*)","./.\\1",vcf@gt[which(gl_max_mod[,c] < 0.99),c])  
}

vcf@gt <- vcf@gt[,c(1,which(ind_miss < 0.5)+1)]
vcf@gt <- vcf@gt[snps_keep,]
vcf@fix <- vcf@fix[snps_keep,]
#vcf@fix <- vcf@fi[,trout_adults]
write.vcf(vcf,file="NFS19_parentage_trout_miss0.98_gt0.01.vcf.gz")

#make another filtered vcf
vcf <- read.vcfR("./NFS19_parentage_0.9_0.15_filt.vcf")
snps_keep <- which(snp_miss < 0.03)

vcf@gt <- vcf@gt[,c(1,trout_adults+1)]
for (c in 2:ncol(vcf@gt)){
  vcf@gt[which(gl_max_mod[,c] < 0.99),c] <- gsub("^./.(:.*)","./.\\1",vcf@gt[which(gl_max_mod[,c] < 0.99),c])  
}

vcf@gt <- vcf@gt[,c(1,which(ind_miss < 0.5)+1)]
vcf@gt <- vcf@gt[snps_keep,]
vcf@fix <- vcf@fix[snps_keep,]
write.vcf(vcf,file="NFS19_parentage_trout_miss0.97_gt0.01.vcf.gz")

#and another
vcf <- read.vcfR("./NFS19_parentage_0.9_0.15_filt.vcf")
snps_keep <- which(snp_miss < 0.05)

vcf@gt <- vcf@gt[,c(1,trout_adults+1)]
for (c in 2:ncol(vcf@gt)){
  vcf@gt[which(gl_max_mod[,c] < 0.99),c] <- gsub("^./.(:.*)","./.\\1",vcf@gt[which(gl_max_mod[,c] < 0.99),c])  
}

vcf@gt <- vcf@gt[,c(1,which(ind_miss < 0.5)+1)]
vcf@gt <- vcf@gt[snps_keep,]
vcf@fix <- vcf@fix[snps_keep,]
write.vcf(vcf,file="NFS19_parentage_trout_miss0.95_gt0.01.vcf.gz")

#and another
vcf <- read.vcfR("./NFS19_parentage_0.9_0.15_filt.vcf")
snps_keep <- which(snp_miss < 0.1)

vcf@gt <- vcf@gt[,c(1,trout_adults+1)]
for (c in 2:ncol(vcf@gt)){
  vcf@gt[which(gl_max_mod[,c] < 0.99),c] <- gsub("^./.(:.*)","./.\\1",vcf@gt[which(gl_max_mod[,c] < 0.99),c])  
}

vcf@gt <- vcf@gt[,c(1,which(ind_miss < 0.5)+1)]
vcf@gt <- vcf@gt[snps_keep,]
vcf@fix <- vcf@fix[snps_keep,]
write.vcf(vcf,file="NFS19_parentage_trout_miss0.9_gt0.01.vcf.gz")

#make lifehist file
vcf <- read.vcfR("./NFS19_miss0.9_maf0.15_gt0.01_ld_filt.vcf",nrows=1)
names <- colnames(vcf@gt)[-1]
names <- unname(sapply(names, function(x) sub(".*aln_(WCR2019.*).sorted.bam$","\\1",x)))

lifehist <- read.table("./lifehist.txt",header=F,sep=",",stringsAsFactors = F)
lifehist_filt <- lifehist[which(lifehist$V1 %in% names),]
lifehist_filt$V1 == names
write.table(lifehist_filt,file="./lifehist_filt.txt",sep=",",row.names = F,quote=F)

#after running plink commands on vcf files, do the following:
#lfmm format to make geno
library(LEA)
vcf2lfmm("./NFS19_miss0.98_maf0.15_gt0.01_ld_filt.vcf")
temp <- read.lfmm("./NFS19_miss0.98_maf0.15_gt0.01_ld_filt.lfmm")
vcf <- read.vcfR("./NFS19_miss0.98_maf0.15_gt0.01_ld_filt.vcf",nrows=1)
names <- colnames(vcf@gt)[-1]
names <- unname(sapply(names, function(x) sub(".*aln_(WCR2019.*).sorted.bam$","\\1",x)))
temp <- cbind(names,temp)
write.table(temp,"./NFS19_sequoia_0.98_0.15_0.01_ld_filt.txt",quote=F,col.names = F,row.names=F)

vcf2lfmm("./NFS19_miss0.97_maf0.15_gt0.01_ld_filt.vcf")
temp <- read.lfmm("./NFS19_miss0.97_maf0.15_gt0.01_ld_filt.lfmm")
vcf <- read.vcfR("./NFS19_miss0.97_maf0.15_gt0.01_ld_filt.vcf",nrows=1)
names <- colnames(vcf@gt)[-1]
names <- unname(sapply(names, function(x) sub(".*aln_(WCR2019.*).sorted.bam$","\\1",x)))
temp <- cbind(names,temp)
write.table(temp,"./NFS19_sequoia_0.97_0.15_0.01_ld_filt.txt",quote=F,col.names = F,row.names=F)

vcf2lfmm("./NFS19_miss0.95_maf0.15_gt0.01_ld_filt.vcf")
temp <- read.lfmm("./NFS19_miss0.95_maf0.15_gt0.01_ld_filt.lfmm")
vcf <- read.vcfR("./NFS19_miss0.95_maf0.15_gt0.01_ld_filt.vcf",nrows=1)
names <- colnames(vcf@gt)[-1]
names <- unname(sapply(names, function(x) sub(".*aln_(WCR2019.*).sorted.bam$","\\1",x)))
temp <- cbind(names,temp)
write.table(temp,"./NFS19_sequoia_0.95_0.15_0.01_ld_filt.txt",quote=F,col.names = F,row.names=F)

vcf2lfmm("./NFS19_miss0.9_maf0.15_gt0.01_ld_filt.vcf")
temp <- read.lfmm("./NFS19_miss0.9_maf0.15_gt0.01_ld_filt.lfmm")
vcf <- read.vcfR("./NFS19_miss0.9_maf0.15_gt0.01_ld_filt.vcf",nrows=1)
names <- colnames(vcf@gt)[-1]
names <- unname(sapply(names, function(x) sub(".*aln_(WCR2019.*).sorted.bam$","\\1",x)))
temp <- cbind(names,temp)
write.table(temp,"./NFS19_sequoia_0.9_0.15_0.01_ld_filt.txt",quote=F,col.names = F,row.names=F)
