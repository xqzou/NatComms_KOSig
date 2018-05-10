library(Biostrings)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(reshape2)
library(proxy) # Library of similarity/dissimilarity measures for 'dist()'
library(NMF)
library(fastICA)
library(pracma)

knockoutlist<- read.table("./knockoutlist.txt", sep = "\t", header = T, as.is = T)
muttype_freq_template <- read.table("./muttype_freq_template.txt", sep = "\t", header = T, as.is = T)

MeanDistGenomes <- function(genomes_96muttype){
  allDist <- NULL
  dist_indel <- as.matrix(proxy::dist(t(genomes_96muttype), diag=FALSE, upper=FALSE, method="Euclidean",auto_convert_data_frames = FALSE))
  dist_indel[lower.tri(dist_indel,diag=TRUE)]=NA
  dist_indel <- as.data.frame(as.table(dist_indel))
  dist_indel=na.omit(dist_indel)
  dist_indel <- dist_indel[order(dist_indel[,1],dist_indel[,2]),]
  return(mean(dist_indel$Freq))
}

SdDistGenomes <- function(genomes_96muttype){
  allDist <- NULL
  dist_indel <- as.matrix(proxy::dist(t(genomes_96muttype), diag=FALSE, upper=FALSE, method="Euclidean",auto_convert_data_frames = FALSE))
  dist_indel[lower.tri(dist_indel,diag=TRUE)]=NA
  dist_indel <- as.data.frame(as.table(dist_indel))
  dist_indel=na.omit(dist_indel)
  dist_indel <- dist_indel[order(dist_indel[,1],dist_indel[,2]),]
  return(sd(dist_indel$Freq))
}



gen_muttype_selected <- function(CTsubs,selectedID){
  CTsubs <- CTsubs[CTsubs$VariantID %in% selectedID[,1],]
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs$muttype <- paste0(CTsubs$Ref, CTsubs$Alt, CTsubs$pre_context, CTsubs$rear_context)
  muttype_freq <- data.frame(table(CTsubs$muttype))
  names(muttype_freq) <- c("mut_type", "freq")
  muttype_freq$ref <-substr((muttype_freq$mut_type),1,1)
  muttype_freq$alt <-substr((muttype_freq$mut_type),2,2)
  muttype_freq$pre_context <-substr((muttype_freq$mut_type),3,3)
  muttype_freq$rear_context <-substr((muttype_freq$mut_type),4,4)
  muttype_freq$mutation <- paste0(muttype_freq$ref,">",muttype_freq$alt)
  muttype_freq$context <- paste0(muttype_freq$pre_context,muttype_freq$ref,muttype_freq$rear_context)
  muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(muttype_freq)
  
}

quick_pointplot <- function(inputdata,xcol,ycol,outputfilename){
  
  pdf(file=outputfilename, onefile=TRUE,height=5,width=6, useDingbats=FALSE)
  p <- ggplot(sample_coreindel_percentage,aes(x=xcol,y=ycol))+geom_point(size=5)
  p <- p+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p
  dev.off()
  
}
gen_muttype <- function(CTsubs){
  CTsubs_copy <- CTsubs
  
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs$muttype <- paste0(CTsubs$Ref, CTsubs$Alt, CTsubs$pre_context, CTsubs$rear_context)
  muttype_freq <- data.frame(table(CTsubs$muttype))
  names(muttype_freq) <- c("mut_type", "freq")
  muttype_freq$ref <-substr((muttype_freq$mut_type),1,1)
  muttype_freq$alt <-substr((muttype_freq$mut_type),2,2)
  muttype_freq$pre_context <-substr((muttype_freq$mut_type),3,3)
  muttype_freq$rear_context <-substr((muttype_freq$mut_type),4,4)
  muttype_freq$mutation <- paste0(muttype_freq$ref,">",muttype_freq$alt)
  muttype_freq$context <- paste0(muttype_freq$pre_context,muttype_freq$ref,muttype_freq$rear_context)
  muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(muttype_freq)
  
}

gen_muttype_selected <- function(CTsubs,selectedID){
  CTsubs <- CTsubs[CTsubs$VariantID %in% selectedID[,1],]
  CTsubs_copy <- CTsubs
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Alt <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Alt)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$pre_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$rear_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$rear_context <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$pre_context)))
  CTsubs[CTsubs$Ref %in% c("G","A"),]$Ref <- as.character(complement(DNAStringSet(CTsubs_copy[CTsubs_copy$Ref %in% c("G","A"),]$Ref)))
  
  CTsubs$muttype <- paste0(CTsubs$Ref, CTsubs$Alt, CTsubs$pre_context, CTsubs$rear_context)
  muttype_freq <- data.frame(table(CTsubs$muttype))
  names(muttype_freq) <- c("mut_type", "freq")
  muttype_freq$ref <-substr((muttype_freq$mut_type),1,1)
  muttype_freq$alt <-substr((muttype_freq$mut_type),2,2)
  muttype_freq$pre_context <-substr((muttype_freq$mut_type),3,3)
  muttype_freq$rear_context <-substr((muttype_freq$mut_type),4,4)
  muttype_freq$mutation <- paste0(muttype_freq$ref,">",muttype_freq$alt)
  muttype_freq$context <- paste0(muttype_freq$pre_context,muttype_freq$ref,muttype_freq$rear_context)
  muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
  
  #  write.table(muttype_freq,"muttype_freq.txt", sep="\t", row.names=F, col.names=T, quote=F)
  return(muttype_freq)
  
}

plot_muttype <- function(muttype_freq, outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  filename <- paste0(outputname, "_plot.pdf")
  pdf(file=filename, onefile=TRUE,width=8,height=5)
  p <- ggplot(data=muttype_freq, aes(x=muttype_freq[,1], y=muttype_freq[,2],fill=mutation,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab(colnames(muttype_freq)[1])+ylab(colnames(muttype_freq)[2])
  #+scale_x_discrete(breaks= c(muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1]),labels = as.character(muttype_freq[,"mutation"]))
  p <- p+coord_cartesian(ylim=c(0, max(muttype_freq[,2])))
  p <- p+scale_x_discrete(labels = as.character(muttype_freq[,"context"]))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}

muttype_bysamples_relfreq <- function(allsubs,outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/a_1242/14_subs_parent08_child0/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  samplelist <- data.frame(table(allsubs$Sample))
  names(samplelist) <- c("Sample","Freq")
  allsample_muttype_freq <- NULL
  for(i in 1:dim(samplelist)[1]){
    print(i)
    sample_subs <- allsubs[allsubs$Sample==as.character(samplelist[i,"Sample"]),]
    muttype_freq <- gen_muttype(sample_subs)
    muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
    muttype_freq[is.na(muttype_freq)] <- 0
    muttype_freq$Sample <- samplelist[i,"Sample"]
    muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
    allsample_muttype_freq <- rbind(allsample_muttype_freq,muttype_freq)
  }
  write.table(allsample_muttype_freq,paste0(outputname, "_muttype_freqbysample.txt"), sep="\t", row.names=F, col.names=T, quote=F)
  
  filename <- paste0(outputname, "_plot_bysample.pdf")
  pdf(file=filename, onefile=TRUE,width=25,height=15)
  p <- ggplot(data=allsample_muttype_freq, aes(x=mut_type, y=relfreq,fill=mutation,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("mutation type")+ylab("relative frequency")
  #+scale_x_discrete(breaks= c(muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1]),labels = as.character(muttype_freq[,"mutation"]))
  p <- p+coord_cartesian(ylim=c(0, 0.10))
  p <- p+scale_x_discrete(labels = as.character(allsample_muttype_freq[,"context"]))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~Sample,ncol=8)
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}

muttype_bysamples_freq <- function(allsubs,outputname){
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  muttype_freq_template <- read.table("/nfs/cancer_archive04/xz3/a_1242/14_subs_parent08_child0/muttype_freq_template.txt", sep = "\t", header = T, as.is = T)
  samplelist <- data.frame(table(allsubs$Sample))
  names(samplelist) <- c("Sample","Freq")
  allsample_muttype_freq <- NULL
  for(i in 1:dim(samplelist)[1]){
    print(i)
    sample_subs <- allsubs[allsubs$Sample==as.character(samplelist[i,"Sample"]),]
    muttype_freq <- gen_muttype(sample_subs)
    muttype_freq <- merge(muttype_freq_template,muttype_freq[,c("mut_type","freq")],by="mut_type",all=T)
    muttype_freq[is.na(muttype_freq)] <- 0
    muttype_freq$Sample <- samplelist[i,"Sample"]
    muttype_freq$relfreq <- muttype_freq$freq/sum(muttype_freq$freq)
    allsample_muttype_freq <- rbind(allsample_muttype_freq,muttype_freq)
  }
  write.table(allsample_muttype_freq,paste0(outputname, "_muttype_freqbysample.txt"), sep="\t", row.names=F, col.names=T, quote=F)
  print(max(allsample_muttype_freq$freq))
  filename <- paste0(outputname, "_plot_bysample.pdf")
  pdf(file=filename, onefile=TRUE,width=25,height=15)
  p <- ggplot(data=allsample_muttype_freq, aes(x=mut_type, y=freq,fill=mutation,width=0.5))+ geom_bar(position="dodge", stat="identity")+xlab("mutation type")+ylab("relative frequency")
  #+scale_x_discrete(breaks= c(muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1]),labels = as.character(muttype_freq[,"mutation"]))
  # p <- p+coord_cartesian(ylim=c(0, 300))
  p <- p+scale_x_discrete(labels = as.character(allsample_muttype_freq[,"context"]))
  p <- p+scale_fill_manual(values=mypalette)
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  p <- p+facet_wrap(~Sample,ncol=8,scale="free_y")
  #expand_limits(x = muttype_freq[1,1], y = 0) 
  #p <- p+annotate("text", x=muttype_freq[1,1]:muttype_freq[dim(muttype_freq)[1],1], y=muttype_freq[,"mut_type"]+ max(muttype_freq[,"mut_type"])/50, label = as.character(muttype_freq[,"mutation"]), size=3)
  print(p)
  dev.off()
  
}
plotbar <- function(basismatrix,nbasis,gtitle,outputname){
  basismatrix <- data.frame(basismatrix)
  names(basismatrix) <- paste0("sig",seq(nbasis))
  basismatrix$id <- row.names(basismatrix)
  basismatrix_melt <- melt(basismatrix,id.vars = c("id"))
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=8,width=10, useDingbats=FALSE)
  p <- list()
  for(i in 1:nbasis){
    print(i)
    basismatrix_melt_i <- basismatrix_melt[basismatrix_melt$variable==paste0("sig",i),]
    p[[i]] <- ggplot(basismatrix_melt_i,aes(x=id,y=value), environment = environment())+geom_bar(position="dodge", stat="identity")
    p[[i]] <- p[[i]]+scale_x_discrete(limits = as.character(basismatrix$id))+ggtitle(gtitle)
    p[[i]] <- p[[i]]+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
                           axis.text.y=element_text(size=10,colour = "black"),
                           axis.title.x = element_text(size=15),
                           axis.title.y = element_text(size=15),
                           plot.title = element_text(size=10),
                           panel.grid.minor.x=element_blank(),
                           panel.grid.major.x=element_blank(),
                           panel.grid.major.y = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           panel.background = element_rect(fill = "white"),
                           panel.border = element_rect(colour = "black", fill=NA))
  }
  # print(p[[1]])
  #print(p[[2]])
  #print(p[[3]])
  do.call(grid.arrange,c(p,ncol=2))
  
  dev.off()
  
}

cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}



transitfunc <- function(m,n){
  colSums(m)*n
  return(colSums(m)*n)
}
normalizefunc <- function(m){
  return(m / colSums(m)[col(m)])
}

bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}

plotexposure <- function(coefmatrix,nbasis,outputname){
  sig_distri <- data.frame(t(coefmatrix)[,1:nbasis])
  sig_distri$rowsum <- rowSums(sig_distri[,1:nbasis])
  for(i in 1:nbasis){
    sig_distri[,i] <- sig_distri[,i]/sig_distri[,"rowsum"]
  }
  sig_distri$sample <- rownames(sig_distri)
  sig_distri_melt <- melt(sig_distri, id=c("sample","rowsum"))
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=10,width=15, useDingbats=FALSE)
  p <- ggplot(sig_distri_melt,aes(x=sample,y=value,fill=variable))+geom_bar(stat="identity")+scale_colour_brewer(palette="Spectral")
  p <- p+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  print(p)
  dev.off()
  
  
}
plotbasis <- function(basismatrix,muttype,nbasis,gtitle,outputname,h,w) {
  basismatrix <- data.frame(basismatrix)
  names(basismatrix) <- paste0("sig",seq(nbasis))
  basismatrix$MutationType <- muttype
  basismatrix$Mutation <- substr(basismatrix$MutationType,3,5)
  basismatrix <- basismatrix[order(basismatrix$Mutation),]
  basismatrix_melt <- melt(basismatrix,id.vars = c("MutationType", "Mutation"))
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  p <- list()
  for(i in 1:nbasis){
    print(i)
    basismatrix_melt_i <- basismatrix_melt[basismatrix_melt$variable==paste0("sig",i),]
    p[[i]] <- ggplot(basismatrix_melt_i,aes(x=MutationType,y=value,fill=Mutation), environment = environment())+geom_bar(position="dodge", stat="identity")+scale_fill_manual(values=mypalette)
    p[[i]] <- p[[i]]+scale_x_discrete(limits = as.character(basismatrix$MutationType))+ggtitle(paste0("sig",i))
    p[[i]] <- p[[i]]+theme(axis.text.x=element_text(angle=90, vjust=0.5,size=5),
                           axis.text.y=element_text(size=10,colour = "black"),
                           axis.title.x = element_text(size=15),
                           axis.title.y = element_text(size=15),
                           plot.title = element_text(size=10),
                           panel.grid.minor.x=element_blank(),
                           panel.grid.major.x=element_blank(),
                           panel.grid.major.y = element_blank(),
                           panel.grid.minor.y = element_blank(),
                           panel.background = element_rect(fill = "white"),
                           panel.border = element_rect(colour = "black", fill=NA))
  }
  # print(p[[1]])
  #print(p[[2]])
  #print(p[[3]])
  do.call(grid.arrange,c(p,nrow=nbasis))
  
  dev.off()
  
}

