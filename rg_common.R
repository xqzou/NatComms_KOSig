library(Biostrings)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(reshape2)

library(proxy) # Library of similarity/dissimilarity measures for 'dist()'
knockoutlist<- read.table("/nfs/cancer_archive04/xz3/a_1242/15_subs/knockoutlist.txt", sep = "\t", header = T, as.is = T)

bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}

plotbasis_rg <- function(basismatrix,rgsubtype,rgtype,mypalette,nbasis,gtitle,outputname,h,w) {
  basismatrix <- data.frame(basismatrix)
  names(basismatrix) <- paste0("sig",seq(nbasis))
  basismatrix$MutationType <- rgsubtype
  basismatrix$Mutation <- rgtype
  basismatrix <- basismatrix[order(basismatrix$Mutation),]
  basismatrix_melt <- melt(basismatrix,id.vars = c("MutationType", "Mutation"))
  rg_labels <- c("Deletion 0-10K","Deletion 10K-1M","Deletion > 1M","Inversion 0-10K","Inversion 10K-1M","Inversion > 1M","Tandem-dup. 0-10K","Tandem-dup. 10K-1M","Tandem-dup. > 1M", "Translocation") 
  
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  #mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  
  p <- list()
  for(i in 1:nbasis){
    print(i)
    basismatrix_melt_i <- basismatrix_melt[basismatrix_melt$variable==paste0("sig",i),]
    p[[i]] <- ggplot(basismatrix_melt_i,aes(x=MutationType,y=value,fill=Mutation), environment = environment())+geom_bar(position="dodge", stat="identity")+scale_fill_manual(values=mypalette)
    p[[i]] <- p[[i]]+scale_x_discrete(limits = as.character(basismatrix$MutationType),labels=rg_labels)+ggtitle(paste0("sig",i))+scale_y_continuous(limits=c(0,0.7),breaks=(seq(0,0.7,0.2)),labels=percent)
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
