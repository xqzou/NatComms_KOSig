library(Biostrings)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)
library(reshape2)
library(RColorBrewer)
library(proxy) # Library of similarity/dissimilarity measures for 'dist()'
knockoutlist<- read.table("./knockoutlist.txt", sep = "\t", header = T, as.is = T)

bootstrapGenomesfun <- function(genomes){
  
  return(apply(genomes, 2, function(x) rmultinom(1, sum(x), x)))
}


plotbasis_indel <- function(basismatrix,indelsubtype,indeltype,mypalette,mypositions,mylabels,nbasis,gtitle,outputname,h,w) {
  basismatrix <- data.frame(basismatrix)
  names(basismatrix) <- paste0("sig",seq(nbasis))
  basismatrix$MutationType <- indelsubtype
  basismatrix$Mutation <- indeltype
  basismatrix <- basismatrix[order(basismatrix$Mutation),]
  basismatrix_melt <- melt(basismatrix,id.vars = c("MutationType", "Mutation"))
  pdf(file=paste0(outputname,".pdf"), onefile=TRUE,height=h,width=w, useDingbats=FALSE)
  #mypalette <- c("sky blue","black", "red", "grey", "seagreen2", "lightpink2")
  p <- list()
  for(i in 1:nbasis){
    print(i)
    basismatrix_melt_i <- basismatrix_melt[basismatrix_melt$variable==paste0("sig",i),]
    p[[i]] <- ggplot(basismatrix_melt_i,aes(x=MutationType,y=value,fill=Mutation), environment = environment())+geom_bar(position="dodge", stat="identity")+scale_fill_manual(values=mypalette)
    p[[i]] <- p[[i]]+scale_x_discrete(limits = mypositions,labels=mylabels)+ggtitle(outputname)
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

