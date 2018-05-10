source("./rg_common.R")
BootstrapDiff_final_meanclone <- function(muttype_profile,cnum, subnum, bsnum,h,w,outputname){
  
  
  parentclones <- muttype_profile[,c(2,10,18,26,34,42,50,58,66)]
  parentGenome_percentage <- rowSums(parentclones)/sum(parentclones)
  
  bootstrapiClone_all <- NULL
  for(i in 1:dim(parentclones)[2]){
    clone_i <-  parentclones[,(i)]
    if(sum(clone_i)>0){
      RepsubiClone <- matrix(rep(clone_i,bsnum/9),ncol = bsnum/9)
      bootstrapiClone <- apply(RepsubiClone, 2, function(x) rmultinom(1, sum(x), x))
      bootstrapiClone_all <- cbind(bootstrapiClone_all,bootstrapiClone)
    }
  }
  
  parentclones_sig <- rowMeans(bootstrapiClone_all)/sum(rowMeans(bootstrapiClone_all))
  sel_parents_percentage_all <- NULL
  for(j in 1:bsnum){
    sel_idx <- sample(1:dim(bootstrapiClone_all)[2],7,replace=T)
    sel_parents <- bootstrapiClone_all[,sel_idx]
    sel_parents_percentage <- rowSums(sel_parents)/sum(sel_parents)
    sel_parents_percentage_all <- cbind(sel_parents_percentage_all,sel_parents_percentage)
  }
  sel_parents_percentage_final <- rowMeans(sel_parents_percentage_all)
  sel_parents_percentage_all_diff <- sel_parents_percentage_all-parentGenome_percentage
  cloneFrobeniusDist <- apply(sel_parents_percentage_all_diff,2, function(x) norm(as.matrix(x),"f"))
  cloneFrobeniusDist_threshold <- quantile(cloneFrobeniusDist, probs = 0.99)
  print(cloneFrobeniusDist_threshold)
  
  clonecentroidDist_all <- NULL
  p1_clonecentroidDist_all <- NULL
  for(i in 1:dim(knockoutlist)[1]){
    hap1_parentchild96 <- muttype_profile[,c(1,(2+TotalcloneNum*(i-1)):(1+TotalcloneNum*i))]
    subclones <- hap1_parentchild96[,(2+cnum):(1+cnum+subnum)]
    subclones_mutnum <- colSums(subclones)
    
    subclone.centroid <- rowSums(subclones)/sum(subclones)
    subclonesum_centroid <-rowSums(subclones)
    
    clonecentroidDist <- round(norm(as.matrix(parentGenome_percentage-subclone.centroid ),"f"),digits=4)
    clonecentroidDist_all <- c(clonecentroidDist_all,clonecentroidDist)
    p1_clonecentroidDist_all <- c(p1_clonecentroidDist_all,round(length(which(cloneFrobeniusDist>clonecentroidDist))/length(cloneFrobeniusDist),digits = 4))
    
    
    # Distribution of subclones
    
    bootstrapiSubClone_all <- NULL
    for(k in 1:dim(subclones)[2]){
      clone_i <-  subclones[,(k)]
      if(sum(clone_i)>0){
        RepsubiClone <- matrix(rep(clone_i,bsnum/7),ncol = bsnum/7)
        bootstrapiSubClone <- apply(RepsubiClone, 2, function(x) rmultinom(1, sum(x), x))
        bootstrapiSubClone_all <- cbind(bootstrapiSubClone_all,bootstrapiSubClone)
      }
    }
    
    
    sel_subclones_percentage_all <- NULL
    for(j in 1:bsnum){
      sel_idx <- sample(1:dim(bootstrapiSubClone_all)[2],7,replace=T)
      sel_subclones <- bootstrapiSubClone_all[,sel_idx]
      sel_subclones_percentage <- rowSums(sel_subclones)/sum(sel_subclones)
      sel_subclones_percentage_all <- cbind(sel_subclones_percentage_all,sel_subclones_percentage)
    }
    sel_subclones_percentage_final <- rowMeans(sel_subclones_percentage_all)
    sel_subclones_percentage_all_diff <- sel_subclones_percentage_all-subclone.centroid
    SubcloneFrobeniusDist <- apply(sel_subclones_percentage_all_diff,2, function(x) norm(as.matrix(x),"f"))
    SubcloneFrobeniusDist_threshold <- quantile(SubcloneFrobeniusDist, probs = 0.99)
    print(SubcloneFrobeniusDist_threshold)
    
    p1_SubclonecentroidDist_all <- round(length(which(SubcloneFrobeniusDist>clonecentroidDist))/length(SubcloneFrobeniusDist),digits = 4)
    
    totalerror <- data.frame(cbind(cloneFrobeniusDist,SubcloneFrobeniusDist),row.names = NULL)
    names(totalerror) <- c("clone","subclone")
    totalerror_melt <- melt(totalerror)
    names(totalerror_melt) <- c("flag","FrobeniusDist")
    
    binw <- 0.0005
    mypalette <- c("lightpink2","seagreen2")
    pdf(file=paste0(outputname,"_",knockoutlist[i,1],".pdf"), onefile=TRUE,width = w,height = h)
    g1 <-ggplot(totalerror_melt, aes(x=FrobeniusDist,fill=flag)) + geom_histogram(binwidth=binw,alpha=.5, position="identity")+scale_fill_manual(values=mypalette)+ggtitle(knockoutlist[i,1])
    
    g1 <- g1+ geom_vline(aes(xintercept=cloneFrobeniusDist_threshold), colour="red", linetype="dashed")
    g1 <- g1+ geom_vline(aes(xintercept=SubcloneFrobeniusDist_threshold), colour="green", linetype="dashed")
    
    g1 <- g1+ annotate("segment", x = clonecentroidDist, xend = clonecentroidDist, y = 200, yend = 0, colour="blue", size=1, arrow=arrow())
    g1 <-g1 +theme(axis.text.x=element_text(size=10,colour = "black"),
                   axis.text.y=element_text(size=10,colour = "black"),
                   plot.title = element_text(size=10),
                   panel.grid.minor.x=element_blank(),
                   panel.grid.major.x=element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   panel.border = element_rect(colour = "black", fill=NA))
    print(g1)
    dev.off()
    
  }
  
  
  
  
  
}

BootstrapSig_quantile_final_rg_bg <- function(muttype_profile,parent_muttype_profile,cnum, subnum, bsnum,h,w,outputfilename){
  
  parentclones <- parent_muttype_profile
  parentGenome_percentage <-   rowSums(parentclones)/sum(parentclones)
  
  
  bootstrapiClone_all <- NULL
  for(i in 1:dim(parentclones)[2]){
    clone_i <-  parentclones[,(i)]
    if(sum(clone_i)>0){
      RepsubiClone <- matrix(rep(clone_i,bsnum/8),ncol = bsnum/8)
      bootstrapiClone <- apply(RepsubiClone, 2, function(x) rmultinom(1, sum(x), x))
      bootstrapiClone_all <- cbind(bootstrapiClone_all,bootstrapiClone)
    }
  }
  
  parentclones_sig <- rowMeans(bootstrapiClone_all)/sum(rowMeans(bootstrapiClone_all))
  # norm(as.matrix(parentclones_sig-parentGenome_percentage),"F")
  sel_parents_percentage_all <- NULL
  for(j in 1:bsnum){
    sel_idx <- sample(1:bsnum,7,replace=T)
    sel_parents <- bootstrapiClone_all[,sel_idx]
    sel_parents_percentage <- rowSums(sel_parents)/sum(sel_parents)
    sel_parents_percentage_all <- cbind(sel_parents_percentage_all,sel_parents_percentage)
  }
  sel_parents_percentage_final <- rowMeans(sel_parents_percentage_all)
  #sel_parents_percentage_all_diff <- sel_parents_percentage_all-sel_parents_percentage_final
  sel_parents_percentage_all_diff <- sel_parents_percentage_all-parentGenome_percentage
  cloneFrobeniusDist <- apply(sel_parents_percentage_all_diff,2, function(x) norm(as.matrix(x),"f"))
  cloneFrobeniusDist_threshold <- quantile(cloneFrobeniusDist, probs = 0.99)
  
  
  hap1_parentchild96_percentage <- muttype_profile
  hap1_parentchild96_percentage[,2:(1+cnum+subnum)] <- hap1_parentchild96_percentage[,2:(1+cnum+subnum)]/colSums(hap1_parentchild96_percentage[,2:(1+cnum+subnum)])[col(hap1_parentchild96_percentage[,2:(1+cnum+subnum)])]
  hap1_parentchild96_percentage <- replace(hap1_parentchild96_percentage, is.na(hap1_parentchild96_percentage), 0)
  
  
  subclones <- muttype_profile[,(2+cnum):(1+cnum+subnum)]
  subclones_mutnum <- colSums(subclones)
  
  
  subclone_centroid <- rowSums(subclones)/sum(subclones)
  subclonesum_centroid <-rowSums(subclones)
  subclonesum <- sum(subclonesum_centroid)
  clonecentroidDist <- norm(as.matrix(parentGenome_percentage-subclone_centroid ),"f")
  
  
  
  #  each subclone profile with total subclone mutations
  subBootstrapNum <- round(bsnum/dim(subclones)[2])
  bootstrapsubiClone_all <- NULL
  for(k in 1:dim(subclones)[2]){
    subclone_i <-  subclones[,k]
    if(sum(subclone_i)>0){
      RepsubiClone <- matrix(rep(subclone_i,subBootstrapNum),ncol = subBootstrapNum)
      bootstrapsubiClone <- apply(RepsubiClone, 2, function(x) rmultinom(1, sum(x), x))
      bootstrapsubiClone_all <- cbind(bootstrapsubiClone_all,bootstrapsubiClone)
    }
  }
  
  sel_subclones_percentage_all <- NULL
  sel_subclones_all <- NULL
  for(j in 1:dim(bootstrapsubiClone_all)[2]){
    sel_idx <- sample(1:dim(bootstrapsubiClone_all)[2],7,replace=T)
    sel_subclones <- bootstrapsubiClone_all[,sel_idx]
    sel_subclones_all <- cbind(sel_subclones_all,rowSums(sel_subclones))
    sel_subclones_percentage <- rowSums(sel_subclones)/sum(sel_subclones)
    sel_subclones_percentage_all <- cbind(sel_subclones_percentage_all,sel_subclones_percentage)
  }
  #sel_subclones_percentage_final <- rowMeans(sel_subclones_percentage_all)
  #sel_parents_percentage_all_diff <- sel_parents_percentage_all-sel_parents_percentage_final
  sel_subclones_percentage_all_diff <- sel_subclones_percentage_all-subclone_centroid
  subcloneFrobeniusDist <- apply(sel_subclones_percentage_all_diff,2, function(x) norm(as.matrix(x),"f"))
  subcloneFrobeniusDist_threshold <- quantile(subcloneFrobeniusDist, probs = 0.99)
  selected_bootstrapsubClone=sel_subclones_all[,which(subcloneFrobeniusDist<=subcloneFrobeniusDist_threshold)]
  
  bootstrapsubClone_threshold <- apply(sel_subclones_percentage_all, 1, quantile, probs=0.99)
  
  #bootstrapClone_max <-apply(selected_bootstrapsubClone, 1, max)
  
  #bootstrapClone_min <- ifelse(bootstrapClone_min>bootstrapClone_threshold[1,],bootstrapClone_min,bootstrapClone_threshold[1,]) 
  #bootstrapClone_max <- ifelse(bootstrapClone_max<bootstrapClone_threshold,bootstrapClone_max,bootstrapClone_threshold) 
  
  #a=sum(parentGenome_revised[which((parentGenome_revised>final_cutoff))])
  # a_percentage <- a/sum(parentGenome_revised)
  if(length(which((parentGenome_percentage<=bootstrapsubClone_threshold)))==96 | (clonecentroidDist<=subcloneFrobeniusDist_threshold)){
    print(paste0("clone is not distint from subclones"))
    stop_flag <- 1  # clone is not distint from subclones
  } else {
    p <- 0
    stop_flag <- 0
    while(p<=subclonesum & stop_flag==0){
      print(paste0("p:",p,"/",subclonesum))
      parentExposure_count <- subclonesum-p
      KnockoutExposure_count <- p
      
      reconstructclonecentroid <- NULL
      reconstructclonetoparentDist <- NULL
      KnockoutExposure_all <- NULL
      reconstructclone_all <- NULL
      bs_parentExposureGenome_all <- NULL
      bs_reclonecentroidDist_all <- NULL
      for(nIter in 1:100){
        print(nIter)
        #bs_parentGenome <- rmultinom(1,sum(parentclones), parentGenome)
        bs_parentExposureGenome <- rmultinom(1,parentExposure_count, parentGenome_percentage)
        #sel_idx <- sample(1:bsnum,7,replace=T)
        #parentGenome_percentage <- rowSums(bootstrapiClone_all[,sel_idx])
        #bs_parentExposureGenome <- rmultinom(1,parentExposure_count, parentGenome_percentage)
        
        KnockoutExposure <- subclonesum_centroid-bs_parentExposureGenome # residue 
        
        while(length(which(KnockoutExposure<0))>0){
          KnockoutExposure[which(KnockoutExposure<0)] <- 0
          KnockoutExposure <- KnockoutExposure-rmultinom(1,(sum(KnockoutExposure)-KnockoutExposure_count), KnockoutExposure)
          
        }
        
        reconstructclone <- bs_parentExposureGenome + KnockoutExposure
        reconstructclone_all <- cbind(reconstructclone_all,reconstructclone)
        #length(which((reconstructclone>final_cutoff)))
        
        #reconstructclonecentroidDist <- c(reconstructclonecentroidDist,norm(as.matrix(reconstructclone/sum(reconstructclone)-subclone.centroid ),"f"))
        # print(paste0("reconstructclonecentroidDist:",reconstructclonecentroidDist))
        
        #reconstructclonetoparentDist <- c(reconstructclonetoparentDist,norm(as.matrix(reconstructclone/sum(reconstructclone)-parentGenome_percentage ),"f"))
        # print(paste0("reconstructclonetoparentDist:",reconstructclonetoparentDist))
        
        KnockoutExposure_all <- cbind(KnockoutExposure_all,KnockoutExposure)
        bs_parentExposureGenome_all <- cbind(bs_parentExposureGenome_all,bs_parentExposureGenome)
        
        bs_reclonecentroidDist <- norm(as.matrix(reconstructclone/sum(reconstructclone)-subclone_centroid ),"f")
        bs_reclonecentroidDist_all <- c(bs_reclonecentroidDist_all,bs_reclonecentroidDist)
      }
      
      reconstructclone_all_percentage <- reconstructclone_all/colSums(reconstructclone_all)[col(reconstructclone_all)]
      if(max(apply((bootstrapsubClone_threshold-reconstructclone_all_percentage), 2, min))<0 & min(bs_reclonecentroidDist_all)>cloneFrobeniusDist_threshold){
        p <- p+1
        # print(paste0(max(apply((bootstrapClone_max-reconstructclone_all_percentage), 2, min))))
        #  print(paste0(min(bs_reclonecentroidDist_all),"/",max(subcloneFrobeniusDist_threshold,cloneFrobeniusDist_threshold)))
      } else {
        # reconstructclonetoparentDist[which(reconstructclonecentroidDist==reconstructclonecentroidDist_min)]
        reconstructDist_distinct <- which(apply((bootstrapsubClone_threshold-reconstructclone_all_percentage), 2, min)>=0 & bs_reclonecentroidDist_all<=cloneFrobeniusDist_threshold)
        reconstructDist_distinct_num <- length(reconstructDist_distinct)
        #print(paste0("number of reconstructDist_distinct: ",reconstructDist_distinct_num))
        if(reconstructDist_distinct_num>5){
          print(paste0("number of reconstructDist_distinct: ",reconstructDist_distinct_num))
          stop_flag <- 1
          print(paste0("parentExposurePercentage:",parentExposure_count/subclonesum))
          KnockoutExposure_mean <- rowMeans(KnockoutExposure_all[,reconstructDist_distinct])
          parentExposure_mean <- rowMeans(bs_parentExposureGenome_all[,reconstructDist_distinct])
          #reconstructclonetoparentDist[which(reconstructclonecentroidDist==min(reconstructDist_distinct$toCentroid))]
          # KnockoutExposure_all[,which(reconstructclonetoparentDist==max(reconstructDist_distinct$toParent))]
          knockoutBasis <- KnockoutExposure_mean/sum(KnockoutExposure_mean)
          parentBasis <- parentExposure_mean/sum(parentExposure_mean)
          write.table(knockoutBasis,paste0(outputfilename,"_knockout_96basis_quantile500_final2.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
          write.table(KnockoutExposure_mean,paste0(outputfilename,"_knockoutExposure_quantile500_final2.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
          write.table(parentBasis,paste0(outputfilename,"_parent_96basis_quantile500_final2.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
          write.table(parentExposure_mean,paste0(outputfilename,"_parentExposure_quantile500_final2.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
          
          rg10type <- c("Deletion","Deletion","Deletion","Invertion","Invertion","Invertion","Tandem-dup","Tandem-dup","Tandem-dup","translocation")
          rg10type_palette <- c("red","sky blue","seagreen2","grey")
          # plotbasis_rg(knockoutExposure,parentchildGenome$rgsubtype,rg10type,rg10type_palette,1,paste0("rank=",2),paste0(knockoutlist[i,1],"_10exposure"),5,5)
          plotbasis_rg(knockoutBasis,parentchildGenome$rgsubtype,rg10type,rg10type_palette,1,paste0("rank=",2),paste0(outputfilename,"_knockout_10basis_mypalette"),5,5)
          plotbasis_rg(parentBasis,parentchildGenome$rgsubtype,rg10type,rg10type_palette,1,paste0("rank=",2),paste0(outputfilename,"_parent_10basis_mypalette"),5,5)
          
          #     plotbasis_rg(knockoutBasis,parentchildGenome$MutationType,1,paste0("rank=",2),paste0(outputfilename,"_knockout_96basis_quantile500_final2"),2,8)
          #      plotbasis_rg(parentBasis,parentchildGenome$MutationType,1,paste0("rank=",2),paste0(outputfilename,"_parent_96basis_quantile500_final2"),2,8)
          
        }
        p <- p+1
      }
      
      
    }
    
  }
  
  
}


hap1_parentchild10 <- read.table("./Sample_10rgtypes.txt",sep = "\t",header = T,as.is = T)
CloneNum <- 1
SubcloneNum <- 7
TotalcloneNum <- SubcloneNum + CloneNum
BootstrapNum <-100
################################################################
# Differentiate parental clone profiles and knockout profiles
################################################################

BootstrapDiff_final_meanclone(hap1_parentchild10,CloneNum,SubcloneNum,BootstrapNum,3,4,"bootstrap_10typessample_rgs")


start_col <- 2
end_col <- dim(hap1_parentchild10)[2]-1
CloneNum <- 1
SubcloneNum <- 7
TotalcloneNum <- SubcloneNum + CloneNum
BootstrapNum <- 100
parentclones <- hap1_parentchild10[,c(2,10,18,26,34,42,50,58,66)]

for(i in 1:dim(knockoutlist)[1]){
  print(i)
  knockoutgene <- knockoutlist[i,1]
  
  parentchildGenome <- hap1_parentchild10[,c(1,(start_col+TotalcloneNum*(i-1)):(1+TotalcloneNum*i))]
  BootstrapSig_quantile_final_rg_bg(parentchildGenome,parentclones,CloneNum,SubcloneNum,BootstrapNum,3,4,paste0("rgsig_final_",knockoutlist[i,1]))
}
