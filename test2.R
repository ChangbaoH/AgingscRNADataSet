#############纠正
##膀胱
{
  filea <- "E:/衰老数据库结果_0/膀胱/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }

}
##肠道组织
{
  filea <- "E:/衰老数据库结果_0/肠道组织/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-mid
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-mid/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-mid/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-mid/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-mid/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-mid/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  #########################基因表达pattern
  {
    myfilepath <- paste(filea,"基因表达pattern/patten/DOWN/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'UP',x = alltypefiles)
    file.rename(alltypefiles,newname)

    myfilepath <- paste(filea,"基因表达pattern/patten/UP/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'UP',replacement = 'DOWN',x = alltypefiles)
    file.rename(alltypefiles,newname)
  }
  setwd(paste(filea,"差异分析1/mid-old/",sep = ""))
}
##肺
{
  filea <- "E:/衰老数据库结果_0/肺/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##肝脏
{
  filea <- "E:/衰老数据库结果_0/肝脏/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##睾丸
{
  filea <- "E:/衰老数据库结果_0/睾丸/"
  ##########youth-mid
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-mid/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-mid/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-mid/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-mid/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-mid/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##骨骼肌
{
  filea <- "E:/衰老数据库结果_0/骨骼肌/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##骨髓
{
  filea <- "E:/衰老数据库结果_0/骨髓/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##脑
{
  filea <- "E:/衰老数据库结果_0/脑/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-mid
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-mid/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-mid/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-mid/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-mid/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-mid/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  #########################基因表达pattern
  {
    myfilepath <- paste(filea,"基因表达pattern/pattern/DOWN/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'UP',x = alltypefiles)
    file.rename(alltypefiles,newname)
    
    myfilepath <- paste(filea,"基因表达pattern/pattern/UP/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'UP',replacement = 'DOWN',x = alltypefiles)
    file.rename(alltypefiles,newname)
  }
  setwd(paste(filea,"差异分析1/mid-old/",sep = ""))
}
##皮肤
{
  filea <- "E:/衰老数据库结果_0/皮肤/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-mid
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-mid/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-mid/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-mid/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-mid/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-mid/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  #########################基因表达pattern
  {
    myfilepath <- paste(filea,"基因表达pattern/patten/DOWN/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'UP',x = alltypefiles)
    file.rename(alltypefiles,newname)
    
    myfilepath <- paste(filea,"基因表达pattern/patten/UP/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'UP',replacement = 'DOWN',x = alltypefiles)
    file.rename(alltypefiles,newname)
  }
  setwd(paste(filea,"差异分析1/mid-old/",sep = ""))
}
##前列腺
{
  filea <- "E:/衰老数据库结果_0/前列腺/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##肾
{
  filea <- "E:/衰老数据库结果_0/肾/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##食道
{
  filea <- "E:/衰老数据库结果_0/食道/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##视网膜
{
  filea <- "E:/衰老数据库结果_0/视网膜/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##胃
{
  filea <- "E:/衰老数据库结果_0/胃/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##心脏
{
  filea <- "E:/衰老数据库结果_0/心脏/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
}
##血液
{
  filea <- "E:/衰老数据库结果_0/血液/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-mid
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-mid/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-mid/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-mid/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-mid/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-mid/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-supold
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-supold/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-supold/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-supold/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-supold/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-supold/",temp[i,3],"_DEG_youth_supold.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-supold/",temp[i,3],"_DEG_youth_supold.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-supold/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########mid-supold
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-supold/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-supold/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-supold/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-supold/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-supold/",temp[i,3],"_DEG_mid_supold.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-supold/",temp[i,3],"_DEG_mid_supold.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-supold/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########old-supold
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/old-supold/",sep = ""))){
      dir.create(paste(filea,"差异分析1/old-supold/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/old-supold/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/old-supold/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/old-supold/",temp[i,3],"_DEG_old_supold.txt",sep=""))
      png(paste0(filea,"差异分析1/old-supold/",temp[i,3],"_DEG_old_supold.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/old-supold/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  #########################基因表达pattern
  {
    myfilepath <- paste(filea,"基因表达pattern/patten/DOWN/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'UP',x = alltypefiles)
    file.rename(alltypefiles,newname)
    
    myfilepath <- paste(filea,"基因表达pattern/patten/UP/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'UP',replacement = 'DOWN',x = alltypefiles)
    file.rename(alltypefiles,newname)
  }
  setwd(paste(filea,"差异分析1/mid-old/",sep = ""))
}
##胰腺
{
  filea <- "E:/衰老数据库结果_0/胰腺/"
  ##########mid-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/mid-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/mid-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/mid-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/mid-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.txt",sep=""))
      png(paste0(filea,"差异分析1/mid-old/",temp[i,3],"_DEG_mid_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/mid-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-mid
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-mid/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-mid/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-mid/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-mid/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-mid/",temp[i,3],"_DEG_youth_mid.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-mid/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  ##########youth-old
  {
    if(!dir.exists(paste(filea,"差异分析1/",sep = ""))){
      dir.create(paste(filea,"差异分析1/",sep = ""))
    }
    if(!dir.exists(paste(filea,"差异分析1/youth-old/",sep = ""))){
      dir.create(paste(filea,"差异分析1/youth-old/",sep = ""))
    }
    filenames <- list.files(paste(filea,"差异分析/youth-old/",sep = ""))
    temp <- as.data.frame(array(NA,dim = c(length(filenames),3)))
    temp[,1] <- filenames
    for (i in 1:length(filenames)) {
      temp[i,2] <-strsplit(x=temp[i,1],split='[.]')[[1]][2]
      temp[i,3] <-strsplit(x=temp[i,1],split='[_]')[[1]][1]
    }
    temp <- temp[which(temp[,2]=="txt"),]
    for (i in 1:dim(temp)[1]) {
      deg1 <- read.table(paste(filea,"差异分析/youth-old/",temp[i,1],sep = "")) 
      deg2 <- deg1;deg2[,2] <- -1*deg1[,2];deg2[,3]<-deg1[,4];deg2[,4]<-deg1[,3];
      for (j in 1:dim(deg2)[1]) {
        if(deg1[j,6]=="UP"){deg2[j,6]="DOWN"}
        else if(deg1[j,6]=="DOWN"){deg2[j,6]="UP"}
        else{}
      }
      write.table(deg2,file=paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.txt",sep=""))
      png(paste0(filea,"差异分析1/youth-old/",temp[i,3],"_DEG_youth_old.png",sep=""),width = 1500,height = 850)
      p<- ggplot(data=DEG,aes(x=avg_log2FC, y=-log10(p_val),color=change)) +
        geom_point(alpha=0.4, size=4) +guides(colour = guide_legend(override.aes = list(size=8)))+
        xlab("log2 fold change") + ylab("-log10(P.value)") +
        ggtitle(this_tile) +
        scale_colour_manual(values = c('blue','black','red'))+
        theme_bw(base_size = 30)+
        theme(axis.title.x =element_text(size=40), axis.title.y=element_text(size=40),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              axis.text.x = element_text(size = 40,  color = "black"),
              axis.text.y = element_text(size = 40,  color = "black"),
              legend.title= element_text(size=40),
              legend.text = element_text(size=40),
              plot.title = element_text(hjust = 0.5)) ## corresponding to the levels(res$change)
      print(p)
      dev.off()
    }
    #########################富集
    myfilepath <- paste(filea,"富集分析1/GO/youth-old/",sep = "")
    setwd(myfilepath)
    
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'temp',x = alltypefiles)
    file.rename(alltypefiles,newname)
    newname1 = gsub(pattern = 'UP',replacement = 'DOWN',x = newname)
    file.rename(newname,newname1)
    newname2 = gsub(pattern = 'temp',replacement = 'UP',x = newname1)
    file.rename(newname1,newname2)
  }
  #########################基因表达pattern
  {
    myfilepath <- paste(filea,"基因表达pattern/patten/DOWN/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'DOWN',replacement = 'UP',x = alltypefiles)
    file.rename(alltypefiles,newname)
    
    myfilepath <- paste(filea,"基因表达pattern/patten/UP/",sep = "")
    setwd(myfilepath)
    alltypefiles = dir(myfilepath)
    newname = gsub(pattern = 'UP',replacement = 'DOWN',x = alltypefiles)
    file.rename(alltypefiles,newname)
  }
  setwd(paste(filea,"差异分析1/mid-old/",sep = ""))
}

