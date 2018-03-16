
for(i in 1:10){#change here
  setwd(paste("/Users/changjiang/Documents/trees_loris/trees_random_",i,sep=""))
  files<-list.files(pattern = "*.tree.nf")
  write.table(files,"files_list.txt",quote=FALSE,row.names = FALSE, col.names = FALSE)
}

for(i in 2:10){#change here
  setwd(paste("/Users/changjiang/Documents/trees_loris/trees_random_",i,sep=""))
  files<-read.table('files_list.txt')
  keep<-read.table(paste('name_maxout_minM',i,'.txt',sep=""),header = TRUE)
  remove<-setdiff(files$V1,keep$tree_name)
  file.remove(remove)
}

library(reshape2)
for(i in 1:10){#change here
  setwd(paste("/Users/changjiang/Documents/trees_loris/trees_random_",i,sep=""))
  #setwd("~/Documents/trees_real")
  name_M<-read.table(paste("name_maxout_minM",i,".txt",sep=""),header = TRUE)
  #name_M<-read.table("name_maxout_minM.txt",header = TRUE)
  
  fun_beta<-function(x){
    tree<-read.table(x) #dump.1.3.amd.QY.case39-2136.tree.nf
    f.ave<-mean(tree$V5)#change here
    w.ave<-mean(tree$V4)#change here
    beta<-c(16,1/0.5,1,1/2,1/16)*(f.ave/w.ave)#change here
    return(beta)
  }
  
  beta<-apply(as.matrix(name_M[,1]),1,fun_beta)
  beta<-t(beta)
  name_M[,4:8]<-beta
  
  write.table(name_M,"name_m_beta.txt",quote=FALSE,row.names = FALSE, col.names = FALSE)
}

#change the tree format, move computation weight from column 5 to 4
setwd("~/Documents/trees_real")
name_M<-read.table("name_maxout_minM.txt",header = TRUE)

fun_change<-function(x){
  tree<-read.table(x)
  tree$V6<-tree$V5
  tree$V5<-tree$V4
  tree$V4<-tree$V6
  write.table(tree[,c(1:5)],x,quote=FALSE,row.names = FALSE, col.names = FALSE)
}

apply(as.matrix(name_M[,1]),1,fun_change)

