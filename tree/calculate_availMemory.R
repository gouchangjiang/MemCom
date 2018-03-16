#-----------------------------------------------
rm(list=ls())
for(index in 2:10){#changehere
  setwd(paste('/Users/changjiang/Documents/trees_loris/trees_random_',index,sep=""))
  maxout_minMem<-read.table(paste('name_maxout_minM',index,'.txt',sep=""),header=TRUE,stringsAsFactors = FALSE)
  temp<-lapply(1:nrow(maxout_minMem),function(x){seq(length=5,from=maxout_minMem[x,2],to=maxout_minMem[x,3])})
  result<-matrix(data = NA,nrow=nrow(maxout_minMem),ncol = 7,byrow = TRUE,dimnames = NULL)
  for(i in 1:nrow(result)){
    for(j in 3:ncol(result)){
      result[i,j]=temp[[i]][j-2]
    }
  }
  result[,1]=maxout_minMem[,1]
  result[,2]=1 #used for store beta
  write.table(result,paste("name_beta_M",index,".txt",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)
}

write.table(maxout_minMem[,c(1,3)],paste("name_Mavail_",index,".txt",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE)
#lapply(temp,write,"Available_Memory.txt",append=TRUE,ncolumns=265)
#———————————————————————————————————————————————

setwd("/Users/changjiang/Documents/real_trees_nf/")
new.folder<-"/Users/changjiang/Documents/trees_worked_nf"
file.copy(maxout_minMem[,1],new.folder)

list.all<-list.files(,pattern = "tree_*")
list.all<-list.all[1:3000]
dif<-setdiff(list.all,maxout_minMem$tree_name)
file.remove(dif)
