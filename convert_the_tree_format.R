rm(list = ls())

#----------------------------------
changeSequence<-function(filename){
  raw.data<-read.table(filename)
  raw.data[,c(1,2)]<-t(apply(raw.data[,c(1,2)],1,function(x){nrow(raw.data)+1-x}))
  raw.data[1,2]<-0
  out.data<-raw.data[order(raw.data$V1),]
  write.table(out.data,file = filename,quote=FALSE,row.names = FALSE,col.names = FALSE)
}

for(i in 1:10){ #change here
  setwd(paste("/Users/changjiang/Documents/trees_loris/trees_random_",i,sep=""))
  files<-list.files(pattern = "tree_*")
  files<-as.array(files)
  apply(files,1,changeSequence)
}

#--------------------------------