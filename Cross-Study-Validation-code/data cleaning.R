setwd("/Users/tangzhaoying/Desktop/data/wk12")

##################
#input the id and GB_ACC
GPL6=read.csv('GPL6602.csv',header=F)
GPL5=read.csv('GPL570.csv',header=F)


##################
#read gene expression data 10899
gse10899=read.csv('GSE10899.csv')
gse10899=data.frame(gse10899)

#merge GPL with gene expressiond data
library(sqldf)
gse10899_merge=sqldf("select gse10899.*, GPL6.V2 from gse10899, GPL6 where gse10899.ID_REF=GPL6.V1")
gse10899_merge$ID_REF=gse10899_merge$V2
gse10899_merge$V2=NULL



######################
#read gene expression data 14479
gse14479=read.csv('GSE14479.csv')
gse14479=data.frame(gse14479)

#merge GPL with gene expressiond data
gse14479_merge=sqldf("select gse14479.*, GPL5.V2 from gse14479, GPL5 where gse14479.ID_REF=GPL5.V1")
gse14479_merge$ID_REF=gse14479_merge$V2
gse14479_merge$V2=NULL


#choose the common gene in both dataset
common=toupper(intersect(gse14479_merge$ID_REF,gse10899_merge$ID_REF))

gse10899=gse10899_merge[gse10899_merge$ID_REF%in%common,]
gse14479=gse14479_merge[gse14479_merge$ID_REF%in%common,]

#transform two datasets
GSE10899=data.frame(t(gse10899[,-1]))
colnames(GSE10899)=gse10899[,1]
type=c(0,1,1,0,0,1,1,1,1,0)
GSE10899=cbind(type,GSE10899)
save(GSE10899,file='GSE10899_common.Rdata')

GSE14479=t(gse14479[,-1])
colnames(GSE14479)=gse14479[,1]
GSE14479=data.frame(GSE14479)
type=c(rep(1,16),rep(0,9))
GSE14479=cbind(type,GSE14479)
save(GSE14479,file='GSE14479_common.Rdata')
