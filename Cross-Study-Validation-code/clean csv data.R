#CSV dataset
#JUDE_AML dataset
setwd("/Users/tangzhaoying/Desktop/data/wk11")
org_aml=read.csv('JUDE_AML.csv',header=F)
delet=seq(from=1,to=159,by=2)
org1=org_aml[,delet]#delete useless information
org2=org1[,-1]#delete probe

t_aml=t(org2)#transfer the matrix
colnames(t_aml)=org1[,1]#add probes as column name

jude_aml=t_aml
save(jude_aml,file='jude_aml.rdata')


#JUDE ALL dataset
org_all=read.csv('JUDE_ALL.csv',header=F)
delet=seq(from=1,to=87,by=2)
all_org1=org_all[,delet]#delete useless information
all_org2=all_org1[,-1]#delete probe
t_all=t(all_org2)#transfer the matrix
colnames(t_all)=all_org1[,1]#add probes as column name
jude_all=t_all
save(jude_all,file='jude_all.rdata')

#combined all and aml together
jude_aml=data.frame(jude_aml)
jude_all=data.frame(jude_all)
jude=rbind(jude_aml,jude_all)

#add types
aml_obs=seq(1,79)
all_obs=seq(1,43)
obs=c(aml_obs,all_obs)
jude$type=obs

save(jude,file='jude.rdata')


#GSE10899
gse10899=read.csv('GSE10899.csv')
gse10899=data.frame(gse10899)
GSE10899=data.frame(t(gse10899[,-1]))
colnames(GSE10899)=gse10899[,1]

type=c(0,1,1,0,0,1,1,1,1,0)
GSE10899=cbind(type,GSE10899)
GSE10899=GSE10899[,1:36847]
save(GSE10899,file='GSE10899.Rdata')

#GSE14417
gse14417=read.csv('GSE14417.csv')
gse14417=data.frame(gse14417)
GSE14417=t(gse14417[,-1])
colnames(GSE14417)=gse14417[,1]
GSE14417=data.frame(GSE14417)

#delete health people
GSE14417=GSE14417[c(1:16,25:33),]
type=c(rep(1,16),rep(0,9))
GSE14417=cbind(type,GSE14417)
GSE14417=GSE14417[,1:25627]
save(GSE14417,file='GSE14417.rdata')


#GSE14479
gse14479=read.csv('GSE14479.csv')
gse14479=data.frame(gse14479)
GSE14479=t(gse14479[,-1])
colnames(GSE14479)=gse14479[,1]
GSE14479=data.frame(GSE14479)
type=c(rep(1,16),rep(0,9))
GSE14479=cbind(type,GSE14479)
GSE14479=GSE14479[,1:54676]
save(GSE14479,file='GSE14479.Rdata')

#JUDE
load('JUDE.Rdata')
type=jude$type
jude=cbind(type,jude[,-12626])
save(jude,file='jude.rdata')


#Golub
load('golub.rdata')
golub=data.frame(golub)
names(golub)[1]='type'
save(golub,file='golub.rdata')
