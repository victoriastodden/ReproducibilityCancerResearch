#Simulation of Golub paper


#load data from website
#source("http://bioconductor.org/biocLite.R")
#biocLite("golubEsets")
library(golubEsets)
data("Golub_Train")
golub=Golub_Train
golub.expression = exprs(golub)
dim(golub.expression)


#transform the matrix, let sample be in rows and gene be in columns
golub.expression.trans = t(golub.expression)
dim(golub.expression.trans)

#add AML,ALL 0,1 into dataset
golub.pheno=pData(golub)
attach(golub.pheno)
leuk.type = (ALL.AML == "AML")
table(leuk.type)

all.aml.expression = data.frame(cbind(leuk.type, golub.expression.trans))
dim(all.aml.expression)
#View(all.aml.expression)

#Test data
data("Golub_Test")
test.golub=Golub_Test
test.golub.expression = exprs(test.golub)
dim(test.golub.expression)


#transform the matrix, let sample be in rows and gene be in columns
test.golub.expression.trans = t(test.golub.expression)


#add AML,ALL 0,1 into dataset
test.golub.pheno=pData(test.golub)
attach(test.golub.pheno)
test.leuk.type = (ALL.AML == "AML")
table(test.leuk.type)

test.all.aml.expression = as.matrix(cbind(test.leuk.type, test.golub.expression.trans))
dim(test.all.aml.expression)

all.aml.expression=data.frame(all.aml.expression)
test.all.aml.expression=data.frame(test.all.aml.expression)

# First method, without k-means -------------------------------------------

#calculate mean,std of columns for different genes for the sample in class 1 and class2, respectively.
aml.mean.expression=colMeans(all.aml.expression[leuk.type == 1,], na.rm="TRUE")
aml.std.expression=apply(all.aml.expression[leuk.type == 1,], 2, sd)

all.mean.expression=colMeans(all.aml.expression[leuk.type == 0,], na.rm="TRUE")
all.std.expression=apply(all.aml.expression[leuk.type == 0,], 2, sd)

#calculate P use (mu1-mu2)/(sigma1+sigma2)
P=NULL
for (i in 2:length(all.mean.expression)){
  p=abs((aml.mean.expression[i]-all.mean.expression[i])/(aml.std.expression[i]+all.std.expression[i]))
  P=c(P,p)
}



#5 gene

gene_5=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),5)+1)])
test.gene_5=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),5)+1)])

#Parametric NaiveBayes (use package)
#install.packages("e1071")
library(e1071)


NBfit.para=naiveBayes(as.factor(gene_5$leuk.type)~.,data = gene_5)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_5)

#test error
#test error
table(test.gene_5$test.leuk.type,ytest.pred.NB.para)

#   ytest.pred.NB.para
#   0  1
# 0 19  1
# 1  3 11
#30 samples were predicted correctly

#10 gene
gene_10=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),10)+1)])
test.gene_10=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),10)+1)])

#Parametric NaiveBayes (use package)

NBfit.para=naiveBayes(as.factor(gene_10$leuk.type)~.,data = gene_10)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_10)

#test error
table(test.gene_10$test.leuk.type,ytest.pred.NB.para)
#ytest.pred.NB.para
#     0  1
# 0  19  1
# 1  1 13


#20 gene
gene_20=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),20)+1)])
test.gene_20=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),20)+1)])

#Parametric NaiveBayes (use package)

NBfit.para=naiveBayes(as.factor(gene_20$leuk.type)~.,data = gene_20)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_20)

#test error
table(test.gene_20$test.leuk.type,ytest.pred.NB.para)

# ytest.pred.NB.para
#    0  1
# 0 19  1
# 1  1 13


#30 gene
gene_30=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),30)+1)])
test.gene_30=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),30)+1)])

#Parametric NaiveBayes (use package)

NBfit.para=naiveBayes(as.factor(gene_30$leuk.type)~.,data = gene_30)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_30)

#test error
table(test.gene_30$test.leuk.type,ytest.pred.NB.para)

#50 gene
gene_50=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),50)+1)])
test.gene_50=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),50)+1)])

#Parametric NaiveBayes (use package)

NBfit.para=naiveBayes(as.factor(gene_50$leuk.type)~.,data = gene_50)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_50)

#test error
table(test.gene_50$test.leuk.type,ytest.pred.NB.para)

#70 gene
gene_70=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),70)+1)])
test.gene_70=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),70)+1)])

#Parametric NaiveBayes (use package)

NBfit.para=naiveBayes(as.factor(gene_70$leuk.type)~.,data = gene_70)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_70)

#test error
table(test.gene_70$test.leuk.type,ytest.pred.NB.para)


#90 gene
gene_90=data.frame(all.aml.expression[,c(1,head(order(P,decreasing = T),90)+1)])
test.gene_90=data.frame(test.all.aml.expression[,c(1,head(order(P,decreasing = T),90)+1)])

#Parametric NaiveBayes (use package)

NBfit.para=naiveBayes(as.factor(gene_90$leuk.type)~.,data = gene_90)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.gene_90)

#test error
table(test.gene_90$test.leuk.type,ytest.pred.NB.para)
# ytest.pred.NB.para
#    0  1
# 0 20  0
# 1  1 13



# second, with k-means ----------------------------------------------------
#transpose the matrix for genes clustering
#clust.all.aml.expression=data.frame(cbind(golub.expression,P))

#5 gene
# library(stats)
# gene_5_cluster=kmeans(clust.all.aml.expression[,1:38],centers=5,iter.max = 20, nstart = 1,
#                       algorithm = c("Hartigan-Wong"))
# 
# table(gene_5_cluster$cluster)
# 
# BN_5=rbind(cluster_5[cluster_5$gene_5_cluster.cluster==1,][which.max(cluster_5[cluster_5$gene_5_cluster.cluster==1,39]),],
#            cluster_5[cluster_5$gene_5_cluster.cluster==2,][which.max(cluster_5[cluster_5$gene_5_cluster.cluster==2,39]),],
#            cluster_5[cluster_5$gene_5_cluster.cluster==3,][which.max(cluster_5[cluster_5$gene_5_cluster.cluster==3,39]),],
#            cluster_5[cluster_5$gene_5_cluster.cluster==4,][which.max(cluster_5[cluster_5$gene_5_cluster.cluster==4,39]),],
#            cluster_5[cluster_5$gene_5_cluster.cluster==5,][which.max(cluster_5[cluster_5$gene_5_cluster.cluster==5,39]),])
# 
# BN_5_t=data.frame(t(BN_5[1:38]),leuk.type)
# NBfit.para=naiveBayes(as.factor(as.numeric(BN_5_t$leuk.type))~.,data = BN_5_t)
# names(NBfit.para)
# NBfit.para$apriori
# NBfit.para$levels
# 
# colnames(BN_5_t)
# 
# 
# #test data
# 
# #problems here
# test.BN_5_t=test.all.aml.expression[,c(colnames(test.all.aml.expression) %in% colnames(BN_5_t),7130)]
# 
# ytest.pred.NB.para=predict(NBfit.para,test.BN_5_t)
# 
# #test error
# table(test.BN_5_t$test.leuk.type,ytest.pred.NB.para)
# 

#5 gene
gene_5_cluster=kmeans(golub.expression,centers=5,iter.max = 20, nstart = 1,
                                            algorithm = c("Hartigan-Wong"))

clust=data.frame(cbind(golub.expression,gene_5_cluster$cluster))

#select one gene in cluster1 us SNR ranking



#calculate P for cluster1, and find gene1 name
cluster1=clust[clust$V39==1,]
cluster1$V39=NULL
newclust1=data.frame(cbind(t(cluster1),leuk.type))


clust1.aml.mean.expression=colMeans(newclust1[leuk.type == 1,], na.rm="TRUE")
clust1.aml.std.expression=apply(newclust1[leuk.type == 1,], 2, sd)

clust1.all.mean.expression=colMeans(newclust1[leuk.type == 0,], na.rm="TRUE")
clust1.all.std.expression=apply(newclust1[leuk.type == 0,], 2, sd)


P=NULL
for (i in 1:length(clust1.all.mean.expression)-1){
  p=abs((clust1.aml.mean.expression[i]-clust1.all.mean.expression[i])/(clust1.aml.std.expression[i]+clust1.all.std.expression[i]))
  P=c(P,p)
}

gene1=colnames(newclust1)[head(order(P,decreasing = T),1)]

#calculate P for cluster2, and find gene2 name
cluster2=clust[clust$V39==2,]
cluster2$V39=NULL
newclust2=data.frame(cbind(t(cluster2),leuk.type))

clust2.aml.mean.expression=colMeans(newclust2[leuk.type == 1,], na.rm="TRUE")
clust2.aml.std.expression=apply(newclust2[leuk.type == 1,], 2, sd)

clust2.all.mean.expression=colMeans(newclust2[leuk.type == 0,], na.rm="TRUE")
clust2.all.std.expression=apply(newclust2[leuk.type == 0,], 2, sd)


P=NULL
for (i in 1:length(clust2.all.mean.expression)-1){
  p=abs((clust2.aml.mean.expression[i]-clust2.all.mean.expression[i])/(clust2.aml.std.expression[i]+clust2.all.std.expression[i]))
  P=c(P,p)
}

gene2=colnames(newclust2)[head(order(P,decreasing = T),1)]


#calculate P for cluster3, and find gene3 name
cluster3=clust[clust$V39==3,]
cluster3$V39=NULL
newclust3=data.frame(cbind(t(cluster3),leuk.type))

clust3.aml.mean.expression=colMeans(newclust3[leuk.type == 1,], na.rm="TRUE")
clust3.aml.std.expression=apply(newclust3[leuk.type == 1,], 2, sd)

clust3.all.mean.expression=colMeans(newclust3[leuk.type == 0,], na.rm="TRUE")
clust3.all.std.expression=apply(newclust3[leuk.type == 0,], 2, sd)


P=NULL
for (i in 1:length(clust3.all.mean.expression)-1){
  p=abs((clust3.aml.mean.expression[i]-clust3.all.mean.expression[i])/(clust3.aml.std.expression[i]+clust3.all.std.expression[i]))
  P=c(P,p)
}

gene3=colnames(newclust3)[head(order(P,decreasing = T),1)]

#calculate P for cluster4, and find gene4 name
cluster4=clust[clust$V39==4,]
cluster4$V39=NULL
newclust4=data.frame(cbind(t(cluster4),leuk.type))

clust4.aml.mean.expression=colMeans(newclust4[leuk.type == 1,], na.rm="TRUE")
clust4.aml.std.expression=apply(newclust4[leuk.type == 1,], 2, sd)

clust4.all.mean.expression=colMeans(newclust4[leuk.type == 0,], na.rm="TRUE")
clust4.all.std.expression=apply(newclust4[leuk.type == 0,], 2, sd)


P=NULL
for (i in 1:length(clust4.all.mean.expression)-1){
  p=abs((clust4.aml.mean.expression[i]-clust4.all.mean.expression[i])/(clust4.aml.std.expression[i]+clust4.all.std.expression[i]))
  P=c(P,p)
}

gene4=colnames(newclust4)[head(order(P,decreasing = T),1)]

#calculate P for cluster5, and find gene5 name
cluster5=clust[clust$V39==5,]
cluster5$V39=NULL
newclust5=data.frame(cbind(t(cluster5),leuk.type))

clust5.aml.mean.expression=colMeans(newclust5[leuk.type == 1,], na.rm="TRUE")
clust5.aml.std.expression=apply(newclust5[leuk.type == 1,], 2, sd)

clust5.all.mean.expression=colMeans(newclust5[leuk.type == 0,], na.rm="TRUE")
clust5.all.std.expression=apply(newclust5[leuk.type == 0,], 2, sd)


P=NULL
for (i in 1:length(clust5.all.mean.expression)-1){
  p=abs((clust5.aml.mean.expression[i]-clust5.all.mean.expression[i])/(clust5.aml.std.expression[i]+clust5.all.std.expression[i]))
  P=c(P,p)
}

gene5=colnames(newclust5)[head(order(P,decreasing = T),1)]

Gene_5_name=c(gene1,gene2,gene3,gene4,gene5)

#5 gene_cluster

Gene_5=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_5_name])
test.Gene_5=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_5_name])

#Parametric NaiveBayes (use package)
#install.packages("e1071")
#library(e1071)


NBfit.para=naiveBayes(as.factor(Gene_5$leuk.type)~.,data = Gene_5)
names(NBfit.para)
NBfit.para$apriori
NBfit.para$levels


#Predict for test data

ytest.pred.NB.para=predict(NBfit.para,test.Gene_5)

#test error
#test error
table(test.Gene_5$test.leuk.type,ytest.pred.NB.para)

# ytest.pred.NB.para
# FALSE TRUE
# FALSE    20    0
# TRUE      2   12


#implement function for cluster SNR
k=10
gene_cluster=kmeans(golub.expression,centers=k,iter.max = 100, nstart = 1,
                      algorithm = c("Hartigan-Wong"))
clust=data.frame(cbind(golub.expression,gene_cluster$cluster))

newclust=vector("list",k)

for (i in 1:k){
  newclust[[i]]=clust[clust$V39==i,]
  newclust[[i]]$V39=NULL
  newclust[[i]]=data.frame(cbind(t( newclust[[i]]),leuk.type))
}


clust.aml.mean.expression=vector("list",k)
clust.aml.std.expression=vector("list",k)
clust.all.mean.expression=vector("list",k)
clust.all.std.expression=vector("list",k)

for (i in 1:k){
  clust.aml.mean.expression[[i]]=colMeans(newclust[[i]][leuk.type == 1,], na.rm="TRUE")
  clust.aml.std.expression[[i]]=apply(newclust[[i]][leuk.type == 1,], 2, sd)
  clust.all.mean.expression[[i]]=colMeans(newclust[[i]][leuk.type == 0,], na.rm="TRUE")
  clust.all.std.expression[[i]]=apply(newclust[[i]][leuk.type == 0,], 2, sd)
}

P=vector("list",k)
for (j in 1:k){
  for (i in 1:length(clust.all.mean.expression[[j]])-1){
    p=abs((clust.aml.mean.expression[[j]][i]-clust.all.mean.expression[[j]][i])/(clust.aml.std.expression[[j]][i]+clust.all.std.expression[[j]][i]))
    P[[j]]=c(P[[j]],p)
  }
}

Gene=vector("list",k)
Gene_name=NULL
for (i in 1:k){
  Gene[[i]]=colnames(newclust[[i]])[head(order(P[[i]],decreasing = T),1)]
  Gene_name=c(Gene_name,Gene[[i]])
}



#5 gene
Gene_5=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_5=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_5$leuk.type)~.,data = Gene_5)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_5)

#test error
table(test.Gene_5$test.leuk.type,ytest.pred.NB.para)

# ytest.pred.NB.para
# FALSE TRUE
# FALSE    20    0
# TRUE      2   12

#10 gene
Gene_10=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_10=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_10$leuk.type)~.,data = Gene_10)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_10)

#test error
table(test.Gene_10$test.leuk.type,ytest.pred.NB.para)
#ytest.pred.NB.para
# FALSE TRUE
# FALSE    20    0
# TRUE      7    7

#20 gene
Gene_20=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_20=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_20$leuk.type)~.,data = Gene_20)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_20)

#test error
table(test.Gene_20$test.leuk.type,ytest.pred.NB.para)
ytest.pred.NB.para
# FALSE TRUE
# FALSE    19    1
# TRUE      2   12

#30 gene
Gene_30=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_30=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_30$leuk.type)~.,data = Gene_30)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_30)

#test error
table(test.Gene_30$test.leuk.type,ytest.pred.NB.para)

# ytest.pred.NB.para
# FALSE TRUE
# FALSE    19    1
# TRUE      1   13


#50 gene
Gene_50=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_50=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_50$leuk.type)~.,data = Gene_50)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_50)

#test error
table(test.Gene_50$test.leuk.type,ytest.pred.NB.para)

# ytest.pred.NB.para
# FALSE TRUE
# FALSE    20    0
# TRUE      1   13

#70 gene
Gene_70=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_70=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_70$leuk.type)~.,data = Gene_70)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_70)

#test error
table(test.Gene_70$test.leuk.type,ytest.pred.NB.para)
# ytest.pred.NB.para
# FALSE TRUE
# FALSE    19    1
# TRUE      1   13

#90 gene
Gene_90=cbind(leuk.type,all.aml.expression[,colnames(all.aml.expression)%in%Gene_name])
test.Gene_90=cbind(test.leuk.type,test.all.aml.expression[,colnames(test.all.aml.expression)%in%Gene_name])
NBfit.para=naiveBayes(as.factor(Gene_90$leuk.type)~.,data = Gene_90)

#Predict for test data
ytest.pred.NB.para=predict(NBfit.para,test.Gene_90)

#test error
table(test.Gene_90$test.leuk.type,ytest.pred.NB.para)

# ytest.pred.NB.para
# FALSE TRUE
# FALSE    19    1
# TRUE      1   13