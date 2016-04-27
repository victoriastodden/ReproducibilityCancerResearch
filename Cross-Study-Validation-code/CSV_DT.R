setwd("/Users/tangzhaoying/Desktop/data/wk11")

# Train on Golub data -----------------------------------------------------
load('golub.rdata')
golub=data.frame(golub)

library(randomForest)
rfModel= randomForest(y=golub$type, x=golub[,2:7130],mtry = 84, ntree=10000,do.trace=T)

sortedImp_gini = sort(rfModel$importance[,1], decreasing=T);
sortedImpt_scaled = sort(importance(rfModel, scale = T)[,1], decreasing = T )

ImpGene=names(sortedImp_gini)[1:1038]
train_Imp=golub[,c(colnames(golub)%in%ImpGene)]
type=golub$type
train_Imp=cbind(type,train_Imp)
train_Imp$type=as.factor((train_Imp$type))

library(RWeka)
DT=J48(as.factor(type)~., data = train_Imp)

#test on GSE10899
load('GSE10899.RDATA')
#colnames(GSE10899)[colSums(is.na(GSE10899)) > 0] 
rfModel= randomForest(y=GSE10899$type, x=GSE10899[,2:36847],mtry = 192, ntree=10000,do.trace=F)

sortedImp_gini = sort(rfModel$importance[,1], decreasing=T);
sortedImpt_scaled = sort(importance(rfModel, scale = T)[,1], decreasing = T )

ImpGene=names(sortedImp_gini)[1:1038]
train_Imp=GSE10899[,c(colnames(GSE10899)%in%ImpGene)]
type=GSE10899$type
train_Imp=cbind(type,train_Imp)
train_Imp$type=as.factor((train_Imp$type))

DT_pred_99=predict(DT, newdata = train_Imp)
table(test_Imp$leuk.type, DT_pred)