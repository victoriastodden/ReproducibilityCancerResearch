
# Implementing Logistic Regression ----------------------------------------

library(golubEsets)
data("Golub_Train")
golub=Golub_Train
golub.expression = exprs(golub)
dim(golub.expression)


#transform the matrix, let sample be in rows and gene be in columns
golub.expression.trans = t(golub.expression)
dim(golub.expression.trans)
golub.pheno=pData(golub)
attach(golub.pheno)
leuk.type = (ALL.AML == "AML")
table(leuk.type)

all.aml.expression = as.matrix(cbind(leuk.type, golub.expression.trans))
dim(all.aml.expression)
train=data.frame(all.aml.expression)


#Use gini index select important genes
library(randomForest)
rfModel= randomForest(y=train$leuk.type, x=train[,2:7130],mtry = 84, ntree=10000,do.trace=T)

sortedImp_gini = sort(rfModel$importance[,1], decreasing=T);
sortedImpt_scaled = sort(importance(rfModel, scale = T)[,1], decreasing = T )

ImpGene=names(sortedImp_gini)[1:3051]
train_Imp=train[,c(colnames(train)%in%ImpGene)]
train_Imp=cbind(leuk.type,train_Imp)
#train_Imp$leuk.type=as.factor((train_Imp$leuk.type))
##prepare for test dataset
data("Golub_Test")
test.golub=Golub_Test
test.golub.expression = exprs(test.golub)
dim(test.golub.expression)


#transform the matrix, let sample be in rows and gene be in columns
test.golub.expression.trans = t(test.golub.expression)
dim(test.golub.expression.trans)
test.golub.pheno=pData(test.golub)
attach(test.golub.pheno)
test.leuk.type = (ALL.AML == "AML")
table(test.leuk.type)

test.all.aml.expression = as.matrix(cbind(test.leuk.type, test.golub.expression.trans))
dim(test.all.aml.expression)
test=data.frame(test.all.aml.expression)

test_Imp=test[,c(colnames(test)%in%ImpGene)]
test_Imp=cbind(test.leuk.type,test_Imp)
#test_Imp$test.leuk.type=as.factor(test_Imp$test.leuk.type)        


#Use random subset building GeneLogit model
library(foreign)  

# Reading data  
read.data = function()
{
  dataset = train_Imp
  y = dataset[,1]
  x = as.matrix(dataset[,-1])
  print(dim(x))
  
  dimnames(x) = NULL   #remove the variable name
  print(dim(x))
  list(x=x, y=y)
}

# Read data
obj = read.data()
x = obj$x
y = obj$y
rm(obj)


########################################

#localFDR.cal(x, y, v=100)      #see Liao et al, Bioinformatics 20, 2694-2701
#tau = model.estimation(x, y)
result =  pena.logit(x, y, q=20, tau= .628)



read.data = function()
{
  dataset = test_Imp
  y = dataset[,1]
  x = as.matrix(dataset[,-1])
  print(dim(x))
  
  dimnames(x) = NULL 
  print(dim(x))
  list(x=x, y=y)
}

dataset.test = read.data()
x = dataset.test$x
y = dataset.test$y
rm(dataset.test)

b = result$b
r1 = result$gene.selected
xx = x[, r1]
xx = cbind(rep(1, nrow(xx)), xx)
temp = as.vector(xx %*% b) + log(20/14) - log(27/11)

fitted = plogis(temp)

plot(y, fitted)

mean((y-fitted)^2) 

ratio = sum(y)/sum(fitted)
fitted_adjusted = fitted*ratio
plot(y, fitted_adjusted)

a=ifelse(fitted >= 0.5, 1, 0)
a
table(test_Imp$test.leuk.type,a)
# a
# 0  1
# FALSE 20  0
# TRUE   0 14