require(randomForest)
require(pROC)
require(MASS)
library(caret)

names <- c("id", "clumpThickness", "uniformityOfCellSize", "uniformityOfCellShape",
           "marginalAdhesion", "singleEpithelialCellSize", "bareNuclei",
           "blandChromatin", "normalNucleoli", "mitoses", "class")
data <- read.table("./breast-cancer-wisconsin.data", sep=",", col.names = names)

set.seed(20)

####### Missing Values ###########
# - Como existen 16 datos que presentan missing values en la variable barnuclei y el total de datos es 699,
#   se opta por eliminar estos datos.
data.original <- data

data.n = nrow(data)
data.m = ncol(data)

for (row in 1:data.n) {
  for (col in 1:data.m) {
    if (data[row, col] == "?") {
      data[row, col] <- NA
    }
  }
}
data$bareNuclei <- as.integer(data$bareNuclei)
data <- na.omit(data)

smp_size = floor(0.632*nrow(data))

train_ind = sample(seq_len(nrow(data)),size = smp_size)

features <- data[ , 2:11]
features$class <- factor(features$class, levels=c(2,4), labels=c("benigna","cancerosa"))

train.features <- features[train_ind,]
test.features <- features[-train_ind,]




######### Random Forest ###########

set.seed(20)
features.rf <- randomForest(class ~ ., data=features, importance=TRUE, proximity=TRUE)

print(features.rf)

#Importancia: este nos permite aplicar el principio de parcimonia; regularizacion.; este nos permite ver que tanto afecta la variable en el proceso de clasificainon
round(importance(features.rf),2)

#MeanDecreaseAccuracy;  Que tanto incide que se encuentre o no la variable de clasificacion.
#MeanCecreaseGini:      Tiene que ver con las impourezas de cada nodo al hacer las divicion de variables
varImpPlot(features.rf)

#Proximidad
features.mds <- cmdscale(1 - features.rf$proximity, eig=TRUE)

#escalamiento clasico multidimencional (la que explicao el profe en catedra)

op <- par(pty="s")
pairs(cbind(features[1:9],features.mds$points),cex=0.5,gap=0,
      col=c("red","green")[as.numeric(features$class)],
      main="Breastcancer Wisconsin: Predictos and MDS of Proximity Based on RandomForest")

par(op)

print(features.mds$GOF)
MDSplot(features.rf,features$class)

set.seed(20)
features.rf100 <- randomForest(class ~ ., data=features,ntree =100, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf100)

set.seed(20)
features.rf200 <- randomForest(class ~ ., data=features,ntree =200, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf200)

set.seed(20)
features.rf600 <- randomForest(class ~ ., data=features,ntree =600, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf600)

features.rf100 <- randomForest(class ~ ., data=features,ntree =100, mtry=4, importance=TRUE, proximity=TRUE)
print(features.rf100)

plot(features.rf600)

parcoord(features[,1:9], var.label=TRUE,col=c("green","red")[as.numeric(features$class)])
legend("bottomright",legend=c("cancerosa","benigna"),fill=2:4)



###### PRUEBAS ######

### Quitando variables mitoses

features.mod1 <- features[,-9]

set.seed(20)
features.rf100 <- randomForest(class ~ ., data=features.mod1, ntree = 100, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf100)
round(importance(features.rf100),2)
varImpPlot(features.rf100)

set.seed(20)
features.rf200.2 <- randomForest(class ~ ., data=features.mod1, ntree = 200, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf200.2)

set.seed(20)
features.rf200.4 <- randomForest(class ~ ., data=features.mod1, ntree = 200, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf200.4)

set.seed(20)
features.rf600 <- randomForest(class ~ ., data=features.mod1, ntree = 600, mtry = 2, importance=TRUE, proximity=TRUE)
print(features.rf600)
round(importance(features.rf600),2)
varImpPlot(features.rf600)

### quitando adhesion con 100

features.mod2 <- features.mod1[,-4]

set.seed(20)
features.rf100 <- randomForest(class ~ ., data=features.mod2, ntree = 100, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf100)
round(importance(features.rf100),2)
varImpPlot(features.rf100)

### quitando normaNuclei y single con 100
features.mod2 <- features.mod2[,-7]
features.mod2 <- features.mod2[,-4]
set.seed(20)
features.rf100 <- randomForest(class ~ ., data=features.mod2, ntree = 100, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf100)


#Importancia: este nos permite aplicar el principio de parcimonia; regularizacion.; este nos permite ver que tanto afecta la variable en el proceso de clasificainon
round(importance(features.rf100),2)

#MeanDecreaseAccuracy;  Que tanto incide que se encuentre o no la variable de clasificacion.
#MeanCecreaseGini:      Tiene que ver con las impourezas de cada nodo al hacer las divicion de variables
varImpPlot(features.rf100)



features.rf100 <- randomForest(class ~ ., data=features.mod1,ntree =100, importance=TRUE, proximity=TRUE)
print(features.rf100)

features.rf600 <- randomForest(class ~ ., data=features.mod1,ntree =600, importance=TRUE, proximity=TRUE)
print(features.rf600)

features.rf100 <- randomForest(class ~ ., data=features.mod1,ntree =100, mtry=2, importance=TRUE, proximity=TRUE)
print(features.rf100)

features.rf100 <- randomForest(class ~ ., data=features.mod1,ntree =100, mtry=4, importance=TRUE, proximity=TRUE)
print(features.rf100)


features.mod1 <- features[,-3]
features.mod1 <- features.mod1[,-8]

set.seed(20)
features.rf100 <- randomForest(class ~ ., data=features.mod1, ntree = 200, mtry = 2, importance=TRUE, proximity=TRUE)
print(features.rf100)
round(importance(features.rf100),2)
varImpPlot(features.rf100)


pred <- predict(features.rf100, test.features)
table(test.features[,10], pred, dnn = c("benigna", "cancerosa"))

probs <- predict(features.rf100, test.features, type = "prob")
library(ROCR)
pred2 <- prediction(probs[, 2], test.features[, 10])
perf <- performance(pred2, "tpr", "fpr")
plot(perf)


pred3 <- prediction(features.rf100$votes[, 2], features[, 10])
perf3 <- performance(pred3, "tpr", "fpr")
plot(perf3)
