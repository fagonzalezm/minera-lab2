require(randomForest)
require(pROC)
require(MASS)

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

features <- data[,2:11]
features$class <- factor(features$class, levels=c(2,4), labels=c("benigna","cancerosa"))

data_set_size <- floor(nrow(features)/2)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(features), size = data_set_size)
# Assign the data to the correct sets
training <- features[indexes,]
validation1 <- features[-indexes,]

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


features.mod1 <- features[,-9]

set.seed(20)
features.rf100 <- randomForest(class ~ ., data=features.mod1, ntree = 100, mtry = 4, importance=TRUE, proximity=TRUE)
print(features.rf100)
round(importance(features.rf100),2)
varImpPlot(features.rf100)


library(ROCR)

indexes <- sample(1:nrow(features), size = data_set_size)

prediction_for_table <- predict(features.rf100,validation1[,-9])
table(observed=validation1[,9],predicted=prediction_for_table)

prediction_for_roc_curve <- predict(features.rf100,validation1[,-9],type="prob")

# Use pretty colours:
pretty_colours <- c("#F8766D","#00BA38")
# Specify the different classes 
classes <- levels(validation1$class)


for (i in 1:2)
{
  # Define which observations belong to class[i]
  true_values <- ifelse(validation1[,9]==classes[i],1,0)
  # Assess the performance of classifier for class[i]
  pred <- prediction(prediction_for_roc_curve[,i],true_values)
  perf <- performance(pred, "tpr", "fpr")
  if (i==1)
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i]) 
  }
  else
  {
    plot(perf,main="ROC Curve",col=pretty_colours[i],add=TRUE) 
  }
  # Calculate the AUC and print it to screen
  auc.perf <- performance(pred, measure = "auc")
  print(auc.perf@y.values)
}
