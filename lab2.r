require(randomForest)
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

######### Random Forest ###########

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

features.rf100 <- randomForest(class ~ ., data=features,ntree =100, importance=TRUE, proximity=TRUE)
print(features.rf100)

features.rf600 <- randomForest(class ~ ., data=features,ntree =600, importance=TRUE, proximity=TRUE)
print(features.rf600)

features.rf100 <- randomForest(class ~ ., data=features,ntree =100, mtry=3, importance=TRUE, proximity=TRUE)

print(features.rf100)

plot(features.rf600)

parcoord(features[,1:9], var.label=TRUE,col=c("red","green")[as.numeric(features$class)])
legend("bottomright",legend=c("cancerosa","benigna"),fill=2:4)

