### Machine Learning for Cancer Genomics 
##  K-NN & Random Forest Methods with Leukemia Cancer 
## Bioconductor Installation & Data Installation 



# Acute Lymphoblastic Leukemia Data installation 


# Segregating the B-Cell Cancer Type 
B_Cell <- grep("^B",as.character(ALL$BT))

# Philadelphia Chromsome Fusion 

gene_Phila <- which(as.character(ALL$mol.biol) %in% c())
ALL_Phila <- ALL[,intersect(gene_Phila,B_Cell)]

# Following Step ONE to filter gene with little variation



ALL_Phila <- nsFilter(ALL_Phila,var.cutoff = )$eset
View(ALL_Phila)
class(ALL_Phila)

# STEP 2 , Features Standardization 
# Define the function to measure the IQR for each raw 

rowIQR <- function(eSet) {
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet,floor(*numSamp))
  upQ <- rowQ(eSet,ceiling(*numSamp))
  upQ-lowQ
}

standardize <- function(x) (x-rowMedians(x))/rowIQR(x)
ALL_Phila_exprs <- standardize(exprs(ALL_Phila))
View(ALL_Phila_exprs)

#STEP 3 , Selecting a DIstance & Similarity Metrics 
# Euclidean DIstance 
eucD <- dist(exprs(ALL_Phila))
eucM <- as.matrix(eucD)
dim(eucM)

eucM<- dist(exprs(ALL_Phila),method = "manhattan")
eucMat <- as.matrix(eucM)

# Color Paletter for the distance visualization 
install.packages("RColorBrewer")
library(RColorBrewer)

hmcol <- colorRampPalette(brewer.pal(10,"RdBu"))(256)

heatmap(eucM,
        symm = TRUE,col=hmcol,
        distfun = as.dist)

# Localizing the genes of interets in the proximity 
BiocManager::install("bioDist")
library(bioDist)

closest.top("39120_at",eucM,1)

# Spearman correlation method 
spearD <- spearman.dist(exprs(ALL_Phila))
spearM <- as.matrix(spearD)
heatmap(spearM,
        symm = TRUE,col=hmcol,
        distfun = as.dist)
closest.top("39120_at",spearM,1)

# TAU Kendall DIstance Metrics 
tau_D <- tau.dist(exprs(ALL_Phila))
tau_M <- as.matrix(tau_D)

heatmap(tau_M,
        symm = TRUE,col=hmcol,
        distfun = as.dist)
closest.top("39120_at",tau_M,1)

# Hamming Mutual Information 

Mutual_D <- MIdist(exprs(ALL_Phila))
Mutual_M <- as.matrix(Mutual_D)
heatmap(Mutual_M,
        symm = TRUE,col=hmcol,
        distfun = as.dist)
closest.top("39120_at",Mutual_M,1)

# Kullback_Liebler Distance 
KL_D <- KLdist.matrix(exprs(ALL_Phila))
KL_M <- as.matrix(KL_D)
heatmap(KL_M,
        symm = TRUE,col=hmcol,
        distfun = as.dist)
closest.top("39120_at",KL_M,1)

# Multi_Dimensional Projection 

install.packages("MASS")
library(MASS)
mass <- sammon(tau_M,trace = FALSE)
mass
plot(mass$points,
     col=ifelse(ALL_Phila$mol.biol=="BCR/ABL","black","red"),
     xlab="Dimension1",ylab = "Dimension2")
install.packages("rgl")
library(rgl)       

mass_3D <- sammon(tau_M,
                  k=3,
                  trace = FALSE)
plot3d(mass$points,
     col=ifelse(ALL_Phila$mol.biol=="BCR/ABL","black","red"),
     xlab="Dimension1",ylab = "Dimension2")

# T_Statistics 
tstat <- rowttests(ALL_Phila,"mol.biol")
tstat_order <- order(tstat$p.value, decreasing = TRUE)
ALL_Phila_mass <- ALL_Phila[tstat_order[1:60],]
dist_Man <- dist(exprs(ALL_Phila_mass),method = "euclidean")
mass_tstat <- sammon(dist_Man,trace = FALSE)
plot(mass$points,
     col=ifelse(ALL_Phila$mol.biol=="BCR/ABL","black","red"),
     xlab="Dimension1",ylab = "Dimension2")


### Supervised Machine Learning 
##  Data Spliting 

NEGS <- which(ALL_Phila$mol.biol=="NEG")
BCR <- which(ALL_Phila$mol.biol=="BCR/ABL")
sample_1 <- sample(NEGS,22,replace = FALSE)
sample_2 <- sample(BCR,22, replace = FALSE)
Training_set <- c(sample_1,sample_2)
Test_Set <- setdiff(1:80,Training_set)

# Machine Learning using the "MLInterfaces" 
BiocManager::install("MLInterfaces")
library(MLInterfaces)

# KNN Algorithm 
knn1 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               trainInd = Training_set)
confuMat(knn1)
error_rate_knn <- round((confuMat(knn1)[1,2]+confuMat(knn1)[2,1])/sum(confuMat(knn1))*100,2)

# Diagnol Linear Descriminant Analysis Algorithm 
dlda1 <- MLearn(formula=mol.biol~.,
                data = ALL_Phila,
                .method = dldaI,
                trainInd = Training_set)
confuMat(dlda1)
error_rate_dlda1 <- round((confuMat(dlda1)[1,2]+confuMat(dlda1)[2,1])/sum(confuMat(dlda1))*100,2)



### Performance Improvements & Tuning 
##  Manual / Statistical Methods 
#   T-Statistics Method on the training set 


training_ttets <- rowttests(ALL_Phila[,Training_set],"mol.biol")
order_training_ttest <- order(abs(training_ttets), decreasing = TRUE)
training_selected <- featureNames(ALL_Phila)[order_training_ttest[1:50]]

# KNN with the reduced genes 

knn2 <- MLearn(formula = mol.biol~.,
                data = ALL_Phila[training_selected,],
                .method = knnI(k=1,l=0),
                trainInd = Training_set)
confuMat(knn2)
error_rate_knn2 <- round((confuMat(knn2)[1,2]+confuMat(knn2)[2,1])/sum(confuMat(knn2))*100,2)

# Diagnal Linear Desriminant Analysis 
dlda2 <- MLearn(formula=mol.biol~.,
                data = ALL_Phila[training_selected,],
                .method = dldaI,
                trainInd = Training_set)
confuMat(dlda2)
error_rate_dlda2 <- round((confuMat(dlda2)[1,2]+confuMat(dlda2)[2,1])/sum(confuMat(dlda2))*100,2)


### Cross Validation 
##  Cross Validation applying xvalSpec
#   KNN Algorithm 

knn3 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               xvalSpec("LOO"))
               
confuMat(knn3)
# KNN with the feature selection using cross validation 

knn4 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               xvalSpec("LOO",fsFun =fs.absT(50)))
confuMat(knn4)

table(unlist(fsHistory(knn4)))
# KNN with the feature selection with lower selected features 
knn5 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               xvalSpec("LOO",fsFun =fs.absT(5)))
confuMat(knn5)
table(unlist(fsHistory(knn5)))

### Ensemble Methods & Random Forest Algorithms 
#   Install Package the random Forest 


Rf1 <- MLearn(mol.biol~.,
              ALL_Phila,
              randomForestI,
              Training_set,
              ntree=1200,
              mtry=50,
              importance=TRUE)
confuMat(Rf1,"train")
confuMat(Rf1,"test")

Rf2 <- MLearn(mol.biol~.,
              ALL_Phila,
              randomForestI,
              Training_set,
              ntree=1200,
              mtry=5,
              importance=TRUE)


confuMat(Rf2,"train")
confuMat(Rf2,"test")    # Much lower error rate 

# Feature Selection using Random Forest 
BiocManager::install("hgu95av2.db")
library(hgu95av2.db)

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

getVarImp(Rf1)
plot(getVarImp(Rf1),n=15,plat="hgu95av2",toktype="SYMBOL")

getVarImp(Rf2)
plot(getVarImp(Rf2),n=15,plat="hgu95av2",toktype="SYMBOL")

MeanDecreaseGini()
