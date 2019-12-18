source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
library(Biobase)
library(cluster)
library(mclust)
library(flexclust)
library(kohonen)
library(kernlab)
library(tidyr)
library(pracma)
source("~/Downloads/dkmpp/R/dkmpp.r")
source("~/Downloads/dkmpp/R/minmaxnorm.r")


##############################################################
## Extract data

gset <- getGEO("GSE51082", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL97", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

test_data = read.csv("~/Downloads/geo2r.cgi.txt", sep = "")
features <- as.character(test_data$ID[1:2000])

geo1 <- exprs(gset) 
geo1 <- geo1[features,]

geo11 <- pData(gset)
SampleGroup <- geo11$source_name_ch1

###############################################
## K means
AR = rep(NA,100)
for (i in 1:100){
  k1 = kmeans(t(geo1),centers = 6)
  AR[i] = adjustedRandIndex(SampleGroup,k1$cluster)
}
mean(AR)

###############################################
## PAM s

k2 = pam(t(geo1),6, metric = "euclidean")
AR_P = adjustedRandIndex(SampleGroup,k2$clustering)

###############################################
## hclust *4
AR_h = rep(NA,4)
euc.dist <- dist(t(geo1))
#euc.dist <- as.dist(1 - corMat)
h1 = hclust(euc.dist, method = "single")
AR_h[1] = adjustedRandIndex(SampleGroup,cutree(h1, k = 6))

h2 = hclust(euc.dist, method = "complete")
AR_h[2] = adjustedRandIndex(SampleGroup,cutree(h2, k = 6))

h3 = hclust(euc.dist, method = "average")
AR_h[3] = adjustedRandIndex(SampleGroup,cutree(h3, k = 6))

h4 = hclust(euc.dist, method = "centroid")
AR_h[4] = adjustedRandIndex(SampleGroup,cutree(h4, k = 6))

###############################################
## Neural Gas from flexclust
AR_NG = rep(NA,100)
for (i in 1:100){
  k1 = cclust(t(geo1),k = 6, method = "neuralgas")
  AR_NG[i] = adjustedRandIndex(SampleGroup,k1@cluster)
}
mean(AR_NG)

###############################################
## SOM from kohonen package 

som_grid <- somgrid(xdim = 5, ydim=5, topo="hexagonal")
AR_som <- rep(NA,100)
for (i in 1:100){
  som_model <- som(t(geo1), 
                   grid=som_grid, 
                   rlen=500, 
                   alpha=c(0.05,0.01), 
                   keep.data = TRUE )
  som_cluster <- cutree(hclust(dist(som_model$data[[1]])), 6)
  AR_som[i] <- adjustedRandIndex(SampleGroup,som_cluster)
}

mean(AR_som)


###############################################
## Spectral Clustering from kernlab (*100)
AR_sc <- rep(NA,100)
for (i in 1:100){
  sc <- specc(t(geo1), centers=6)
  AR_sc[i] <- adjustedRandIndex(SampleGroup,sc)
}
mean(AR_sc)

###############################################
## K means++
library(LICORS)
AR_Spp = rep(NA,100)
for (i in 1:100){
  k1 = kmeanspp(t(geo1),k = 6)
  AR_Spp[i] = adjustedRandIndex(SampleGroup,k1$cluster)
}
mean(AR_Spp)
###############################################
## dkm ++
library(vegan)
AR_dkmp = rep(NA,100)
my_data = minmaxnorm(t(geo1))
for (i in 1:100){
  k1 <- dkmpp(my_data,k = 6)
  AR_dkmp[i] = adjustedRandIndex(SampleGroup,k1)
}
mean(AR_dkmp)

###############################################

depthfunc <- function(points1, points2) {
  distanceMatrix <- matrix(NA, nrow=dim(points1)[1], ncol=dim(points2)[1])
  for(i in 1:nrow(points2)) {
    distanceMatrix[,i] <- 1/sqrt(1+rowSums(t(t(points1)-points2[i,])^2))
  }
  distanceMatrix
}
K_meansdep <- function(x, centers, distFun, nItter) {
  clusterHistory <- vector(nItter, mode="list")
  centerHistory <- vector(nItter, mode="list")
  
  for(i in 1:nItter) {
    distsToCenters <- distFun(x, centers)
    clusters <- apply(distsToCenters, 1, which.max)
    centers <- apply(x, 2, tapply, clusters, mean)
    # Saving history
    clusterHistory[[i]] <- clusters
    centerHistory[[i]] <- centers
  }
  
  list(clusters=clusterHistory, centers=centerHistory)
}

## test
ktest = t(geo1)
centers <- ktest[sample(nrow(ktest), 5),]
res <- K_meansdep(ktest, centers, depthfunc, 3)
AR_dep = adjustedRandIndex(SampleGroup,res$clusters[[1]])
mean(AR_dkmp)


