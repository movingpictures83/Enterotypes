# Clear R workspace
#rm(list=ls(all=TRUE))

#setwd("/Users/vanessa/Documents/Work/Scripts/Enterotypes/")
#setwd("/lclhome/tcickovs/PluMA/plugins/Respirotypes")
KLD <- function(x,y) sum(x * log(x/y))
#Jensen-Shannon divergence
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

#Partitioning around medoids (PAM) clustering algorithm
library(cluster)
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
#  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

#Noise removal. Not always necessary
noise.removal <- function(pcframe, percent=0.01, top=NULL){
  pcframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  #print(percent)
  return(Matrix_1)
}

input <- function(inputfile) {
  pc <<- read.table(inputfile, header = TRUE, row.names=1, dec=".", sep="\t");
}


run <- function() {

#Data loading


pc.dist=dist.JSD(pc)
k=3
pc.cluster=pam.clustering(pc.dist, k)
print("CLUSTER")
print(summary(pc.cluster))
obs.silhouette=mean(silhouette(pc.cluster, pc.dist)[,3])

pc.denoized=noise.removal(pc, percent=0.01)

#Between-class analysis (BCA)
library(ade4)
obs.pca=dudi.pca(data.frame(t(pc)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(pc.cluster), scannf=F, nf=k-1) 
s.class(obs.bet$ls, fac=as.factor(pc.cluster), grid=F, col=c(4,2,3))
#s.class(obs.bet$ls, fac=as.factor(pc.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))


#Principal Coordinate Analysis (PCoA)
obs.pcoa=dudi.pco(pc.dist, scannf=F, nf=k)
s.class(obs.pcoa$li, fac=as.factor(pc.cluster), grid=F, col=c(3,2,4))
#s.class(obs.pcoa$li, fac=as.factor(pc.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

#General info
result <- pam(as.dist(pc.dist), k, diss=TRUE)
print("RESULT")
print(summary(result))

fac=as.factor(colnames(pc))
print(colnames(pc))
s.class(obs.bet$ls, fac, grid=F, clabel=0.5, col=pc.cluster)
s.class(obs.pcoa$li, fac, grid=F, clabel=0.5, col=pc.cluster)
print("PCOA")
print(summary(obs.pcoa))
print("BET")
print(summary(obs.bet))
}

output <- function(outputfile) {

}
