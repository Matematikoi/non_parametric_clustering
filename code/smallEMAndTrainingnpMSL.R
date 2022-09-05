library("mixtools")
library("dplyr")

totalIterations <- 250
fileName <- paste('filename_with_gene_expressions.csv', sep ="")
clusterSize <- 17; 
maxIterationSmallEM <- 5
amountOfSmallEM <- 1000

set.seed(240)


invisible(eval(parse(text=commandArgs(TRUE))))


data <- as.matrix(read.csv(fileName, row.names = 1));


smallEM <- function (data, clusterSize, maxIterationSmallEM){
  return (
    repnormmixEM(
      t(data),
      k = clusterSize,
      verb = TRUE,
      maxit = maxIterationSmallEM,
      )
  )
}

initializeWithSmallEM <- function (data,clusterSize,maxIterationSmallEM){
  bestClustering <- NULL
  for (i in seq(amountOfSmallEM)){
    currentCluster <- smallEM (data, clusterSize, maxIterationSmallEM);
    if (is.null(bestClustering)){
      bestClustering <- currentCluster
    }
    if(bestClustering$loglik > currentCluster$loglik){
      bestClustering <- currentCluster
    }
  }
  return (bestClustering);
}

getMidPoints <- function (data, post, clusterSize){
  midPoints  <- matrix(ncol = ncol(data))
  for (i in seq(clusterSize)){
    acu <- integer(ncol(data))
    cnt <- 0 
    for (j in seq(nrow(data))){
      if(which.max(post[j,]) == i ){
        acu <- acu + data[j,]
        cnt <- cnt + 1
      }
    }
    if (cnt > 0 ){
      new_point <- acu/cnt
    }else{
      new_point <- data[which.max(post[,i]),] 
    }
    midPoints <- rbind(midPoints, new_point)
    
  }
  return (midPoints[-1,])
}

start_time <- Sys.time()
startingClustering <- initializeWithSmallEM (data,clusterSize,maxIterationSmallEM)
midPoints <- getMidPoints(t(startingClustering$x),startingClustering$posterior,clusterSize)
nonDuplicatedCenters <- midPoints[!duplicated(midPoints),]
print(dim(nonDuplicatedCenters))
npClustering <- npMSL(
  data,
  nonDuplicatedCenters, 
  samebw = FALSE ,
  verb = TRUE,
  maxiter = totalIterations
)
end_time <- Sys.time()

end_time - start_time

posteriors <- npClustering$posteriors
cluster <- apply(posteriors, 1 , which.max)
name <- paste(
  cpm,
  npAlgorithm,
  "it",
  totalIterations,
  "k",
  length(unique(cluster)),
  "original-k",
  clusterSize,
  ".csv",
  sep = "_"
)
write.csv(cluster, file = paste("results/clusters",name,sep ="_"), row.names = row.names(data))
write.csv(posteriors, file = paste("results/posteriors",name,sep ="_"), row.names = row.names(data))


