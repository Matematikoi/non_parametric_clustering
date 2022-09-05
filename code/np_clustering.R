library("mixtools")
library("dplyr")
cpm <- '1_5'
npAlgorithm <- "npMSL"
totalIterations <- 10000 #500
clusterSize <- 17; #change to 17
maxIterationSmallEM <- 5
amountOfSmallEM <- 10 #change to 10 or something like that
numberOfCores <- 10


set.seed(240)


invisible(eval(parse(text=commandArgs(TRUE))))

name = "data/1_5_cpm_mouse_log.csv"

#read the data and combine it. 
data <- read.csv(name, row.names = 1, header= TRUE)

conds <- c(
  "kidney",
  "kidney",
  "liver",
  "liver",
  "lung",
  "lung",
  "short_interstine",
  "short_interstine"
)

# data <- data[sample(nrow(data), size = 1000),]




getMidPoints <- function (cpm){
  return (readRDS("smallEMdata/npCluster_smallEM__1_5_npMSL_.RDS"))
}

start_time <- Sys.time()
nonDuplicatedCenters <- getMidPoints(cpm)

npClustering <- npMSL(
  data,
  nonDuplicatedCenters,
  samebw = FALSE ,
  verb = TRUE,
  maxiter = totalIterations,
  blockid= conds
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
write.csv(cluster, file = paste("results/clusters_21_",name,sep ="_"), row.names = row.names(data))
write.csv(posteriors, file = paste("results/posteriors_21_",name,sep ="_"), row.names = row.names(data))


