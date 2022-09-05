library(HTSCluster)
library(HTSFilter)
library(Biobase)

set.seed(123)
name = "data/raw_counts_mouse.csv"

#read the data and combine it. 
data <- read.csv(name, row.names = 1, header= TRUE)

max_index = 35
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

data_filter <- HTSFilter(data, conds, norm="TMM")

time_start <- Sys.time();
pmm <- PoisMixClusWrapper(y=data_filter$filteredData, gmin=1, gmax = max_index, conds=conds, split.init=TRUE)
time_end <- Sys.time();

for (index in 2:max_index){
  new_index = paste("g=",index, sep = "")
  write.csv(
    pmm$all.results[[new_index]]$labels,
    file = paste("results_poisson/poisson_cluster",index,"label.csv", sep = "_"),
    row.names = row.names(data_filter[["filteredData"]])
  )
  write.csv(
    pmm$all.results[[new_index]]$probaPost,
    file = paste("results_poisson/poisson_cluster",index,"posterior.csv", sep = "_"),
    row.names = row.names(data_filter[["filteredData"]])
  )
}
saveRDS(pmm, file = "results_poisson/poisson_info.RDS")


