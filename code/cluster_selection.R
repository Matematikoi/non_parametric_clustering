library("mixtools")

set.seed(8)
name <- "data/1_5_cpm_mouse_log.csv"
max_selection <- 20
max_iteration <- 500

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

start_time <- Sys.time()
cluster_selection <- repnormmixmodel.sel( 
  t(data), 
  k = max_selection, 
  verb = TRUE, 
  epsilon = 1e-02,
  maxit = max_iteration
)
end_time <- Sys.time()

write.csv(cluster_selection, file = "cluster_selection/mouse_cluster_selection.csv")
