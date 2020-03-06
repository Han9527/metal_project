metabolites <- read.table("./tcdb_heatmap.txt",
                            header = T,
                            stringsAsFactors = F,
                            row.names = 1)
mean_loop <- as.data.frame(rownames(metabolites))
colnames(mean_loop)[1] <- "GeneID"
v <- rownames(metabolites)
length(metabolites)
for (i in seq(1, length(metabolites), by = 3)) {
  v1 <- metabolites[i]
  v2 <- metabolites[i+1]
  v3 <- metabolites[i+2]
  d <- data.frame(v1, v2, v3)
  mean_m <- data.frame(apply(d, 1, mean)) 
  metabolites_mean <- cbind(v, mean_m)
  colnames(metabolites_mean)[1] <- "GeneID"
  colnames(metabolites_mean)[2] <- substring(colnames(metabolites)[i], 1, 2)
  mean_loop <- dplyr::inner_join(mean_loop,metabolites_mean,by="GeneID") 
  if (i+2 > 18)
    break
  print(i)
  write.table(mean_loop, 
              file = "tcdb_mean", 
              sep = "\t", 
              quote = FALSE,
              row.names = FALSE)
}


