library(dplyr)

f_counts <- read.table("genes.counts.matrix",header = T, stringsAsFactors = F, sep = "\t")
colnames(f_counts)[1] <- "GeneID"
f_id <- read.table("gene_only_tcdb.txt",stringsAsFactors = F, sep = "\t", header = T)

f_trans_exp <- dplyr::inner_join(f_id, f_counts, by = "GeneID")
colnames(f_trans_exp)[1] <- ""
write.table(f_trans_exp, file = "genes.counts.matrix.tcdb",
            sep = "\t", row.names = F, quote = F)
