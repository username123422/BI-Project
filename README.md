# BI-Project
## In this project we aim to use RNA abundance (from RNA sequencing data), translation (measured using ribosome profiling), and mRNA sequence to predict protein abundace. We will use the data from the different experiments for 1 cell line (the HeLa cell line). 
## Load Data



##Normalize using count per million reads


rna_cl_train_data <- rna_cl_train
str(rna_cl_train)

colSums(rna_cl_train[,-1])

counts_per_million <- function (count_matrix) { 
  total_reads = colSums(count_matrix)  
  normalized_matrix = count_matrix
  for (gene in 1:nrow(count_matrix)) { 
    normalized_matrix [gene, ] = 1000000 * count_matrix[gene , ] / total_reads
  }
  return (normalized_matrix )
}
counts_per_million <- function (count_matrix) { 
  return (t(apply(count_matrix, 1, function(x){1000000*x/colSums(count_matrix)})))
}

counts_per_million(rna_cl_train)

normalized_counts = counts_per_million(rna_cl_train[,-1]) 

## Subset Data 
