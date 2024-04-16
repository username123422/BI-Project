# BI-Project
## In this project we aim to use RNA abundance (from RNA sequencing data), translation (measured using ribosome profiling), and mRNA sequence to predict protein abundace. We will use the data from the different experiments for 1 cell line (the HeLa cell line). 
## Load Data



## Normalize  RNA abundance (Measured from RNA Seq Data) using count per million reads
##This normalizes the entire RNA Abundance data set

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

## Normalize Translation (Measured using ribosome profiling) using count per million reads

ribo_cl_train_data <- ribo_cl_train
str(ribo_cl_train)

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

counts_per_million(ribo_cl_train)

normalized_counts_ribo = counts_per_million(ribo_cl_train[,-1]) 

##test : We did this test to verify that we did our normilization correctly by just doing the normilization on 1 experiment individually and comparing the values to those in the normalized_counts_rna
test_df = rna_cl_train[,1:5]
test_df$nor_GSM546922 = test_df$GSM546922 / sum(test_df$GSM546922) * 1000000

## Loading Experiment data

experiment_data <- read.csv("human_infor_train.csv")

str(experiment_data)

## Subset Data 
##We want to extract the experiments from the HeLa Cell line from each data set

##Datasets

##normalized_counts_rna
##normalized_counts_ribo
##experiment_data

##Trying to create a dataset with just the experiments from the HeLa cell line (from the data set with the experiments and cell lines)

heLa_experiments <- unique(experiment_data$experiment_alias[experiment_data$cell_line == "HeLa"])

