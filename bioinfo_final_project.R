# BI-Project
## In this project we aim to use RNA abundance (from RNA sequencing data), translation (measured using ribosome profiling), and mRNA sequence to predict protein abundace. We will use the data from the different experiments for 1 cell line (the HeLa cell line). 
## Load Data



## Data Investigation and loading relevant functions for normalization

rna_cl_train_data <- rna_cl_train
str(rna_cl_train)


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


ribo_cl_train_data <- ribo_cl_train
str(ribo_cl_train)

## Subset Data 
##We want to extract the experiments from the HeLa Cell line from each data set
column_names <- c("GSM2100596","GSM2100597","GSM2100598",
                  "GSM2100599",
                  "GSM2100600",
                  "GSM2100601",
                  "GSM2100602",
                  "GSM2734576",
                  "GSM2734577",
                  "GSM2734578",
                  "GSM2734579",
                  "GSM2734580",
                  "GSM2734581",
                  "GSM2734582",
                  "GSM2734583",
                  "GSM2734584",
                  "GSM2734585",
                  "GSM2817679",
                  "GSM2817680",
                  "GSM2817681",
                  "GSM2817682",
                  "GSM2817683",
                  "GSM2817684",
                  "GSM2817685",
                  "GSM2817686",
                  "GSM2817687",
                  "GSM2817688",
                  "GSM2817689",
                  "GSM2825128",
                  "GSM2825129",
                  "GSM2825131",
                  "GSM2825132",
                  "GSM3098458",
                  "GSM3098459",
                  "GSM3098460",
                  "GSM3098461",
                  "GSM3098462",
                  "GSM3098463",
                  "GSM3098464",
                  "GSM3098465",
                  "GSM3098466",
                  "GSM3098467",
                  "GSM3098468",
                  "GSM3098469",
                  "GSM3098494",
                  "GSM3098495",
                  "GSM3098496",
                  "GSM3098497",
                  "GSM3098498",
                  "GSM3098499",
                  "GSM3098500",
                  "GSM3098501",
                  "GSM3098502",
                  "GSM3098503",
                  "GSM3098504",
                  "GSM3098505",
                  "GSM3168233",
                  "GSM3168235",
                  "GSM3168236",
                  "GSM3562582",
                  "GSM3562583",
                  "GSM3562584",
                  "GSM3562585",
                  "GSM3562586",
                  "GSM3562587",
                  "GSM3944606",
                  "GSM3944607",
                  "GSM3944608",
                  "GSM3944609",
                  "GSM3944610",
                  "GSM3944611",
                  "GSM3944612",
                  "GSM3944613",
                  "GSM3944614",
                  "GSM3944615",
                  "GSM3944616",
                  "GSM3944617",
                  "GSM3944618",
                  "GSM3944619",
                  "GSM3944620",
                  "GSM4154182",
                  "GSM4154183",
                  "GSM4154184",
                  "GSM5291033",
                  "GSM5291034",
                  "GSM5291035",
                  "GSM5291036",
                  "GSM5291037",
                  "GSM5291038",
                  "GSM546922",
                  "GSM546924",
                  "GSM546926",
                  "GSM546928",
                  "GSM546930",
                  "GSM546920")

experiment_rna_subset <- rna_cl_train_data[,column_names]
experiment_ribo_subset <- ribo_cl_train_data[,column_names]

#Now we normalize the relevant experimental data using counts per a million reads

#For RNA abundance
colSums(experiment_rna_subset[,-1])

counts_per_million(experiment_rna_subset)
normalized_rna_counts = counts_per_million(experiment_rna_subset[,-1])

#And translational data
colSums(experiment_ribo_subset[,-1])

counts_per_million(experiment_ribo_subset)
normalized_rna_counts = counts_per_million(experiment_ribo_subset[,-1])