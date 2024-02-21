#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Testing sbatch script on multiple datasets and multiple trees ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Index for dataset and phylogenetic tree
args <- commandArgs(T)

print(args[1])

df <- read.csv(paste0(args[1],".csv"))

head(df)