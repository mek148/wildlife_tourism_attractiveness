#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Testing sbatch script on multiple datasets and multiple trees ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

args <- commandArgs(T)

print(args[1])
print(args[2])

df <- read.csv(paste0(getwd(),"/data/",args[1],".csv"))

trees <- read.tree(paste0(getwd(),"/data/sample_100.tre"))
tree <- trees[args[2]]