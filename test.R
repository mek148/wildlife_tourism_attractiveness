#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Testing sbatch script on multiple datasets and multiple trees ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

args <- commandArgs(T)

print(args[1])

df <- read.csv(paste0(getwd(),"/data/",args[1],".csv"))

head(df)

write.csv(df, paste0(getwd(),"/data/new-",args[1],".csv"))