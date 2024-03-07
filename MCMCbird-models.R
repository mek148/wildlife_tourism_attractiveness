#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c~~~~~~#
## MCMC HURDLE MODEL ON NUMBER OF GUIDES MENTIONING BIRDS ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Libraries
library(ape)
library(MCMCglmm)
library(caper)
library(geiger)
library(parallel)

rm(list = ls())

######################### PREPARE DATA #########################

# Command argument for dataset and phylogenetic tree, as you are looping through 5 datasets and 100 trees
args <- commandArgs(T)
species.df <- read.csv(paste0(getwd(),"/data/",args[1],".csv"))
tree <- read.tree(paste0(getwd(),"/data/sample_100.tre"))[[as.numeric(args[2])]]

# Convert extinction risk and migration to continuous variable
species.df$ext_risk <- ifelse(species.df$status == "LC", 1, ifelse(species.df$status == "NT", 2, 
                     ifelse(species.df$status == "VU", 3, ifelse(species.df$status == "EN", 4, 5))))
species.df$migration <- ifelse(species.df$migration == "Not a Migrant", 1, 
                             ifelse(species.df$migration == "Nomadic", 2, 
                                    ifelse(species.df$migration == "Altitudinal Migrant", 3, 4)))

# Scale covariates
species.df[c("scaled_ED", "scaled_mass", "scaled_range", "scaled_ext", "scaled_ps", "scaled_dist", "scaled_migration", "scaled_elab", "scaled_div")] <- scale(species.df[c("ED_score", "logMass", "logRange", "ext_risk", "ps", "min_dist", "migration", "colour_elab", "colour_div")], center = T, scale = T)

# Order levels 
species.df$tl <- factor(species.df$tl, levels = c("Herbivore", "Omnivore", "Carnivore"))
species.df$tp <- factor(species.df$tp, levels = c("Nocturnal", "Diurnal"))
species.df$colony <- factor(species.df$colony, levels = c("Not colonial", "Colonial"))
species.df$habitat1 <- factor(species.df$habitat1, levels = c("Open vegetation", "Bare", "Forest", "Shrubland", "Artificial", "Aquatic", "Other"))

######################### RUN MODEL #########################

# Use scientific name as row names
rownames(species.df) <- species.df$scientific_name

# Check names match using function from {geiger} and remove non-matching species
nameCheck <- name.check(tree, species.df)
pruned.tree <- drop.tip(tree, nameCheck$tree_not_data) ## Drops tips from multiple trees that are not present in dataset

# Rename first column 
colnames(species.df)[1] <- "species"

# Calculate the inverse of the matrices of phylogenetic correlations 
inv.phylo <- inverseA(pruned.tree, nodes = "ALL", scale = T)

# Poisson model exploring factors that influence number of mentions in guidebooks
prior <- list(R = list(V = diag(2), nu = 0.002, fix = 2), 
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2), 
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2),
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

# Run multiple chains using 100 different trees and five imputed datasets to obtain 500 models incorporating phylogenetic uncertainty and uncertainty in the imputation procedure. Zero-altered Poisson to account for large number of zeros. 
set.seed(1) 

mod <- MCMCglmm(
      source ~ trait *
        (
          scaled_mass + scaled_ED + scaled_ps +  scaled_range * scaled_dist +
            scaled_range * scaled_ext + scaled_elab + scaled_div + tl + tp + colony +
            habitat1 + scaled_migration
        ),
      random = ~ idh(at.level(trait, 1)):species + idh(at.level(trait, 2)):species +
        idh(at.level(trait, 1)):Order + idh(at.level(trait, 2)):Order,
      ginverse = list(species = inv.phylo$Ainv),
      rcov = ~ idh(trait):units,
      data = species.df,
      prior = prior,
      family = "zapoisson",
      nitt = 3500000, 
      thin = 3000, 
      burnin = 10000, 
      verbose = F,
      pl = T,
      pr = T
    )

# Save models
saveRDS(mod, paste0("MCMCbird-models-",args[1],"-",args[2],".rds"))

