#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c~~~~~~#
## MCMC HURDLE MODEL ON NUMBER OF GUIDES MENTIONING BIRDS ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Libraries
library(ape)
library(MCMCglmm)
library(caper)
library(geiger)

rm(list = ls())

# Load data
load(paste0(getwd(),"/input data/",args[1],"-",args[2]))

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
      nitt = 3000000, 
      thin = 3000, 
      burnin = 10000, 
      verbose = F,
      pl = T,
      pr = T
    )

# Save models
saveRDS(mod, paste0("MCMCbird-models-",args[1],"-",args[2],".rds"))

