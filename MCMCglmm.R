#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## MCMC HURDLE MODEL ON NUMBER OF GUIDES MENTIONING BIRDS ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Libraries
library(ape)
library(MCMCglmm)
library(caper)
library(geiger)
library(parallel)

rm(list = ls())

######################### PREPARE DATA #########################

# Imputed bird dataset 
birdData <- read.csv(paste0(getwd(),"/data/imputed_dataset 1 .csv"))

# Bird pylogenies
sampleTrees <- read.tree(paste0(getwd(),"/data/sample_100.tre"))

# Convert extinction risk and migration to continuous variable
birdData$ext_risk <- ifelse(birdData$status == "LC", 1, ifelse(birdData$status == "NT", 2, 
                     ifelse(birdData$status == "VU", 3, ifelse(birdData$status == "EN", 4, 5))))
birdData$migration <- ifelse(birdData$migration == "Not a Migrant", 1, 
                             ifelse(birdData$migration == "Nomadic", 2, 
                                    ifelse(birdData$migration == "Altitudinal Migrant", 3, 4)))

# Scale covariates
birdData[c("scaled_ED", "scaled_mass", "scaled_range", "scaled_ext", "scaled_ps", "scaled_dist", "scaled_migration", "scaled_elab", "scaled_div")] <- scale(birdData[c("ED_score", "logMass", "logRange", "ext_risk", "ps", "min_dist", "migration", "colour_elab", "colour_div")], center = T, scale = T)

# Convert binary variables back to factors
birdData$colony <- ifelse(birdData$colony == 1, "Colonial", "Not colonial")
birdData$tp <- ifelse(birdData$tp == 1, "Diurnal", "Nocturnal")

# Order levels 
birdData$tl <- factor(birdData$tl, levels = c("Herbivore", "Omnivore", "Carnivore"))
birdData$tp <- factor(birdData$tp, levels = c("Nocturnal", "Diurnal"))
birdData$colony <- factor(birdData$colony, levels = c("Not colonial", "Colonial"))
birdData$habitat1 <- factor(birdData$habitat1, levels = c("Open vegetation", "Bare", "Forest", "Shrubland", "Artificial", "Aquatic", "Other"))

######################### RUN MODEL #########################

# Use scientific name as row names
rownames(birdData) <- birdData$scientific_name

# Check names match using function from {geiger} and remove non-matching species
nameCheck <- name.check(sampleTrees[[1]], birdData)
prunedTrees <- lapply(sampleTrees, drop.tip, nameCheck$tree_not_data) ## Drops tips from multiple trees that are not present in dataset

# Rename first column 
colnames(birdData)[1] <- "species"

# Calculate the inverse of the matrices of phylogenetic correlations 
invPhylo <- lapply(1:100, function (l) inverseA(prunedTrees[[l]], nodes = "ALL", scale = T))

# Poisson model exploring factors that influence number of mentions in guidebooks
prior <- list(R = list(V = diag(2), nu = 0.002, fix = 2), 
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2), 
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2),
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

# Run multiple chains using 100 different trees and five imputed datasets to obtain 500 models incorporating phylogenetic uncertainty and uncertainty in the imputation procedure. Zero-altered Poisson to account for large number of zeros. 
set.seed(1) 

zapMCMClist <- mclapply(1:5, function (j) {
  mclapply(1:100, function(i) {
    MCMCglmm(
      source ~ trait *
        (
          scaled_mass + scaled_ED + scaled_ps +  scaled_range * scaled_dist +
            scaled_range * scaled_ext + scaled_elab + scaled_div + tl + tp + colony +
            habitat1 + scaled_migration
        ),
      random = ~ idh(at.level(trait, 1)):species + idh(at.level(trait, 2)):species +
        idh(at.level(trait, 1)):Order + idh(at.level(trait, 2)):Order,
      ginverse = list(species = invPhylo[[i]]$Ainv),
      rcov = ~ idh(trait):units,
      data = birdData[[j]],
      prior = prior,
      family = "zapoisson",
      nitt = 3500,
      thin = 30,
      burnin = 100,
      verbose = F,
      pl = T,
      pr = T
    )
    
  },
  
  mc.cores = 100)
  
})

# Save models
save.image("MCMCbird-models.RData")

####################### MODEL CHECKING #########################

# Extract fixed effects
zapMCMCsol <- lapply(zapMCMClist, function(m) m$Sol)
zapMCMCsol <- lapply(zapMCMCsol, function(m) m[,1:42]) ## remove random effect means
zapMCMCsol <- do.call(mcmc.list, zapMCMCsol)

# Extract variances
zapMCMCvcv <- lapply(zapMCMClist, function(m) m$VCV)
zapMCMCvcv <- lapply(zapMCMCvcv, function(m) m[,1:5]) ## Remove residual variance component of the binary model, which we fixed to 1
zapMCMCvcv <- do.call(mcmc.list, zapMCMCvcv)

# Check convergence of chains
diag.sol <- gelman.diag(zapMCMCsol)  
diag.vcv <- gelman.diag(zapMCMCvcv)  
## The closer this factor is to 1, the better the convergence of our chains. In practice, values below 1.1 can be acceptable and values below 1.02 are good

# Check liability for binary process (second part of matrix)
zapMCMCliab <- lapply(zapMCMClist, function(m) m$Liab)
zapMCMCliab <- do.call(mcmc.list, zapMCMCliab) # Cnvert to single chain

ul <- 25
ll <- -25

liab <- lapply(1, function (i) { sum(zapMCMCliab[[i]][,(ncol(zapMCMCliab[[i]])/2+1):ncol(zapMCMCliab[[i]])] > ll & zapMCMCliab[[i]][,(ncol(zapMCMCliab[[i]])/2+1):ncol(zapMCMCliab[[i]])] < ul)/length(zapMCMCliab[[i]][,(ncol(zapMCMCliab[[i]])/2+1):ncol(zapMCMCliab[[i]])]) })
## 100% of liability estimates fall between -25 and 25, suggesting no numerical problems (see Jarrod's email)

# Plot traces of fixed effects and variances
plotSOL <- plot(zapMCMCsol, ask = F) 
plotVCV <- plot(zapMCMCvcv, ask = F) 
dev.off()

# Check effective sample sizes are adequate (see above)
effSOL <- effectiveSize(zapMCMCsol) # Fixed effects 
effVCV <- effectiveSize(zapMCMCvcv[,1:4]) # Phylo random effects
## Most effective sample sizes are close to 1000, or big enough that >1000 could be achieved with a longer run

# Test for autocorrelation. Values below 0.1 are good.
autocorr.VCV <- autocorr.diag(zapMCMCvcv)
autocorr.Sol <- autocorr.diag(zapMCMCsol)

# Function to calculate heritability (phylogenetic signal)
herit_P <- function(x){
  (x$VCV[,1]/((x$VCV[,1]+x$VCV[,3])))
} ## Poisson part

herit_Z <- function(x){
  (x$VCV[,2]/((x$VCV[,2]+x$VCV[,4])+1))
} ## Zero inflation part

# Report mean heritability for each part of the model
her_P <- mean(as.numeric((lapply(1:2, function(i) mean(herit_P(zapMCMClist[[i]])))))) 
her_Z <- mean(as.numeric((lapply(1:2, function(i) mean(herit_Z(zapMCMClist[[i]]))))))
## These are both fairly high, which I think is encouraging as we should expect a phylogenetic signal in the data

# Model summary
summVCV <- summary(zapMCMCvcv) ## Variance
summSOL <- summary(zapMCMCsol) ## Fixed effects

######################## PREDICTIONS #########################

# Model predictions
predictions <- mclapply(1, function (x) { predict(zapMCMClist[[x]], marginal = NULL, 
          type = "response", posterior = "mean") }, mc.cores = 2)

# Fix political stability constant at maximum
newData <- birdData
newData$scaled_ps <- max(newData$scaled_ps)*rnorm(nrow(newData), 1, 1e-5)

fixed.predictions <- mclapply(1:500, function (i) { 
  predict(zapMCMClist[[i]], marginal = NULL, newdata = newData, 
          type = "response", posterior = "mean") }, mc.cores = 2) 

# Create data frames of predictions
predictionsDF <- as.data.frame(t(do.call(rbind.data.frame, Map('c', predictions))))
fixedDF <- as.data.frame(t(do.call(rbind.data.frame, Map('c', fixed.predictions))))

# Export
write.csv(predictionsDF, "bird-predictions.csv", row.names = F)
write.csv(fixedDF, "bird-fixed-predictions.csv", row.names = F)

# Save image without large list of models
rm(zapMCMClist)
save.image("MCMCbird-results.RData")