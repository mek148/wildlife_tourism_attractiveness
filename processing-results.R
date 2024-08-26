#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c~~~~~~#
## PROCESSING BIRD RESULTS ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Read in data from 500 models
filenames <- list.files(pattern = ".rds")
zapMCMClist <- lapply(filenames, "readRDS")

# extract fixed effects
zapMCMCsol <- lapply(zapMCMClist, function(m) m$Sol)
zapMCMCsol <- lapply(zapMCMCsol, function(m) m[,1:40]) ## remove random effect means
zapMCMCsol <- do.call(mcmc.list, zapMCMCsol)

# extract variances
zapMCMCvcv <- lapply(zapMCMClist, function(m) m$VCV)
zapMCMCvcv <- lapply(zapMCMCvcv, function(m) m[,1:5]) ## remove residual variance component of the binary model, which we fixed to 1
zapMCMCvcv <- do.call(mcmc.list, zapMCMCvcv)

# check convergence of chains
diag.sol <- gelman.diag(zapMCMCsol)  
diag.vcv <- gelman.diag(zapMCMCvcv)  
## the closer this factor is to 1, the better the convergence of our chains
## in practice, values below 1.1 can be acceptable and values below 1.02 are good

# check liability for binary process (second part of matrix)
zapMCMCliab <- lapply(zapMCMClist, function(m) m$Liab)
zapMCMCliab <- do.call(mcmc.list, zapMCMCliab) ## convert to single chain

ul <- 25
ll <- -25

liab <- lapply(1:500, function (i) { sum(zapMCMCliab[[i]][,(ncol(zapMCMCliab[[i]])/2+1):ncol(zapMCMCliab[[i]])] > ll & zapMCMCliab[[i]][,(ncol(zapMCMCliab[[i]])/2+1):ncol(zapMCMCliab[[i]])] < ul)/length(zapMCMCliab[[i]][,(ncol(zapMCMCliab[[i]])/2+1):ncol(zapMCMCliab[[i]])]) })
## 100% of liability estimates fall between -25 and 25, suggesting no numerical problems (see Jarrod's email)

# plot traces of fixed effects and variances
plotSOL <- plot(zapMCMCsol, ask = F) 
plotVCV <- plot(zapMCMCvcv, ask = F) 
dev.off()

# check effective sample sizes are adequate (see above)
effSOL <- effectiveSize(zapMCMCsol) # for the fixed effects 
effVCV <- effectiveSize(zapMCMCvcv[,1:4]) # for the phylo random effects

## most effective sample sizes are close to 1000, or big enough that >1000 could be achieved with a longer run
## be careful of warning message: "some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!" - a quick google suggests I need to reduce the complexity of the model in terms of the number of fixed effects (see here: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q2/023500.html). 

# test for autocorrelation
## values below 0.1 
autocorr.VCV <- autocorr.diag(zapMCMCvcv)
autocorr.Sol <- autocorr.diag(zapMCMCsol)

# function to calculate heritability (phylogenetic signal) for the ZAP model
herit_P <- function(x){
  (x$VCV[,1]/((x$VCV[,1]+x$VCV[,3])))
} ## for the Poisson part

herit_Z <- function(x){
  (x$VCV[,2]/((x$VCV[,2]+x$VCV[,4])+1))
} ## for the zero inflation part

# report mean heritability for each part of the model
her_P <- mean(as.numeric((lapply(1:2, function(i) mean(herit_P(zapMCMClist[[i]])))))) 
her_Z <- mean(as.numeric((lapply(1:2, function(i) mean(herit_Z(zapMCMClist[[i]]))))))

## these are both fairly high, which I think is encouraging as we should expect a phylogenetic signal in the data

# model summary
summVCV <- summary(zapMCMCvcv) ## variance
summSOL <- summary(zapMCMCsol) ## fixed effects

library(parallel)
# model predictions
predictions <- mclapply(1:500, function (i) { predict(zapMCMClist[[i]], marginal = NULL, 
                   type = "response", posterior = "mean") }, mc.cores = 10)

# fix political constant at maximum
# first read in bird datasets using example imputed dataset
species.df <- read.csv(paste0(getwd(),"/data/imputed_dataset1.csv"))

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

new.df <- species.df
new.df$scaled_ps <- max(newData$scaled_ps)*rnorm(nrow(new.df), 1, 1e-5)
## newData$scaled_range <- mean(newData$scaled_range)*rnorm(nrow(newData), 1, 1e-5)
## newData$scaled_dist <- min(newData$scaled_dist)*rnorm(nrow(newData), 1, 1e-5)

fixed.predictions <- mclapply(1:500, function (i) { predict(zapMCMClist[[i]], marginal = NULL, 
                    type = "response", posterior = "mean", newdata = new.df) }, mc.cores = 10)
# create dataframes of predictions
predictions.df <- as.data.frame(t(do.call(rbind.data.frame, Map('c', predictions))))
fixed.df <- as.data.frame(t(do.call(rbind.data.frame, Map('c', fixed.predictions))))

# calculate means and add to bird dataframe
species.df$pred <- predictions
species.df$pred_fixed <- fixed.predictions

# calculate residual values
species.df$resid <- species.df$pred - species.df$source
species.df$resid_fixed <- species.df$pred_fixed - species.df$source

# psuedo R2 using pearsons correlation
cor(species.df$pred, species.df$source, method = "pearson")^2 # 0.64
cor(species.df$pred_fixed, species.df$source, method = "pearson")^2 ## 0.58

# percentage of overlooked species
table(bspecies.df$resid >= 0.5)
(3606/9968)*100 ## 36.18%
table(species.df$resid_fixed >= 0.5)
(6099/9968)*100 ## 61.19%

## check predictions of zeros
## calculate number of zeroes in the actual data
## oz <- sum(species.df$source == 0)

## simulate 1000 times
## simZA <- mclapply(1:500, function (i) {

##  simulate(zapMCMClist[[i]], type = "response", posterior = "mean", nsim = 1000, newdata = species.df)

##  }, 

## mc.cores = 24)

## add up the zeroes for each simulation
## distZerosZA <- apply(simZA, 2, function (x) sum(x == 0))
## histogram of number of zeroes
## hist(distZerosZA)
## add line to show 'true' number
## abline(v = oz, col = "red")

# merge with common names
names <- read.csv("Birdlife_names.csv")
merged <- merge(species.df, names)

# export predictions
write.csv(merged, "new-predictions.csv", row.names = F)

# save image
rm(zapMCMClist,predictions.df,predictions,fixed.predictions,fixed.df)
save.image("MCMCbirds.RData")
