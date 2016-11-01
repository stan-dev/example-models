rm(list = ls())
gc()

modelName <- "TwoCptModel"
# modelName <- "GenTwoCptModel"

## Adjust directories to your settings.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path("deliv", "figure", modelName)
tabDir <- file.path("deliv", "table", modelName)
modelDir <- file.path(projectDir, modelName)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path("tools")
stanDir <- file.path("cmdstan")

# source(file.path(scriptDir, "pkgSetup.R"))
library(rstan)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)

source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "cmdStanTools.R"))

rstan_options(auto_write = TRUE)

set.seed(11191951) ## not required but assures repeatable results

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CL", "Q", "V1", "V2", "ka", "sigma")

## Additional variables to monitor
otherRVs <- c("cObsPred")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

################################################################################################
## run Stan
nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

chains <- 1:nChains
mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init)
           runModel(model = model, data = data,
                    iter = iter, warmup = warmup, thin = thin,
                    init = init, seed = sample(1:999999, 1),
                    chain = chain),
         ##                      adapt_delta = 0.95, stepsize = 0.01),
         model = file.path(modelDir, modelName),
         data = file.path(modelDir, paste0(modelName, ".data.R")),
         init = file.path(modelDir, paste0(modelName, ".init.R")),
         iter = nIter, warmup = nBurn, thin = nThin,
         mc.cores = min(nChains, detectCores()))

fit <- read_stan_csv(file.path(modelDir, modelName, paste0(modelName, chains, ".csv")))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
mcmcDensity(fit, parametersToPlot)
pairs(fit, pars = parametersToPlot)

ptable <- parameterTable(fit, parametersToPlot)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
## posterior predictive distributions
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))
data <- data.frame(data$cObs, data$time[-1])
data <- plyr::rename(data, c("data.cObs" = "cObs", "data.time..1." = "time"))

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)

dev.off()
