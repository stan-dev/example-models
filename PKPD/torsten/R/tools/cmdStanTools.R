## 5/27/2016: v1.0

## functions to run cmdStan

compileModel <- function(model, stanDir = stanDir){
  modelName <- basename(model)
  dir.create(model)
  file.copy(paste(model, "stan", sep = "."), file.path(model, paste(modelName, "stan", sep = ".")),
            overwrite = TRUE)
  model <- file.path(model, modelName)
  system(paste("make --directory=", stanDir, " ", model, sep = ""))
}

runModel <- function(model, data, iter, warmup, thin, init, seed, chain = 1,
                     stepsize = 1, adapt_delta = 0.8, max_depth = 10, refresh = 100, tag=NULL){
  modelName <- basename(model)
  model <- file.path(model, modelName)
  if(! is.null(tag)) output <- paste0(model, "_", tag, "_") else output=model
  system(paste(model, " sample algorithm=hmc engine=nuts",
               " max_depth=", max_depth,
               " stepsize=", stepsize,
               " num_samples=", iter,
               " num_warmup=", warmup, " thin=",  thin,
               " adapt delta=", adapt_delta, 
               " data file=", data,
               " init=", init, " random seed=", seed,
               " output file=", paste(output, chain, ".csv", sep = ""),
               " refresh=", refresh,
               sep = ""))
}

runDiagnose <- function(model, data, init, seed, chain = 1, refresh=100){
  modelName <- basename(model)
  model <- file.path(model, modelName)
  system(paste(model, " diagnose",
               " data file=", data,
               " init=", init, " random seed=", seed, 
               " output file=", paste(model, chain, ".csv", sep = ""),
               " refresh=", refresh,
               sep = ""))
}

# runModelFixed <- function(model, data, iter, warmup, thin, init, seed, chain = 1,
#                           stepsize = 1, adapt_delta = 0.8, max_depth = 10, refresh = 100){
#   modelName <- basename(model)
#   model <- file.path(model, modelName)
#   system(paste(model, " sample algorithm=fixed_param",
#                " num_samples=", iter,
#                " data file=", data,
#                " random seed=", seed,
#                " output file=", paste(model, chain, ".csv", sep = ""),
#                " refresh=", refresh,
#                sep = ""), invisible = FALSE)
# }

# runModelFixed <- function(model, data, iter, warmup, thin, init, seed, chain = 1,
#                           stepsize = 1, adapt_delt = 0.8, max_depth = 10, refresh = 100){
#   modelName <- basename(model)
#   model <- file.path(model, modelName)
#   print(paste0(model, " sample algorithm=fixed_param",
#                " num_samples=1 num_warmup=0",
#                " data file=", data,
#                " random seed=", seed,
#                " output file=", paste(model, chain, ".csv", sep = ""),
#                " refresh=", refresh))
#   
#   system(paste0(model, " sample algorithm=fixed_param",
#                " num_samples=1 num_warmup=0",
#                " data file=", data,
#                " init=", init,
#                " random seed=", seed,
#                " output file=", paste(model, chain, ".csv", sep = ""),
#                " refresh=", refresh), invisible = FALSE)
#   
# }

runModelFixed <- function(model, data, iter, warmup, thin, init, seed, chain = 1,
                          stepsize = 1, adapt_delta = 0.8, max_depth = 10, refresh = 100){
  modelName <- basename(model)
  model <- file.path(model, modelName)
  system(paste(model, " sample algorithm=fixed_param",
               " num_samples=", iter,
               " data file=", data,
               " random seed=", seed,
               " output file=", paste(model, chain, ".csv", sep = ""),
               " refresh=", refresh,
               sep = ""), invisible = FALSE)
}

