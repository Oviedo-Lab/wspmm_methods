
# Set for debugging (if needed): Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
# ... in external terminal: cd to working directory, run "scratch.R"

rm(list = ls())

projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"
sink("output.txt", split = TRUE, append = FALSE, type = "output")

# Set random seed for reproducibility
# ... R only. C++ seed set in its code
ran.seed <- 12349999
set.seed(ran.seed)

# For for WSPmm and snk printing
library(wispack)

# Set bootstrap chunk size
sys_name <- Sys.info()["sysname"]
if (sys_name == "Darwin" || sys_name == "Linux") {
  bs_chunksize <- 2
} else {
  bs_chunksize <- 0
}

# Load demo data
countdata <- read.csv("countdata.csv")

# Note on spatial coordinates: 
#  data from horizontal (axial) slice of mouse cortex, ACx
#  y-axes: laminar, 0 is bottom of L6b, 100 is top of L2/3
#  x-axes: columnar, 0 is most posterior, 100 is most anterior

# Define list of genes to analyze
gene.list <- c("Bcl11b", "Fezf2", "Rorb", "Satb2", "Nxph3", "Cux2", "Rorb") 

# Define cortical layer boundaries
# ... Numbers represent the lower boundary of each layer, in bins
# ... Higher numbers closer to cortical surface
layer.boundary.bins <- c(
  NA,    # Layer 1 (removed from data)
  70.25, # Layer 2/3
  57.00, # Layer 4
  21.50, # Layer 5
  1.00,  # Layer 6a
  0.00   # Layer 6b
)

# Define fixed effects to test
fixed.effect.names <- c("hemisphere", "age")

# Define variables in the dataframe for the model
data.variables = list(
  count = "count",
  bin = "bin", 
  parent = "cortex", 
  child = "gene",
  ran = "mouse",
  fixedeffects = fixed.effect.names
)

# Model settings 
model.settings = list(
  # ... these are global options needed to set up model
  buffer_factor = 0.05,                                 # buffer factor for penalizing distance from structural parameter values
  ctol = 1e-6,                                          # convergence tolerance
  max_penalty_at_distance_factor = 0.01,                # maximum penalty at distance from structural parameter values
  LROcutoff = 2.0,                                      # cutoff for LROcp
  LROwindow_factor = 2.0,                               # window factor for LROcp, larger means larger rolling window
  LROfilter_ws_divisor = 2.0,                           # divisor for filter window size in likelihood ratio outlier detection (bigger is smaller window)
  rise_threshold_factor = 0.8,                          # amount of detected rise as fraction of total required to end run
  max_evals = 1000,                                     # maximum number of evaluations for optimization
  rng_seed = 42,                                        # random seed for optimization (controls bootstrap resamples only)
  inf_warp = 1e3                                        # pseudo infinity value larger than any possible possible parameter value, representing unbound warping
)

# Setting suggestions: 
# - Recommend turning settings on data for which you have strong priors. 
# - If bootstraps not fitting due to falling off boundary (very high boundary penalty and few iterations), 
#    increase max_penalty_at_distance_factor from 0.01 to 0.1 or 1.0, so the gradient descent algorithm has 
#    more information about the boundary edge. 
# - If model is overly biased to finding sharp, small change points and missing obvious smooth, large ones, 
#    try increasing the LROwindow_factor. A larger window will be more tuned to large but gradual transitions,
#    while a smaller window will be more tuned to sharp transitions (whether large or small). 
# - Adjusting the LROcutoff is another way of controlling how the model finds rate transitions. Higher values 
#    mean fewer detected transitions. 

merfish.laminar.model <- wisp(
  # Data to model
  count.data.raw = countdata,
  # Variable labels
  variables = data.variables,
  # Local settings for specific fits, used on R side
  use.median = FALSE,
  MCMC.burnin = 0,
  MCMC.steps = 1e3,
  MCMC.step.size = 0.005,
  MCMC.prior = 0.5,                                     
  bootstraps.num = 1e2,
  converged.resamples.only = FALSE,
  max.fork = bs_chunksize,
  dim.bounds = layer.boundary.bins,
  verbose = TRUE,
  print.child.summaries = TRUE,
  # Global settings for initializing model, passed to C++ side
  model.settings = model.settings
)



newplots <- plot.ratecount(
  merfish.laminar.model,
  pred.type = "pred",
  count.type = "count",
  print.all = TRUE
)

plot.MCMC.walks(merfish.laminar.model)








bs_fitted_params <- merfish.laminar.model$sample.params
tslope_effs_mask <- grepl("tslope", merfish.laminar.model$param.names) & grepl("baseline", merfish.laminar.model$param.names)
tslope_effects <- c(bs_fitted_params[,tslope_effs_mask])
hist(tslope_effects)
mean(tslope_effects)
sd(tslope_effects)


