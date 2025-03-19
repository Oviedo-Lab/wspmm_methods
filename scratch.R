

rm(list = ls())
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

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
  struc_values = c(                                     # values of structural parameters to test
    5.0,   # beta_shape_point
    5.0,   # beta_shape_rate
    1.0    # sd_tslope_effect
  ),  
  buffer_factor = 0.05,                                 # buffer factor for penalizing distance from structural parameter values
  ctol = 1e-6,                                          # convergence tolerance
  max_penalty_at_distance_factor = 0.01,                # maximum penalty at distance from structural parameter values
  LROcutoff = 3.0,                                      # cutoff for LROcp
  LROwindow_factor = 3.0,                               # window factor for LROcp
  tslope_initial = 1.0,                                 # initial value for tslope
  wf_initial = 0.15,                                    # initial value for wfactor
  max_evals = 1000,                                     # maximum number of evaluations for optimization
  rng_seed = 42                                         # random seed for optimization
)

# Setting suggestions: 
# - If bootstraps not fitting due to falling off boundary (very high boundary penalty and few iterations), 
#    increase max_penalty_at_distance_factor from 0.01 to 0.1 or 1.0, so the gradient descent algorithm has 
#    more information about the boundary edge. 
# - If model is overly biased to finding sharp, small change points and missing obvious smooth, large ones, 
#    try increasing the LROwindow_factor.

merfish.laminar.model <- wisp(
  count.data.raw = countdata,
  variables = data.variables,
  use.median = FALSE,
  bootstraps.num = 20,
  converged.resamples.only = FALSE,
  max.fork = bs_chunksize,
  batch.size = bs_chunksize,
  dim.bounds = layer.boundary.bins,
  verbose = TRUE,
  print.child.summaries = TRUE, 
  model.settings = model.settings
)


# Jitter (or even recalculate) the initial parameters at each bootstrap iteration 
# Some kind of recovery of the main fit so that it tries several different starting points 
#   and goes with the best? (Not necessarily only for "recovery" reasons.) 
# The diff ratio algorithm should run a bunch to generate an average with CI. 


