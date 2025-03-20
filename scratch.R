
# Set for debugging (if needed): Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
# ... in external terminal: cd to working directory, run "scratch.R"

rm(list = ls())
Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
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
  bs_chunksize <- 10
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
  LROcutoff = 2.0,                                      # cutoff for LROcp
  LROwindow_factor = 2.0,                               # window factor for LROcp, larger means larger rolling window
  LROfilter_ws_divisor = 2.0,                           # divisor for filter window size in likelihood ratio outlier detection (bigger is smaller window)
  tslope_initial = 1.0,                                 # initial value for tslope
  wf_initial = 0.15,                                    # initial value for wfactor
  max_evals = 1000,                                     # maximum number of evaluations for optimization
  initial_fits = 10,                                    # number of initial fits to perform in search of best initial conditions and best fit
  rng_seed = 42                                         # random seed for optimization (controls bootstrap resamples only)
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

# merfish.laminar.model <- wisp(
#   count.data.raw = countdata,
#   variables = data.variables,
#   use.median = FALSE,
#   bootstraps.num = 0,
#   converged.resamples.only = FALSE,
#   max.fork = bs_chunksize,
#   batch.size = bs_chunksize,
#   dim.bounds = layer.boundary.bins,
#   verbose = TRUE,
#   print.child.summaries = TRUE, 
#   model.settings = model.settings
# )


scan_LRO <- function() {
  
  ct <- c(1.5, 2.0, 2.5, 3.0)
  wf <- c(2.0, 3.0, 4.0)
  wd <- c(2.0, 3.0, 4.0)
  grd <- expand.grid(ct, wf, wd)
  
  time <- Sys.time()
  
  cps <- list()
  res <- list()
  bs <- list()
  for (r in 1:nrow(grd)) {
    
    model.settings = list(
      struc_values = c(                                     # values of structural parameters to test
        5.0,   # beta_shape_point
        5.0,   # beta_shape_rate
        1.0    # sd_tslope_effect
      ),  
      buffer_factor = 0.05,                                 # buffer factor for penalizing distance from structural parameter values
      ctol = 1e-6,                                          # convergence tolerance
      max_penalty_at_distance_factor = 0.01,                # maximum penalty at distance from structural parameter values
      LROcutoff = grd[r,1],                                      # cutoff for LROcp
      LROwindow_factor = grd[r,2],                               # window factor for LROcp, larger means larger rolling window
      LROfilter_ws_divisor = grd[r,3],                           # divisor for filter window size in likelihood ratio outlier detection (bigger is smaller window)
      tslope_initial = 1.0,                                 # initial value for tslope
      wf_initial = 0.15,                                    # initial value for wfactor
      max_evals = 1000,                                     # maximum number of evaluations for optimization
      initial_fits = 10,                                    # number of initial fits to perform in search of best initial conditions and best fit
      rng_seed = 42                                         # random seed for optimization (controls bootstrap resamples only)
    )
    
    m <- wisp(
      count.data.raw = countdata,
      variables = data.variables,
      use.median = FALSE,
      bootstraps.num = 0,
      converged.resamples.only = FALSE,
      max.fork = bs_chunksize,
      batch.size = bs_chunksize,
      dim.bounds = layer.boundary.bins,
      verbose = FALSE,
      print.child.summaries = FALSE, 
      model.settings = model.settings
    )
    
    cps[[r]] <- m$change.points 
    res[[r]] <- m$stats$residuals.log
    bs[[r]] <- m$bs.diagnostics[nrow(m$bs.diagnostics),]
    
    d <- Sys.time() - time
    cat("Finished", r, "of", nrow(grd), "in", d, units(d), "\n")
    time <- Sys.time()
    
  }
  
  return(list(cps = cps, res = res, bs = bs, grd = grd))
  
}

LROscan_res <- scan_LRO()



