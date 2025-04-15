
# Setup ################################################################################################################
# Analysis of MERFISH data with Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)

# Clear global environment
rm(list = ls())
# If using VS Code with httpgd, ensure clean start
httpgd::hgd_close() 
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

# For for WSPmm and snk printing
library(wispack)

# Set random seed for reproducibility
# ... R only. C++ seed set in its code
ran.seed <- 12349999
set.seed(ran.seed)
options(error = recover)
# Set for debugging (if needed): Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
# ... in external terminal: cd to working directory, run "Rscript main.R"

# Report what's happening and send all output to both the console and text file
sink("output.txt", split = TRUE, append = FALSE, type = "output")
snk.report("Analysis of MERFISH data by Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)\n")

# Set file paths and bootstrap chunk size
sys_name <- Sys.info()["sysname"]
if (sys_name == "Darwin") {
  source("merfish_preprocessing.R")
  data_path <- paste0(projects_folder, "MERFISH/data/")
  bs_chunksize <- 10
} else if (sys_name == "Linux") {
  source("~/MERFISH/MERFISH-ACx-Spatial-Density/merfish_preprocessing.R")
  data_path <- "/mnt/c/OviedoLab/MERFISH/MERFISH_Data/All_HDF5/"
  bs_chunksize <- 25
}

# Preprocessing MERFISH data ###########################################################################################

# Load and parse raw visgen data from files
count_data <- make_count_data(data_path)

# Transform coordinates for each mouse into laminar and columnar axes and extract layer boundary estimates
# Note: 
#  laminar axis goes to y, with cortical surface in the positive direction (zero is deepest layer, L6)
#  columnar axis goes to x, with anterior in the positive direction (zero is most posterior)
count_data <- cortical_coordinate_transform(
  count_data = count_data, 
  total_bins = 100,        # Number of bins to use when binning data
  keep_plots = TRUE,       # Keep coordinate transformation plots? 
  verbose = TRUE
)

# Unpack 
layer.boundary.bins <- count_data$layer.boundary.bins
coordinate_transform_plots <- count_data$plots
count_data <- count_data$df

# Check orientation of columnar axes in bin coordinates 
columnar_check <- function() {
  
  # Reshape data
  rcdX <- count_data[, c("mouse", "hemisphere", "x_coord", "x_bins")]
  colnames(rcdX) <- c("mouse", "hemisphere", "cen", "xbin")
  rcdY <- count_data[, c("mouse", "hemisphere", "y_coord", "x_bins")]
  colnames(rcdY) <- c("mouse", "hemisphere", "cen", "xbin")
  axis <- rep("x", nrow(rcdX))
  rcdX <- cbind(rcdX, axis)
  axis <- rep("y", nrow(rcdY))
  rcdY <- cbind(rcdY, axis)
  rcd <- rbind(rcdX, rcdY)
  HxA <- paste0(rcd$hemisphere, "_", rcd$axis)
  HxAxM <- paste0(rcd$hemisphere, "_", rcd$axis, "_", rcd$mouse)
  rcd <- cbind(rcd, HxA, HxAxM)
  
  # Make plot
  transform_plot <- ggplot(rcd) +
    geom_point(aes(x = cen, y = xbin, color = factor(mouse))) +
    facet_wrap(~HxAxM, scales = "free_x") +
    labs(title = "Coordinate transforms") +
    theme_minimal()
  plot(transform_plot)
  
}
columnar_check()

# Fit WSPmm model to MERFISH data ######################################################################################

# Define list of genes to analyze
gene.list <- c("Rorb", "Pvalb", "Gad2", "Vip", "Grik3", "Grm1", "Slc32a1", "Sox6", "Reln", "Npsr1", "Dscaml1", "Calb1") 

# Create count data for WSPmm object, from preprocessed count_data, using laminar axis (y)
count.data.WSPmm.y <- create.count.data.WSPmm(
  df.merfish = count_data,
  bin.dim = "y_bins",
  gene.list = gene.list,
  fixed.effect.names = c("hemisphere","age")
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
  rng_seed = 42                                         # random seed for optimization (controls bootstrap resamples only)
)

merfish.laminar.model <- wisp(
  # Data to model
  count.data.raw = count.data.WSPmm.y,
  # Variable labels
  variables = data.variables,
  # Local settings for specific fits, used on R side
  use.median = FALSE,
  MCMC.burnin = 0,
  MCMC.steps = 1e4,
  MCMC.step.size = 0.005,
  MCMC.prior = 0.5,                                     
  bootstraps.num = 0,
  converged.resamples.only = FALSE,
  max.fork = bs_chunksize,
  dim.bounds = layer.boundary.bins,
  verbose = TRUE,
  print.child.summaries = TRUE,
  # Global settings for initializing model, passed to C++ side
  model.settings = model.settings
)





