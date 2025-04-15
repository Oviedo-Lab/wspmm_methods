
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

# Define list of genes to analyze
gene.list <- c("Rorb", "Pvalb", "Gad2", "Vip", "Grik3", "Grm1", "Slc32a1", "Sox6", "Reln", "Npsr1", "Dscaml1", "Calb1") 

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

# Preprocessing MERFISH data, Allen ####################################################################################

make_ROImasks <- FALSE 
process_slices <- FALSE

# Source code / load data
if (make_ROImasks || process_slices) {
  if (sys_name == "Darwin") {
    source("merfish_preprocessing_Allen.R")
  } else if (sys_name == "Linux") {
    source("~/MERFISH/MERFISH-ACx-Spatial-Density/merfish_preprocessing_Allen.R")
  }
} else {
  if (sys_name == "Darwin") {
    S1_allen_slice_data_annotated <- read.csv("S1_allen_slice_data_annotated.csv")
  } else if (sys_name == "Linux") {
    S1_allen_slice_data_annotated <- read.csv("~/MERFISH/MERFISH-ACx-Spatial-Density/S1_allen_slice_data_annotated.csv")
  }
}

# Define helper functions to format registered Allen data for coordinate transformation

reformat_good_allen_data <- function(
    S1_allen_slice_data_annotated
  ) {
    
    # Rename columns for coordinate transform code
    new_names <- colnames(S1_allen_slice_data_annotated)
    new_names[c(3)] <- "mouse" # ... for now, treat layers as mice
    new_names[c(4)] <- "hemisphere"
    new_names[c(6, 7)] <- c("x_coord", "y_coord")
    colnames(S1_allen_slice_data_annotated) <- new_names
    # ... renumber the mouse column 
    S1_allen_slice_data_annotated$mouse <- as.integer(as.factor(S1_allen_slice_data_annotated$mouse))
    # ... from visual inspection, slices 6, 5, and 4 are the cleanest, so collapse together into one mouse. 
    good_layer_mask <- S1_allen_slice_data_annotated$mouse == 6 | S1_allen_slice_data_annotated$mouse == 5 | S1_allen_slice_data_annotated$mouse == 4
    S1_allen_slice_data_annotated <- S1_allen_slice_data_annotated[good_layer_mask,]
    S1_allen_slice_data_annotated$mouse <- 1
    return(S1_allen_slice_data_annotated)
    
  }

make_allen_slice_plots <- function(
    S1_allen_slice_data_annotated
  ) {
    
    slice_plots <- list() 
    for (m in unique(S1_allen_slice_data_annotated$mouse)) {
      slice_plots[[paste0("slice_plot", m)]] <- plot_results(
        S1_allen_slice_data_annotated, 
        S1_allen_slice_data_annotated[S1_allen_slice_data_annotated$mouse == m, c("x_coord", "y_coord")], 
        m,
        paste0("Registered S1 Allen slice data, slice ", m),
        separate_hemi = TRUE)
    }
    return(slice_plots)
    
  }

# Put data into form expected for coordinate transformations
S1_allen_slice_data_annotated <- reformat_good_allen_data(S1_allen_slice_data_annotated)
count_data_allen <- list(
  count_data = S1_allen_slice_data_annotated,
  slice_plots = make_allen_slice_plots(S1_allen_slice_data_annotated)
)

# Perform coordinate transformations and binnings
count_data_allen <- cortical_coordinate_transform(
  count_data = count_data_allen, 
  total_bins = 100,        # Number of bins to use when binning data
  keep_plots = TRUE,       # Keep coordinate transformation plots? 
  L1_removed = FALSE, 
  nat_left = TRUE, 
  verbose = TRUE
)

# Cleanup 
rm(S1_allen_slice_data_annotated)

# Fit WSPmm model to MERFISH data ######################################################################################

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





