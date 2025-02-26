
# Setup ################################################################################################################
# Analysis of MERFISH data with Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)

# Clear global environment
rm(list = ls())
# If using VS Code with httpgd, ensure clean start
httpgd::hgd_close() 
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

# For for WSPmm 
wispack_path <- paste0(projects_folder, "R_packages/wispack/wispack_1.0.tar.gz")
install.packages(wispack_path, repos = NULL)
library(wispack)

# Set random seed for reproducibility
# ... R only. C++ seed set in its code
ran.seed <- 12349999
set.seed(ran.seed)

# Set for debugging (if needed): Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
# ... in external terminal: cd to working directory, run "Rscript script__main_merfish_analysis.R"

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

# Transform coordinates for each mouse and extract layer boundary estimates
count_data <- transform_and_extract(
  count_data, 
  100,        # Number of bins to use when binning data
  TRUE        # Keep coordinate transformation plots? 
)

# Unpack 
layer_boundary_bins <- count_data$layer_boundary_bins
coordinate_transform_plots <- count_data$plots
count_data <- count_data$df

# Fit WSPmm model to MERFISH data ######################################################################################

# Define list of genes to analyze
gene_list_plasticity <- c(
  "Dlx2", "Grin2a", "Arc", "Nptxr", "Camk2a", "Nr4a2", "Ncdn",
  "Bcl11b", "Fezf2", "Rorb", "Satb2"
)

# Define fixed effects to test
fixed.effect.names <- c("hemisphere", "age")

# Data variables 
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
  struc_values = c(5.0, 5.0, 5.0, 1.0, 1.0, 1.0, 1.0),  # values of structural parameters to test
  buffer_factor = 0.05,                                 # buffer factor for penalizing distance from structural parameter values
  ctol = 1e-6,                                          # convergence tolerance
  max_penalty_at_distance_factor = 0.01,                # maximum penalty at distance from structural parameter values
  LROcutoff = 2.0,                                      # cutoff for LROcp
  tslope_initial = 1.0,                                 # initial value for tslope
  wf_initial = 0.05,                                     # initial value for wfactor
  max_evals = 500                                       # maximum number of evaluations for optimization
)

# Create count data for WSPmm object, from preprocessed count_data, using laminar axis (y)
count_data_WSPmm.y <- create.count.data.WSPmm(
  df.merfish = count_data,
  bin.dim = "y_bins",
  gene.list = gene_list_plasticity,
  fixed.effect.names = c("hemisphere","age"),
  verbose = TRUE
)

merfish_laminar_model <- wisp(
  count.data.raw = count_data_WSPmm.y,
  variables = data.variables,
  bootstraps.num = 1e3,
  converged.resamples.only = FALSE,
  max.fork = bs_chunksize,
  batch.size = bs_chunksize,
  dim.bounds = colMeans(layer_boundary_bins),
  verbose = TRUE,
  print.child.summaries = TRUE, 
  model.settings = model.settings
)

