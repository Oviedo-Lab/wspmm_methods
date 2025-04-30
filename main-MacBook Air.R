
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
ran.seed <- 123
set.seed(ran.seed)
options(error = recover)
# Set for debugging (if needed): Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
# ... in external terminal: cd to working directory, run "Rscript main.R"

# Report what's happening and send all output to both the console and text file
sink("output.txt", split = TRUE, append = FALSE, type = "output")
snk.report("Analysis of MERFISH data by Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)\n")

# Set file paths and bootstrap chunk size
source("merfish_preprocessing.R")
data_path <- paste0(projects_folder, "MERFISH/data_SSp/")
bs_chunksize <- 10

# Define list of genes to analyze
gene.list <- c("Bcl11b", "Fezf2", "Satb2", "Nxph3", "Cux2", "Rorb")  

# Preprocessing MERFISH data ###########################################################################################

# Load and parse raw visgen data from files
count_data <- make_count_data(
  data_path,
  remove_L1 = TRUE,
  ROIname = "Primary somatosensory area",
  raw = FALSE
  )

# Transform coordinates for each mouse into laminar and columnar axes and extract layer boundary estimates
# Note: 
#  laminar axis goes to y, with cortical surface in the positive direction (zero is deepest layer, L6)
#  columnar axis goes to x, with anterior in the positive direction (zero is most posterior)
count_data <- cortical_coordinate_transform(
  count_data = count_data, 
  total_bins = 100,        # Number of bins to use when binning data
  keep_plots = TRUE,       # Keep coordinate transformation plots? 
  nat_right = TRUE,
  verbose = TRUE
)

# Unpack 
layer.boundary.bins <- count_data$layer.boundary.bins
coordinate_transform_plots <- count_data$plots
count_data <- count_data$df

# Save layer boundaries
write.csv(
  layer.boundary.bins,
  file = "layer_boundary_bins.csv",
  row.names = FALSE
)

# Fit WSPmm model to MERFISH data ######################################################################################

# Define fixed effects to test
fixed.effect.names <- c("hemisphere", "age")

# Create count data for WSPmm object, from preprocessed count_data, using laminar axis (y)

count.data.WSPmm <- create.count.data.WSPmm(
  df.merfish = count_data,
  bin.dim = "y_bins",
  gene.list = gene.list,
  fixed.effect.names = fixed.effect.names
)

write.csv(
  count.data.WSPmm, 
  file = "S1_laminar_countdata.csv",
  row.names = FALSE
)

#count.data.WSPmm <- read.csv("S1_laminar_countdata.csv")

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
  LROcutoff = 1.5,                                      # cutoff for LROcp
  LROwindow_factor = 3.0,                               # window factor for LROcp, larger means larger rolling window
  rise_threshold_factor = 0.8,                          # amount of detected rise as fraction of total required to end run in initial slope estimation
  max_evals = 1000,                                     # maximum number of evaluations for optimization
  rng_seed = 42,                                        # random seed for optimization (controls bootstrap resamples only)
  warp_precision = 1e-7                                 # decimal precision to retain when selecting really big number as pseudo infinity for unbound warping
)

# Fit model
laminar.model <- wisp(
  # Data to model
  count.data.raw = count.data.WSPmm,
  # Variable labels
  variables = data.variables,
  # Settings used on R side
  use.median = FALSE,
  MCMC.burnin = 0,
  MCMC.steps = 1e4,
  MCMC.step.size = 0.5,
  MCMC.prior = 10.0, 
  bootstraps.num = 1e4,
  converged.resamples.only = TRUE,
  max.fork = bs_chunksize,
  null.rate = log(2),
  null.slope = 1,
  dim.bounds = colMeans(layer.boundary.bins),
  verbose = TRUE,
  print.child.summaries = TRUE,
  # Setting to pass to C++ model
  model.settings = model.settings
)



# Scratch #####

newplots <- plot.ratecount(
  merfish.laminar.model,
  pred.type = "pred.log",
  count.type = "count.log",
  print.all = TRUE
)

View(merfish.laminar.model[["MCMC.diagnostics"]][c(1:10,9990:10001),])

nll <- unlist(merfish.laminar.model[["MCMC.diagnostics"]][["neg.loglik"]])
pnll <- unlist(merfish.laminar.model[["MCMC.diagnostics"]][["pen.neg.value"]])
nll <- (nll - nll[1]) / nll[1]
pnll <- (pnll - pnll[1]) / pnll[1]
ymin <- min(c(nll, pnll))
ymax <- max(c(nll, pnll))
plot(nll, type = "l", col = "blue", ylim = c(ymin, ymax))
lines(pnll, col = "red")

nll <- unlist(merfish.laminar.model[["MCMC.diagnostics"]][["neg.loglik"]])
pnll <- unlist(merfish.laminar.model[["MCMC.diagnostics"]][["pen.neg.value"]])
idx <- which(merfish.laminar.model$param.names == "beta_Rt_cortex_Rorb_right18_X_Tns/Blk2")
ptrace <- merfish.laminar.model[["sample.params"]][,idx]
ptracemax <- max(ptrace)
ptracemin <- min(ptrace)
nll <- scales::rescale(nll, to = c(ptracemin, ptracemax), from = range(nll))
pnll <- scales::rescale(pnll, to = c(ptracemin, ptracemax), from = range(pnll))
plot(ptrace, type = "l")
lines(pnll, col = "red")
lines(nll, col = "blue")



# Compare bs and MCMC results
# ... nll
nll_mcmc <- laminar.model[["diagnostics.MCMC"]]$neg.loglik
nll_bs <- laminar.model[["diagnostics.bs"]]$neg.loglik
bs_success <- laminar.model[["diagnostics.bs"]]$success.code
bs_mask <- bs_success == 3
bs_mean <- mean(nll_bs[bs_mask])
nll_bs[!bs_mask] <- bs_mean
plot(nll_bs, type = "l", col = "blue", ylim = range(c(nll_mcmc, nll_bs)))
lines(nll_mcmc, col = "red")
# ... param
param_mcmc <- laminar.model[["sample.params.MCMC"]]
param_bs <- laminar.model[["sample.params.bs"]]
n_params <- ncol(param_mcmc)
n_samples <- nrow(param_mcmc)
ttest_results <- rep(NA, n_params)
for (i in 1:n_params) {
  redraw <- 10
  downsample <- 40
  for (j in 1:redraw) {
    this_draw <- sample(1:n_samples, downsample, replace = FALSE)
    ttest_results[i] <- t.test(param_mcmc[this_draw,i], param_bs[this_draw,i])$p.value
  }
  ttest_results[i] <- ttest_results[i] / redraw
}
hist(ttest_results)
plot(ttest_results)

library(ggplot2)

comp_den <- function(this_param) {
  # Put samples into a data frame
  df <- data.frame(
    value = c(param_mcmc[,this_param],param_bs[,this_param]),
    group = factor(rep(c("mcmc", "bs"), each = n_samples))
  )
  
  # Plot
  ggplot(df, aes(x = value, color = group)) +
    geom_density(linewidth = 1.2) +
    theme_minimal() +
    labs(title = "Density Comparison", x = "Value", y = "Density")
  
}
comp_den(1)
for (i in 11:40) {
  print(comp_den(i))
}




tt <- sample(100:200, 500, replace = TRUE)
test <- density(tt)
test

comp_den2 <- function() {
  # Precompute density estimates for each param and group
  dens_list <- lapply(1:ncol(param_mcmc), function(pnum) {
    mcmc_vals <- param_mcmc[, pnum]
    bs_vals   <- param_bs[, pnum]
    
    # Density and centering
    d_mcmc <- density(mcmc_vals)
    d_bs   <- density(bs_vals)
    mcmc_peak <- d_mcmc$x[which.max(d_mcmc$y)]
    bs_peak   <- d_bs$x[which.max(d_bs$y)]
    
    # Return centered density curves
    data.frame(
      x = c(abs(d_mcmc$x - mcmc_peak)+1, abs(d_bs$x - bs_peak)+1),
      y = c(d_mcmc$y+1, d_bs$y+1),
      group = rep(c("mcmc", "bs"), each = length(d_mcmc$x)),
      param = factor(pnum)
    )
  })
  
  df_dens <- do.call(rbind, dens_list)
  
  # Plot precomputed densities
  ggplot(df_dens, aes(x = x, y = y, color = group, group = interaction(group, param))) +
    geom_line(linewidth = 1, alpha = 0.25) +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal() +
    labs(title = "Centered Density Comparison",
         x = "Centered Value (peak = 0)", y = "Density")
}

comp_den2()



comp_den2 <- function() {
  # Precompute density estimates for each param and group
  dens_list <- lapply(1:ncol(param_mcmc), function(pnum) {
    mcmc_vals <- param_mcmc[, pnum]
    bs_vals   <- param_bs[, pnum]
    n_resamples <- 100
    mcmc_sw_stat <- rep(NA, n_resamples)
    bs_sw_stat   <- rep(NA, n_resamples)
    
    sample_size <- 35
    for (i in 1:n_resamples) {
      mcmc_vals_sample <- sample(mcmc_vals, sample_size, replace = TRUE)
      bs_vals_sample   <- sample(bs_vals, sample_size, replace = TRUE)
      mcmc_sw_stat[i] <- shapiro.test(mcmc_vals_sample)$statistic
      bs_sw_stat[i]   <- shapiro.test(bs_vals_sample)$statistic
    }
    
    # Put samples into a data frame
    df <- data.frame(
      value = c(mcmc_sw_stat, bs_sw_stat),
      group = factor(rep(c("mcmc", "bs"), each = n_resamples))
    )
    
    # Plot
    ggplot(df, aes(x = value, color = group)) +
      geom_density(linewidth = 1.2) +
      theme_minimal() +
      labs(title = "Shaprio-Wilk Normality Test Statistic on Resampled Parameters", x = "Value", y = "Density")
  })
  
  df_dens <- do.call(rbind, dens_list)
  
  # Plot precomputed densities
  ggplot(df_dens, aes(x = x, y = y, color = group, group = interaction(group, param))) +
    geom_line(linewidth = 1, alpha = 0.25) +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal() +
    labs(title = "Centered Density Comparison",
         x = "Centered Value (peak = 0)", y = "Density")
}








demo_sigmoid()

demo_warp_plot <- demo_warp() 
print(demo_warp_plot)


