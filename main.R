
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
bs_chunksize <- 2

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
  warp_precision = 1e-7                                 # pseudo infinity value larger than any possible possible parameter value, representing unbound warping
)

# Fit model
merfish.laminar.model <- wisp(
  # Data to model
  count.data.raw = count.data.WSPmm,
  # Variable labels
  variables = data.variables,
  # Local settings for specific fits, used on R side
  use.median = FALSE,
  MCMC.burnin = 0,
  MCMC.steps = 1e4,
  MCMC.step.size = 0.5,
  MCMC.prior = 10.0, 
  bootstraps.num = 0,
  converged.resamples.only = FALSE,
  max.fork = bs_chunksize,
  null.rate = log(2),
  null.slope = 1,
  dim.bounds = colMeans(layer.boundary.bins),
  verbose = TRUE,
  print.child.summaries = TRUE,
  # Global settings for initializing model, passed to C++ side
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










demo_warp <- function(
    w = 2, # warping factor
    point_pos = 60,
    point_neg = 40
  ) {
    
    # Data
    x <- (1:1000)/10
    b <- 100
    y <- x
    y1 <- WSP.warp(x, b, w)
    y2 <- WSP.warp(x, b, -w)
    y_pos <- WSP.warp(point_pos, b, w)
    y_neg <- WSP.warp(point_neg, b, -w)
    
    # Organize into a data frame
    df <- data.frame(
      x = rep(x, 3),
      y = c(y, y1, y2),
      curve = factor(rep(c("w = 0", "w > 0", "w < 0"), each = length(x)))
    )
    df$curve <- relevel(df$curve, ref = "w = 0")
    
    df_segments <- data.frame(
      point_pos = point_pos,
      point_neg = point_neg,
      y_pos = y_pos,
      y_neg = y_neg
    )
    
    # Make the ggplot
    demo_plot <- ggplot(df, aes(x = x, y = y, color = curve)) +
      geom_line(linewidth = 1.5) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "darkgray", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", linewidth = 1) +
      geom_segment(
        data = df_segments,
        aes(x = point_pos, xend = point_pos, y = point_pos, yend = y_pos),
        color = "blue4", linetype = "dashed", linewidth = 0.75) + 
      geom_segment(
        data = df_segments,
        aes(x = point_neg, xend = point_neg, y = point_neg, yend = y_neg),
        color = "red4", linetype = "dashed", linewidth = 0.75) +
      annotate("text", x = 10, y = 95, label = "upper asymptote", size = 7.5, color = "black") +
      annotate("text", x = 90, y = 5, label = "lower asymptote", size = 7.5, color = "black") +
      annotate("text", x = point_pos - 10, y = (y_pos + point_pos)/2, label = expression(varphi * "(z)(b - z)"), size = 7.5, color = "black") +
      annotate("text", x = point_neg + 10, y = (y_neg + point_neg)/2, label = expression(varphi * "(b - z)z"), size = 7.5, color = "black") +
      labs(
        x = "z",
        y = expression(omega * "(z, " * rho * ", b = 100)"),,
        title = "WSP Warping Function",
        color = "Direction"
      ) +
      scale_color_manual(
        values = c("black", "red", "blue"),
        labels = c(expression(rho * " = 0"), expression(rho * " < 0"), expression(rho * " > 0"))
        ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)
      )
    
    # Save at 1156 x 843
    return(demo_plot)
    
  }
demo_warp_plot <- demo_warp() 
print(demo_warp_plot)


