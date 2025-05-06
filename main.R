
# Setup ################################################################################################################
# Analysis of MERFISH data with Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)

# Clear global environment
rm(list = ls())
# If using VS Code with httpgd, ensure clean start
httpgd::hgd_close() 
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

# For for WSPmm and snk ("sink") printing
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

# Make results table for stats
param_stats <- laminar.model[["stats"]][["parameters"]][,-c(5,7)]
param_stats[,2:5] <- round(param_stats[,2:5], 4)
raneff_mask <- grepl("wfactor", param_stats$parameter)
param_stats_ran <- param_stats[raneff_mask,]
param_stats_fix <- param_stats[!raneff_mask,]
param_stats_list <- list() 
for (g in gene.list) {
  # Grab gene and baseline masks
  gene_mask_ran <- grepl(g, param_stats_ran$parameter)
  gene_mask_fix <- grepl(g, param_stats_fix$parameter)
  baseline_mask <- grepl("baseline", param_stats_fix$parameter)
  # Split into lists by gene and parameter type (ran, fix, baseline)
  param_stats_ran_g <- param_stats_ran[gene_mask_ran,]
  param_stats_fix_g <- param_stats_fix[gene_mask_fix & !baseline_mask,]
  param_stats_baseline_g <- param_stats_fix[gene_mask_fix & baseline_mask,]
  # Split parameter names into columns to reorganize data
  split_cols_ran <- do.call(rbind, strsplit(param_stats_ran_g$parameter, "_")) # want cols 2,3 (spatial param, ran level)
  split_cols_fix <- do.call(rbind, strsplit(param_stats_fix_g$parameter, "_")) # want cols 2,5,7 (spatial param, treatment, block)
  split_cols_baseline <- do.call(rbind, strsplit(param_stats_baseline_g$parameter, "_")) # want cols 3,5 (spatial param, block)
  # Reorganize data
  results_cols <- c(2:6)
  col_names <- c("est", "CI.low", "CI.high", "p.adj", "sig")
  param_stats_ran_list_g <- list()
  param_stats_fix_list_g <- list()
  param_stats_baseline_list_g <- list()
  for (ranlvl in unique(split_cols_ran[,3])) {
    ranlvl_mask <- split_cols_ran[,3] == ranlvl
    results <- param_stats_ran_g[ranlvl_mask, results_cols]
    param_type <- split_cols_ran[ranlvl_mask,2]
    level <- rep(paste0("mouse ", ranlvl), length(results[,1]))
    treatment <- rep("random", length(results[,1]))
    param_stats_ran_list_g[[ranlvl]] <- as.data.frame(cbind(treatment, param_type, level, results))
    colnames(param_stats_ran_list_g[[ranlvl]]) <- c("effect", "spatial.param", "index", col_names)
  }
  param_stats_ran_list_g <- as.data.frame(do.call(rbind, param_stats_ran_list_g))
  for (trt in unique(split_cols_fix[,5])) {
    treatment_mask <- split_cols_fix[,5] == trt
    results <- param_stats_fix_g[treatment_mask, results_cols]
    param_type <- split_cols_fix[treatment_mask,2]
    param_type[param_type == "Rt"] <- "rate"
    param_type[param_type == "tpoint"] <- "position"
    param_type[param_type == "tslope"] <- "slope scalar"
    block <- split_cols_fix[treatment_mask,7]
    block <- gsub("Tns/Blk", "", block)
    block[param_type == "rate"] <- paste0("block ", block[param_type == "rate"])
    block[param_type != "rate"] <- paste0("t-point ", block[param_type != "rate"])
    treatment <- rep(trt, length(results[,1]))
    param_stats_fix_list_g[[trt]] <- as.data.frame(cbind(treatment, param_type, block, results))
    colnames(param_stats_fix_list_g[[trt]]) <- c("effect", "spatial.param", "index", col_names)
  }
  param_stats_fix_list_g <- as.data.frame(do.call(rbind, param_stats_fix_list_g))
  for (sp in unique(split_cols_baseline[,3])) {
    sp_mask <- split_cols_baseline[,3] == sp
    results <- param_stats_baseline_g[sp_mask, results_cols]
    if (sp == "Rt") {
      param_type <- rep("rate", length(results[,1]))
    } else if (sp == "tpoint") {
      param_type <- rep("position", length(results[,1]))
    } else {
      param_type <- rep("slope scalar", length(results[,1]))
    }
    treatment <- rep("baseline", length(results[,1]))
    block <- split_cols_baseline[sp_mask,5]
    block <- gsub("Tns/Blk", "", block)
    if (sp == "Rt") {
      block <- paste0("block ", block)
    } else {
      block <- paste0("t-point ", block)
    }
    param_stats_baseline_list_g[[sp]] <- as.data.frame(cbind(treatment, param_type, block, results))
    colnames(param_stats_baseline_list_g[[sp]]) <- c("effect", "spatial.param", "index", col_names)
  }
  param_stats_baseline_list_g <- as.data.frame(do.call(rbind, param_stats_baseline_list_g))
  param_stats_list[[g]] <- as.data.frame(do.call(rbind, list(param_stats_baseline_list_g, param_stats_fix_list_g, param_stats_ran_list_g)))
  # For each gene and parameter type (ran, fix, baseline), will have five numeric value columns: est, low, high, p-adj, significance
}

library(dplyr)
library(knitr)
library(kableExtra)

# Sample data
df <- param_stats_list[["Rorb"]]
rownames(df) <- NULL
write.csv(df, file = "test_stats.csv", row.names = FALSE)






df <- read.csv("test_stats.csv")
effect_lengths <- c()
for (e in unique(df$effect)) {
  effect_lengths <- c(effect_lengths, sum(df$effect == e))
}
names(effect_lengths) <- c(
  "baseline (P12, left)", 
  "Fixed Effects: Hemisphere (right)", 
  "Fixed Effect: Age (P18)", 
  "Fixed Effect: Hemisphere-Age Interaction", 
  "Random Effects"
)
df <- df[,-which(colnames(df) == "effect")]

# Create the table with row grouping
kbl(df, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "") %>%
  group_rows(index = effect_lengths) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"), font_size = 8)





# Replace repeated Group values with LaTeX \multirow
df2 <- df %>%
  group_by(effect) %>%
  mutate(effect = ifelse(row_number() == 1, 
                        paste0("\\multirow{", n(), "}{*}{", effect, "}"), 
                        "")) %>%
  ungroup()

# Render table
kbl(df2, format = "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(latex_options = "hold_position")



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

# OTHER DEMOS!! (PUT INTO DEMO SCRIPT)
demo_sigmoid()
demo_warp_plot <- demo_warp() 
print(demo_warp_plot)


new_plot <- laminar.model[["plots"]][["ratecount"]][["plot_pred_parent_cortex_fixEff_Rorb"]] +
  theme(
    plot.title = element_text(hjust = 0.5, size = 30),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
print(new_plot)

new_plot2 <- new_plot + scale_color_manual(values = c("deeppink", "deeppink4", "deepskyblue", "deepskyblue4"))
print(new_plot2)



# Method for printing all child plots on one figure
plot.child.summary2 <- function(
    wisp.results,
    these.parents = NULL,
    these.childs = NULL,
    verbose = TRUE
) {
  
  if (verbose) {
    snk.report("Printing child summary plots", initial_breaks = 2)
    snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 0)
  }
  
  # Check for needed plots
  if (
    length(wisp.results$plots$ratecount) == 0 || 
    length(wisp.results$plots$parameters) == 0 ||
    length(wisp.results$plots$residuals) == 0
  ) {
    stop("No rate or parameter plots found in wisp.results")
  }
  
  # Specify parent and child levels to summarize
  gvp_lvls <- as.character(wisp.results$grouping.variables$parent.lvls)
  if (length(these.parents) != 0) gvp_lvls <- gvp_lvls[gvp_lvls %in% these.parents]
  gv_lvls <- as.character(wisp.results$grouping.variables$child.lvls)
  if (length(these.childs) != 0) gv_lvls <- gv_lvls[gv_lvls %in% these.childs]
  
  for (gvp_lvl in gvp_lvls) {
    
    if (verbose) snk.report(paste0("Making summary plots for ", gvp_lvl))
    first_print <- TRUE
    
    for (gv_lvl in gv_lvls) {
      
      if (verbose && first_print) {
        cat(gv_lvl)
        first_print <- FALSE
      } else if (verbose) {
        cat(",", gv_lvl)
      }
      wisp.results$plots$parameters <- plot.parameters(
        wisp.results = wisp.results,
        child.lvl = gv_lvl, 
        print.plots = FALSE, 
        verbose = FALSE 
      )
      
      # Grab / make plots
      p1 <- wisp.results$plots$ratecount
      p2 <- wisp.results$plots$parameters
      p3 <- wisp.results$plots$residuals
      
      # Find plots for this parent
      p_mask1 <- grepl(gvp_lvl, names(p1))
      p_mask2 <- grepl(gvp_lvl, names(p2))
      p_mask3 <- grepl(gvp_lvl, names(p3))
      
      # Find plots for this child
      c_mask1 <- grepl(gv_lvl, names(p1))
      c_mask2 <- grepl(gv_lvl, names(p2)) 
      c_mask3 <- grepl(gv_lvl, names(p3))
      
      # Find residual hist and qq plots 
      histqq <- grepl("hist|qq", names(p3))
      
      # Find treatment plots
      iX_mask <- grepl("treatment", names(p2))
      
      # Subset
      p_rates <- p1[p_mask1 & c_mask1]
      p_treatment <- p2[p_mask2 & c_mask2 & iX_mask]
      p_otherparams <- p2[p_mask2 & c_mask2 & !iX_mask]
      p_residuals <- p3[p_mask3 & c_mask3 & histqq]
      
      for (pi in 1:length(p_treatment)) {
        p_treatment[[pi]] <- p_treatment[[pi]]  +
          theme(
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 9)
          )
      }
      for (pi in 1:length(p_otherparams)) {
        p_otherparams[[pi]] <- p_otherparams[[pi]]  +
          theme(
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 10),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 9)
          )
      }
      
      # Combine and print
      resids <- do.call(arrangeGrob, c(p_residuals, ncol = length(p_residuals)))
      rates_and_residuals <- arrangeGrob(ggplotGrob(p_rates[[1]]), resids, ncol = 1)
      treatments <- do.call(arrangeGrob, c(as.list(p_treatment), ncol = length(p_treatment)))
      other_params <- do.call(arrangeGrob, c(as.list(p_otherparams), ncol = length(p_otherparams)))
      params <- arrangeGrob(treatments, other_params, ncol = 1)
      grid.arrange(params)
      return(params)
      p <- arrangeGrob(rates_and_residuals, params, ncol = 2, widths = c(0.4,0.6))
      grid.arrange(p)
      
    }
    
  }
  
}
my_plots <- plot.child.summary2(
  laminar.model,
  these.childs = "Rorb",
  verbose = TRUE
)
