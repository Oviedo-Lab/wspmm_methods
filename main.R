
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
snk.report("Analysis of MERFISH data by Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)", end_breaks = 1)

# Source preprocessing functions, set file path, and bootstrap chunk size
source("merfish_preprocessing.R")
data_path <- paste0(projects_folder, "MERFISH/data_SSp/")
bs_chunksize <- 10

# Define list of genes to analyze
gene.list <- c("Bcl11b", "Fezf2", "Satb2", "Nxph3", "Cux2", "Rorb")  

# Preprocessing MERFISH data ###########################################################################################

# Load and parse raw visgen data from files
count_data <- make_count_data(
    data_path,
    remove_L1 = TRUE,        # If TRUE, removes any points labeled as layer 1
    ROIname = "Primary somatosensory area",
    raw = TRUE               # If FALSE, uses normalized data, and that is not desired
  )

# Transform coordinates for each mouse into laminar and columnar axes and extract layer boundary estimates
# Note: 
#  laminar axis goes to y, with cortical surface in the positive direction (zero is deepest layer, L6)
#  columnar axis goes to x, with anterior in the positive direction (zero is most posterior)
count_data <- cortical_coordinate_transform(
    count_data = count_data, 
    total_bins = 100,        # Number of bins to use when binning data
    keep_plots = TRUE,       # Keep coordinate transformation plots? 
    nat_right = TRUE,        # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    verbose = TRUE
  )

# Unpack 
layer.boundary.bins <- count_data$layer.boundary.bins
coordinate_transform_plots <- count_data$plots
count_data <- count_data$df

# Simple check of transcripts per cell per gene per mouse
extract_transcript_counts <- function() {
    snk.report("Transcript counts per cell per gene per mouse")
    mice <- unique(count_data$mouse)
    counts <- count_data[,which(colnames(count_data) %in% gene.list)]
    counts_all <- count_data[,c(18:738)]
    for (i in mice) {
      snk.print_vec("mouse", i, initial_breaks = 2)
      mask <- count_data$mouse == i
      snk.print_vec("All genes", sum(counts_all[mask,])/sum(mask)/ncol(counts_all))
      snk.print_vec("modeled genes only", sum(counts[mask,])/sum(mask)/ncol(counts))
    }
  }
extract_transcript_counts()

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

# Save pre-processed count data
# ... load with: count.data.WSPmm <- read.csv("S1_laminar_countdata.csv")
write.csv(
    count.data.WSPmm, 
    file = "S1_laminar_countdata.csv",
    row.names = FALSE
  )

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
    LROcutoff = 2.0,                                      # cutoff for LROcp, a multiple of standard deviation
    LROwindow_factor = 1.25,                              # window factor for LROcp, larger means larger rolling window
    rise_threshold_factor = 0.8,                          # amount of detected rise as fraction of total required to end run in initial slope estimation
    max_evals = 1000,                                     # maximum number of evaluations for optimization
    rng_seed = 42,                                        # random seed for optimization (controls bootstrap resamples only)
    warp_precision = 1e-7                                 # decimal precision to retain when selecting really big number as pseudo infinity for unbound warping
  )

# Settings for MCMC walk
MCMC.settings = list(
    MCMC.burnin = 1e2,
    MCMC.steps = 1e4,
    MCMC.step.size = 1.0,
    MCMC.prior = 1.0, 
    MCMC.neighbor.filter = 2
  )

# Fit model
laminar.model <- wisp(
    # Data to model
    count.data = count.data.WSPmm,
    # Variable labels
    variables = data.variables,
    # Settings used on R side
    use.median = FALSE,
    MCMC.settings = MCMC.settings,
    bootstraps.num = 1e4,
    converged.resamples.only = TRUE,
    max.fork = bs_chunksize,
    dim.bounds = colMeans(layer.boundary.bins),
    verbose = TRUE,
    print.child.summaries = TRUE,
    # Setting to pass to C++ model
    model.settings = model.settings
  )

# Save
# ... load with: laminar.model <- readRDS("saved_laminar_model-final-noplots.rds")
saveRDS(laminar.model, file = "saved_laminar_model.rds")

# Make and export figures #####
# ... These functions are a mess, but they're one-offs for this data set and this paper

make_fig_results_ratecount <- function() {
    
    colors4 <- c("firebrick1", "firebrick4", "dodgerblue1", "dodgerblue4")
    decomposed_plots <- plot.decomposition(laminar.model, "Rorb")
    
    n_plots <- length(laminar.model[["plots"]][["ratecount"]])
    p_names <- names(laminar.model[["plots"]][["ratecount"]])
    p_names <- p_names[2:n_plots]
    ratecount_plots <- list()
    length(ratecount_plots) <- n_plots - 1
    names(ratecount_plots) <- p_names
    title_size <- 40 
    axis_size <- 24 
    legend_size <- 20
    ratecount_plots_Rorb_ran <- list()
    mice <- c("mouse 1", "mouse 2", "mouse 3", "mouse 4")
    length(ratecount_plots_Rorb_ran) <- 4
    names(ratecount_plots_Rorb_ran) <- mice
    for (p in 1:4) {
      
      gene_name <- bquote("ROR" * beta * ", " * .(mice[p]))
      ratecount_plots_Rorb_ran[[mice[p]]] <- plot.ratecount(
        wisp.results = laminar.model,
        pred.type = "pred",
        count.type = "count",
        dim.boundaries = unlist(laminar.model[["plots"]][["ratecount"]][["plot_pred_parent_cortex_fixEff_Bcl11b"]][["layers"]][[49]][["data"]]),
        #y.lim = c(0, 215),
        count.alpha.none = 0,
        count.alpha.ran = 0.5,
        pred.alpha.none = 0,
        pred.alpha.ran = 1,
        rans.to.print = as.character(p),
        childs.to.print = c("Rorb")
      )[[1]] +
        labs(title = gene_name, x = NULL, y = NULL) + 
        theme(
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size),
          legend.position = "none"
        ) + 
        scale_color_manual(
          labels = c("Left, P12", "right, P12", "left, P18", "right, P18"),
          values = colors4
        )
      
    }
    for (p in 2:n_plots) {
      
      p_ <- names(laminar.model[["plots"]][["ratecount"]])[p]
      gene_name <- gsub("plot_pred_parent_cortex_fixEff_", "", p_)
      
      # Check if this is Rorb
      this_Rorb <- FALSE
      if (gene_name == "Rorb") this_Rorb <- TRUE
      
      # Set legend position 
      leg_pos <- "none"
      if (this_Rorb) leg_pos = "bottom"
      
      # Reformat gene names
      if (this_Rorb) {
        gene_name <- expression("ROR" * beta)
      } else {
        gene_name <- toupper(gene_name)
      }
      
      # Remake Rorb 
      if (this_Rorb) {
        
        rorb_decomp<- plot.ratecount(
          wisp.results = laminar.model,
          pred.type = "pred",
          count.type = "count",
          dim.boundaries = unlist(laminar.model[["plots"]][["ratecount"]][["plot_pred_parent_cortex_fixEff_Bcl11b"]][["layers"]][[49]][["data"]]),
          count.alpha.none = 0.5,
          count.alpha.ran = 0,
          pred.alpha.none = 1,
          pred.alpha.ran = 0,
          rans.to.print = "none",
          childs.to.print = c("Rorb")
        )
        ratecount_plots[[p_]] <- rorb_decomp[[1]]
        
      } else {
        
        # Recolor plot
        ratecount_plots[[p_]] <- laminar.model[["plots"]][["ratecount"]][[p_]] 
        
      }
      
      # Recolor plot
      ratecount_plots[[p_]] <- ratecount_plots[[p_]] +
        labs(title = gene_name, x = NULL, y = NULL) + 
        theme(
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size),
          legend.position = leg_pos
        ) + 
        scale_color_manual(
          name = "",
          labels = c("Left, P12", "right, P12", "left, P18", "right, P18"),
          values = colors4
        )
      
      if (this_Rorb) {
        # Extract data frame from plot 
        found_P12 <- FALSE 
        this_layer <- 1
        while(!found_P12) {
          df12 <- ratecount_plots[[p_]][["layers"]][[this_layer]][["data"]]
          if(all(df12$ran == "none" & df12$treatment == "ref")) found_P12 <- TRUE
          else this_layer <- this_layer + 1
        }
        found_P18 <- FALSE 
        this_layer <- 1 
        while(!found_P18) {
          df18 <- ratecount_plots[[p_]][["layers"]][[this_layer]][["data"]]
          if(all(df18$ran == "none" & df18$treatment == "ref")) found_P18 <- TRUE
          else this_layer <- this_layer + 1
        }
        # Find t-points
        rise <- 25
        tpoints_ref <- laminar.model$fitted.parameters[grepl("baseline_cortex_tpoint_Rorb", laminar.model$param.names)]
        tpoints_18 <- laminar.model$fitted.parameters[grepl("beta_tpoint_cortex_Rorb_18_X_Tns/Blk", laminar.model$param.names)]
        rorb_markup_tp <- data.frame(
          tp_ref = tpoints_ref,
          tp_18 = tpoints_ref + tpoints_18,
          tpy = rep(-rise, length(tpoints_ref)),
          tpyend = rep(0, length(tpoints_ref))
        )
        # Find rates
        Rates_block3 <- c(
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"], # ref level 
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"] + 
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk3", laminar.model$param.names)], # affect of age
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"] + 
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_right_X_Tns/Blk3", laminar.model$param.names)], + # affect of hemisphere
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"] + 
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk3", laminar.model$param.names)] +    # affect of age
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_right_X_Tns/Blk3", laminar.model$param.names)] + # affect of hemisphere
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_right18_X_Tns/Blk3", laminar.model$param.names)] # affect of hemisphere and age
        )
        Rates_block3 <- exp(Rates_block3) - 1
        rorb_markup_Rt <- data.frame(
          Rts_ref = Rates_block3,
          Rtsx = rep(57, length(Rates_block3)),
          Rtsxend = rep(57 + 4, length(Rates_block3))
        )
        # Find slope
        #  ... the slope M outside of log space equals the slope m inside log space times  
        #       the rate R outside log space plus 1, i.e., M = m * (R + 1)
        #       ... Why? m = dr/dx = d(log(R+1)/dx = dR/dx * 1/(R+1)
        P12_rise <- c(
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk2"] - laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk1"],
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"] - laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk2"],
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk4"] - laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"]
        )
        P18_rise <- c(
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk2"] + 
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk2", laminar.model$param.names)] - 
            (laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk1"] + 
               laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk1", laminar.model$param.names)]),
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"] + 
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk3", laminar.model$param.names)] - 
            (laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk2"] + 
               laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk2", laminar.model$param.names)]),
          laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk4"] + 
            laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk4", laminar.model$param.names)] - 
            (laminar.model$fitted.parameters["baseline_cortex_Rt_Rorb_Tns/Blk3"] + 
               laminar.model$fitted.parameters[grepl("beta_Rt_cortex_Rorb_18_X_Tns/Blk3", laminar.model$param.names)])
        )
        P12_slope_scalar <- c(
          laminar.model$fitted.parameters["baseline_cortex_tslope_Rorb_Tns/Blk1"],
          laminar.model$fitted.parameters["baseline_cortex_tslope_Rorb_Tns/Blk2"],
          laminar.model$fitted.parameters["baseline_cortex_tslope_Rorb_Tns/Blk3"]
        )
        P18_slope_scalar <- c(
          laminar.model$fitted.parameters["baseline_cortex_tslope_Rorb_Tns/Blk1"] + 
            laminar.model$fitted.parameters[grepl("beta_tslope_cortex_Rorb_18_X_Tns/Blk1", laminar.model$param.names)],
          laminar.model$fitted.parameters["baseline_cortex_tslope_Rorb_Tns/Blk2"] + 
            laminar.model$fitted.parameters[grepl("beta_tslope_cortex_Rorb_18_X_Tns/Blk2", laminar.model$param.names)],
          laminar.model$fitted.parameters["baseline_cortex_tslope_Rorb_Tns/Blk3"] + 
            laminar.model$fitted.parameters[grepl("beta_tslope_cortex_Rorb_18_X_Tns/Blk3", laminar.model$param.names)]
        )
        P12_log_slope <- P12_rise*P12_slope_scalar/4
        P18_log_slope <- P18_rise*P18_slope_scalar/4
        P12_slope <- P12_log_slope * (df12$pred[round(tpoints_ref, 0)] + 1)
        P18_slope <- P18_log_slope * (df18$pred[round(tpoints_ref + tpoints_18,0)] + 1)
        P12_run <- rise/P12_slope
        P18_run <- rise/P18_slope
        ratecount_plots[[p_]] <- ratecount_plots[[p_]] +
          geom_segment(
            data = rorb_markup_tp,
            aes(x = tp_ref - P12_run, xend = tp_ref, y = tpy, yend = tpyend),
            color = colors4[1], linetype = "solid", linewidth = 1.5,
            arrow = arrow(length = unit(0.15, "inches"), type = "closed")
          ) +
          geom_segment(
            data = rorb_markup_tp,
            aes(x = tp_18 - P18_run, xend = tp_18, y = tpy, yend = tpyend),
            color = colors4[3], linetype = "solid", linewidth = 1.5,
            arrow = arrow(length = unit(0.15, "inches"), type = "closed")
          ) +
          geom_segment(
            data = rorb_markup_Rt,
            aes(x = Rtsx, xend = Rtsxend, y = Rts_ref, yend = Rts_ref),
            color = c(colors4[1], colors4[3], colors4[2], colors4[4]), linetype = "solid", linewidth = 1.5,
            arrow = arrow(length = unit(0.15, "inches"), type = "closed")
          )
      }
    }
    other_gene_col <- arrangeGrob(
      ratecount_plots[[2]], # Cux2
      ratecount_plots[[6]], # Satb2
      ratecount_plots[[3]], # Fezf2
      ratecount_plots[[4]], # Nxph3
      ratecount_plots[[1]], # Bcl11b
      ncol = 1
    )
    Rorb_ran_block <- arrangeGrob(
      ratecount_plots_Rorb_ran[[1]], 
      ratecount_plots_Rorb_ran[[2]],
      ratecount_plots_Rorb_ran[[3]],
      ratecount_plots_Rorb_ran[[4]], 
      ncol = 2
    )
    Rorb_col <- arrangeGrob(Rorb_ran_block, ratecount_plots[[5]], ncol = 1)
    grid.arrange(Rorb_col, other_gene_col, widths = c(1, 0.5), ncol = 2)
    
    # export
    dev.copy(png, filename = "fig_results_ratecount.png", width = 1500, height = 1215)
    dev.off()
    
  }
make_fig_results_ratecount()

make_fig_residuals <- function() {
    
    # Grab plots 
    hist_plot <- laminar.model[["plots"]][["residuals"]][["all_hist"]]
    qq_plot <- laminar.model[["plots"]][["residuals"]][["all_qq"]]
    
    # Change titles 
    hist_plot <- hist_plot + 
      labs(title = "") 
    qq_plot <- qq_plot + 
      labs(title = "")
    
    grid.arrange(hist_plot, qq_plot, ncol = 2, top = textGrob(
      "Log-Linked Residuals", 
      gp = gpar(fontsize = 30) 
    ))
    
    # export
    dev.copy(png, filename = "fig_residuals.png", width = 1365, height = 605)
    dev.off()
    
  }
make_fig_residuals()

make_MCMCbs_comparison <- function() {
    
    # Grab plots 
    walks_low <- laminar.model[["plots"]][["MCMC"]][["plot.walks.parameters_low"]]
    walks_high <- laminar.model[["plots"]][["MCMC"]][["plot.walks.parameters_high"]]
    walks_nll <- laminar.model[["plots"]][["MCMC"]][["plot.walks.nll"]]
    
    plot_autocor <- laminar.model[["plots"]][["parameter.normality"]][["plot_sample_correlations"]]
    plot_shaprio <- laminar.model[["plots"]][["parameter.normality"]][["plot_comparison_Shaprio"]]
    plot_density <- laminar.model[["plots"]][["parameter.normality"]][["plot_comparison_density"]]
    
    # Resize titles 
    title_size <- 20 
    walks_low <- walks_low + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    walks_high <- walks_high + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    walks_nll <- walks_nll + theme(plot.title = element_text(hjust = 0.5, size = title_size), legend.position = "none")
    plot_autocor <- plot_autocor + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    plot_shaprio <- plot_shaprio + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    plot_density <- plot_density + theme(plot.title = element_text(hjust = 0.5, size = title_size))
    
    # Arrange plots
    grid.arrange(
      walks_low, plot_density,
      walks_high,  plot_shaprio,
      walks_nll, plot_autocor,
      ncol = 2,
      top = textGrob(
        "MCMC and Boostrap Comparison", 
        gp = gpar(fontsize = 24) 
      )
    )
    
    # export
    dev.copy(png, filename = "fig_MCMCbs_comparison.png", width = 1500, height = 950)
    dev.off()
    
  }
make_MCMCbs_comparison()

# Make results table for stats #####
# ... This function is a mess, but it's a one-off for this data set and this paper

library(dplyr)
library(knitr)
library(kableExtra)

make_stat_table <- function() {
    
    param_stats <- laminar.model[["stats"]][["parameters"]][,-c(5,7)]
    param_stats[,2:5] <- round(param_stats[,2:5], 4)
    param_stats <- param_stats[!grepl("wfactor", param_stats$parameter) & !grepl("baseline", param_stats$parameter),]
    param_stats_list <- list() 
    for (g in gene.list) {
      # Grab gene mask
      gene_mask <- grepl(g, param_stats$parameter)
      # Split into lists by gene 
      param_stats_g <- param_stats[gene_mask,]
      # Split parameter names into columns to reorganize data
      split_cols_fix <- do.call(rbind, strsplit(param_stats_g$parameter, "_")) # want cols 2,5,7 (spatial param, treatment, block)
      # Reorganize data
      results_cols <- c(2:6)
      col_names <- c("est", "CI.low", "CI.high", "p.adj", "sig")
      param_stats_list_g <- list()
      for (trt in unique(split_cols_fix[,5])) {
        treatment_mask <- split_cols_fix[,5] == trt
        results <- param_stats_g[treatment_mask, results_cols]
        results[,1:4] <- round(results[,1:4], 3)
        param_type <- split_cols_fix[treatment_mask,2]
        param_type[param_type == "Rt"] <- "r"
        param_type[param_type == "tpoint"] <- "p"
        param_type[param_type == "tslope"] <- "s"
        block <- split_cols_fix[treatment_mask,7]
        block <- gsub("Tns/Blk", "", block)
        if (g == "Rorb") {
          # extend param_type and block by 1
          old_mask <- c()
          param_type_new <- c() 
          block_new <- c()
          last_param_type <- "z"
          last_block <- 0
          for (ti in seq_along(param_type)) {
            t <- param_type[ti]
            b <- as.integer(block[ti])
            if (t == last_param_type || length(block_new) == 0) {
              param_type_new <- c(param_type_new, t)
              block_new <- c(block_new, as.character(b))
              old_mask <- c(old_mask, TRUE)
            } else {
              param_type_new <- c(param_type_new, last_param_type, t)
              block_new <- c(block_new, as.character(last_block+1), as.character(b))
              old_mask <- c(old_mask, FALSE, TRUE)
            }
            if (ti == length(param_type)) {
              param_type_new <- c(param_type_new, t)
              block_new <- c(block_new, as.character(b+1))
              old_mask <- c(old_mask, FALSE)
            }
            last_param_type <- t
            last_block <- b
          }
          param_type <- param_type_new
          block <- block_new
          new_results <- as.data.frame(array(NA, dim = c(length(block), ncol(results))))
          new_results[old_mask,] <- as.data.frame(results[,])
          results <- new_results
        }
        block <- paste0("$", param_type, "_", block, "$")
        treatment <- rep(trt, length(results[,1]))
        param_stats_list_g[[trt]] <- as.data.frame(cbind(treatment, block, results))
        colnames(param_stats_list_g[[trt]]) <- c("effect", "$z$", col_names)
      }
      param_stats_list[[g]] <- as.data.frame(do.call(rbind, param_stats_list_g))
      rownames(param_stats_list[[g]]) <- NULL
      # For each gene and parameter type (ran, fix, baseline), will have five numeric value columns: est, low, high, p-adj, significance
    }
    
    # Sample data
    df <- param_stats_list[["Rorb"]]
    
    effect_lengths <- c()
    for (e in unique(df$effect)) {
      effect_lengths <- c(effect_lengths, sum(df$effect == e))
    }
    names(effect_lengths) <- c(
      "Fixed Effects: Hemisphere (right)", 
      "Fixed Effect: Age (P18)", 
      "Fixed Effect: Hemisphere-Age Interaction"
    )
    df <- df[,-which(colnames(df) == "effect")]
    df$p.adj[df$p.adj > 1] <- "1.000"
    df$p.adj[df$p.adj <= 0] <- "$<0.001$"
    
    for (g in gene.list) {
      if (g == "Rorb") next
      df2 <- param_stats_list[[g]]
      df2 <- df2[,-which(colnames(df2) == "effect")]
      df2$p.adj[df2$p.adj > 1] <- "1.000"
      df2$p.adj[df2$p.adj <= 0] <- "$<0.001$"
      df <- cbind(
        df,
        rep(" ", nrow(df)),
        rep(" ", nrow(df))
      )
      colnames(df) <- c(colnames(df)[1:(ncol(df)-2)], paste0("est.",g), paste("p.adj.",g))
      for (p in 1:effect_lengths[1]) {
        mask <- df[p,1] == df2[,1]
        if (any(mask)) {
          df[df[p,1] == df[,1],c(ncol(df)-1, ncol(df))] <- df2[mask,c("est", "p.adj")]
        }
      }
    }
    
    # Create the table with row grouping
    kbl(df, format = "latex", 
        booktabs = TRUE, escape = FALSE, 
        caption = "Fixed effect estimates.\\label{table:FEestimates}", 
        linesep = "") %>%
      add_header_above(
        c(
          " ", 
          "RORB" = 5,
          "Bcl11b" = 2, "Fezf2" = 2,  "Satb2" = 2,  "Nxph3" = 2,  "Cux2" = 2
        )) %>%
      group_rows(index = effect_lengths) %>%
      kable_styling(latex_options = c("scale_down"), font_size = 8)
    
  }
make_stat_table()

# end sink
sink(file = NULL)

