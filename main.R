
# Setup ################################################################################################################
# Analysis of MERFISH data with Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)

# Clear global environment
rm(list = ls())
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

# For for WSPmm and snk ("sink") printing
# ... version 1.0 used for the preprint (draft 1)
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

# Set file path, and bootstrap chunk size
data_path <- paste0(projects_folder, "_molecular_mechanisms_of_ACx_lateralization/data_SSp/")
bs_chunksize <- 100

# Preprocessing MERFISH data ###########################################################################################

# Define list of genes to analyze
laminar.gene.list <- c("Bcl11b", "Fezf2", "Satb2", "Nxph3", "Cux2", "Rorb")  

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
extract_transcript_counts <- function(gene.list) {
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
extract_transcript_counts(laminar.gene.list)

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
    gene.list = laminar.gene.list,
    fixed.effect.names = fixed.effect.names, 
    context = "cortex"
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
    context = "cortex", 
    species = "gene",
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
    MCMC.steps = 1e3,
    MCMC.step.size = 0.1,
    MCMC.prior = 1.0, 
    MCMC.neighbor.filter = 2
  )

# Fit model
laminar.model <- wisp(
    count.data = count.data.WSPmm,
    variables = data.variables,
    use.median = FALSE,
    bootstraps.num = 1e3,
    converged.resamples.only = TRUE,
    max.fork = bs_chunksize,
    verbose = TRUE,
    plot.settings = list(
      print.plots = TRUE, 
      dim.bounds = colMeans(layer.boundary.bins)
    ),
    MCMC.settings = MCMC.settings,
    model.settings = model.settings
  )

this_gene <- "Rorb"
gene_mask <- grepl(this_gene, names(laminar.model[["plots"]][["ratecount"]]))
n_plots <- length(laminar.model[["plots"]][["ratecount"]])
for (p in 1:n_plots) {
  if (gene_mask[p]) print(laminar.model[["plots"]][["ratecount"]][[p]])
}

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
    title_size <- 30 
    axis_size <- 24 
    legend_size <- 20
    ratecount_plots_Rorb_ran <- list()
    mice <- c("mouse 1", "mouse 2", "mouse 3", "mouse 4")
    length(ratecount_plots_Rorb_ran) <- 4
    names(ratecount_plots_Rorb_ran) <- mice
    for (p in 1:4) {
      
      ratecount_plots_Rorb_ran[[mice[p]]] <- plot.ratecount(
        wisp.results = laminar.model,
        pred.type = "pred",
        count.type = "count",
        dim.boundaries = unlist(laminar.model[["plots"]][["ratecount"]][["plot_pred_context_cortex_fixEff_Bcl11b"]][["layers"]][[49]][["data"]]),
        #y.lim = c(0, 215),
        count.alpha.none = 0,
        count.alpha.ran = 0.5,
        pred.alpha.none = 0,
        pred.alpha.ran = 1,
        rans.to.print = as.character(p),
        speciess.to.print = c("Rorb")
      )[[1]] +
        labs(title = mice[p], x = NULL, y = NULL) + 
        theme(
          plot.title = element_text(hjust = 0.5, size = title_size),
          axis.title = element_text(size = axis_size),
          axis.text = element_text(size = axis_size),
          legend.title = element_text(size = legend_size),
          legend.text = element_text(size = legend_size),
          legend.position = "none"
        ) + 
        scale_color_manual(
          labels = c("left, P12", "right, P12", "left, P18", "right, P18"),
          values = colors4
        )
      
    }
    for (p in 2:n_plots) {
      
      p_ <- names(laminar.model[["plots"]][["ratecount"]])[p]
      gene_name <- gsub("plot_pred_context_cortex_fixEff_", "", p_)
      
      # Check if this is Rorb
      this_Rorb <- FALSE
      if (gene_name == "Rorb") this_Rorb <- TRUE
      
      # Set legend position 
      leg_pos <- "none"
      
      # Reformat gene names
      if (this_Rorb) {
        gene_name <- expression("ROR" * beta)
      } else {
        gene_name <- toupper(gene_name)
      }
      
      # Recolor plot
      ratecount_plots[[p_]] <- laminar.model[["plots"]][["ratecount"]][[p_]] 
      
      # Remake Rorb 
      if (this_Rorb) {
        
        rorb_decomp<- plot.ratecount(
          wisp.results = laminar.model,
          pred.type = "pred",
          count.type = "count",
          dim.boundaries = unlist(laminar.model[["plots"]][["ratecount"]][["plot_pred_context_cortex_fixEff_Bcl11b"]][["layers"]][[49]][["data"]]),
          count.alpha.none = 0.5,
          count.alpha.ran = 0,
          pred.alpha.none = 1,
          pred.alpha.ran = 0,
          rans.to.print = "none",
          speciess.to.print = c("Rorb")
        )
        Rorb_none <- rorb_decomp[[1]]  +
          labs(title = gene_name) + 
          theme(
            plot.title = element_text(hjust = 0.5, size = title_size),
            axis.title = element_text(size = axis_size),
            axis.text = element_text(size = axis_size),
            legend.title = element_text(size = legend_size),
            legend.text = element_text(size = legend_size),
            legend.position = "top"
          ) + 
          scale_color_manual(
            name = "",
            labels = c("left, P12", "right, P12", "left, P18", "right, P18"),
            values = colors4
          )
        
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
          labels = c("left, P12", "right, P12", "left, P18", "right, P18"),
          values = colors4
        )
      
      if (this_Rorb) {
        # Extract data frame from plot 
        found_P12 <- FALSE 
        this_layer <- 1
        while(!found_P12) {
          df12 <- Rorb_none[["layers"]][[this_layer]][["data"]]
          if(all(df12$ran == "none" & df12$treatment == "ref")) found_P12 <- TRUE
          else this_layer <- this_layer + 1
        }
        found_P18 <- FALSE 
        this_layer <- 1 
        while(!found_P18) {
          df18 <- Rorb_none[["layers"]][[this_layer]][["data"]]
          if(all(df18$ran == "none" & df18$treatment == "ref")) found_P18 <- TRUE
          else this_layer <- this_layer + 1
        }
        # Find t-points
        rise <- 10
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
        Rorb_none <- Rorb_none +
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
    
    # First figure
    other_gene_col <- arrangeGrob(
      ratecount_plots[[2]], # Cux2
      ratecount_plots[[6]], # Satb2
      ratecount_plots[[3]], # Fezf2
      ratecount_plots[[5]], # Rorb
      ratecount_plots[[4]], # Nxph3
      ratecount_plots[[1]], # Bcl11b
      ncol = 3
    )
    
    png("fig_results_ratecount.png", width = 1800, height = 1100)
    grid.arrange(other_gene_col, ncol = 1)
    dev.off()
    
    # Second figure
    Rorb_ran_block <- arrangeGrob(
      ratecount_plots_Rorb_ran[[1]], 
      ratecount_plots_Rorb_ran[[2]],
      ratecount_plots_Rorb_ran[[3]],
      ratecount_plots_Rorb_ran[[4]], 
      ncol = 2
    )
    
    png("fig_results_ratecount_Rorb.png", width = 1400, height = 1400)
    grid.arrange(Rorb_none, Rorb_ran_block, heights = c(1, 0.5), ncol = 1)
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
    for (g in laminar.gene.list) {
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
    
    for (g in laminar.gene.list) {
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

# Radial data ##########################################################################################################

preprocess_Droin <- FALSE
if (preprocess_Droin) {
  # Load Droin et al. data
  library(R.matlab)
  files <- list.files(path = "Droin_data", pattern = "\\.mat$", full.names = TRUE)
  Droin_data <- lapply(files, readMat)
  names(Droin_data) <- basename(files)
  
  # Subset to genes of interest
  genes_to_keep <- c("glul", "ass1", "arntl", "dbp", "elovl3", "pck1")
  for (fn in basename(files)) {
    gene_idx <- which(unlist(Droin_data[[fn]][["all.genes"]]) %in% genes_to_keep)
    genes_to_keep <- unlist(Droin_data[[fn]][["all.genes"]])[gene_idx]
    Droin_data[[fn]][["seq.data"]] <- Droin_data[[fn]][["seq.data"]][gene_idx,]
    Droin_data[[fn]][["mat.norm"]] <- Droin_data[[fn]][["mat.norm"]][gene_idx,]
    Droin_data[[fn]][["MeanGeneExp"]] <- Droin_data[[fn]][["MeanGeneExp"]][gene_idx,]
    Droin_data[[fn]][["SE"]] <- Droin_data[[fn]][["SE"]][gene_idx,]
    Droin_data[[fn]][["q.vals"]] <- Droin_data[[fn]][["q.vals"]][gene_idx,]
    rownames(Droin_data[[fn]][["seq.data"]]) <- genes_to_keep
    rownames(Droin_data[[fn]][["mat.norm"]]) <- genes_to_keep
    rownames(Droin_data[[fn]][["MeanGeneExp"]]) <- genes_to_keep
    rownames(Droin_data[[fn]][["SE"]]) <- genes_to_keep
    names(Droin_data[[fn]][["q.vals"]]) <- genes_to_keep
    Droin_data[[fn]][["all.genes"]] <- genes_to_keep
  }
  
  # Create count_data frame for use with WSPmm
  count_data <- data.frame()
  for (i in seq_along(Droin_data)) {
    gene_list <- rownames(Droin_data[[i]][["mat.norm"]])
    n_genes <- length(gene_list)
    n_cells <- ncol(Droin_data[[i]][["mat.norm"]])
    n_bins <- ncol(Droin_data[[i]][["Pmat"]])
    
    Pmat_row_sums <- rowSums(Droin_data[[i]][["Pmat"]])
    denormalized_rates <- Droin_data[[i]][["mat.norm"]] * 1e3
    
    # ... for each gene
    for (g in c(1:n_genes)) {
      
      # ... for each zone (i.e., spatial coordinate, or bin)
      for (b in c(1:n_bins)) {
        
        # Normalized count of gene g in each cell, for mouse i
        sim_gene_counts_by_cell <- rpois(n_cells, denormalized_rates[g,])
        
        # Get probability of bin membership for each cell, for mouse i 
        cell_weight <- rbinom(n = n_cells, size = 1, prob = Droin_data[[i]][["Pmat"]][,b] / Pmat_row_sums)
        
        count_data <- rbind(
          count_data,
          data.frame(
            # add mouse number
            mouse = rep(i, n_cells),
            # extract ZT from file name and add
            ZT = rep(as.numeric(gsub("\\D", "", names(Droin_data)[i])), n_cells),
            # add bin number
            bin = rep(b, n_cells),
            # add gene name
            gene = rep(gene_list[g], n_cells),
            # compute and log-transform the zone-weighted normalized gene counts
            count = sim_gene_counts_by_cell * cell_weight
          )
        )
        
      }
    }
  }
  #if (min(count_data$count) < 0) count_data$count <- count_data$count - min(count_data$count)
  
  # Export
  write.csv(count_data, file = "Droin_radial_count_data_sim.csv", row.names = FALSE)
} else {
  # Import
  count_data <- read.csv("Droin_radial_count_data.csv")
}

# Define fixed effects to test
fixed.effect.names <- c("ZT")

# Define variables in the dataframe for the model
data.variables <- list(
  count = "count",
  bin = "bin", 
  context = "liver", 
  species = "gene",
  ran = "mouse",
  timeseries = "ZT",
  fixedeffects = fixed.effect.names
)

# Model settings 
# ... all settings shown here are defaults
model.settings <- list(
  # ... these are global options needed to set up model
  buffer_factor = 0.05,                                 # buffer factor for penalizing distance from structural parameter values
  ctol = 1e-6,                                          # convergence tolerance
  max_penalty_at_distance_factor = 0.01,                # maximum penalty at distance from structural parameter values
  LROcutoff = 1.5,                                      # cutoff for LROcp, a multiple of standard deviation
  LROwindow_factor = 1.25,                              # window factor for LROcp, larger means larger rolling window
  rise_threshold_factor = 0.8,                          # amount of detected rise as fraction of total required to end run in initial slope estimation
  max_evals = 1000,                                     # maximum number of evaluations for optimization
  rng_seed = 42,                                        # random seed for optimization (controls bootstrap resamples only)
  warp_precision = 1e-7,                                # decimal precision to retain when selecting really big number as pseudo infinity for unbound warping
  round_none = FALSE
)

# Settings for MCMC walk
MCMC.settings <- list(
  MCMC.burnin = 1e2,
  MCMC.steps = 1e2,
  MCMC.step.size = 0.5,
  MCMC.prior = 1.0, 
  MCMC.neighbor.filter = 2
)

# Settings for plotting
plot.settings <- list(
  print.plots = FALSE,
  title_size = 14,
  count_size = 2.5,
  count_jitter = 0.0,
  count.alpha.none = 0.8,
  count.alpha.ran = 0.5,
  pred.alpha.ran = 0.0
)

# Fit model
# ... all settings shown here are defaults
radial.model <- wisp(
  count.data = count_data,
  variables = data.variables,
  bootstraps.num = 1e3,
  #fit_only = TRUE,
  max.fork = bs_chunksize,
  model.settings = model.settings,
  MCMC.settings = list(MCMC.steps = 0),
  plot.settings = plot.settings
)

v1 <- c("#ADD8E6", "#FF7F50", "#9ACD31", "#FFC1CB")
v2 <- c("#ED6986", "#9CA620", "#2AA9A3", "#A28EC2")
new_ratecount <- list()
for (p in seq_along(radial.model[["plots"]][["ratecount"]])) {
  if (p == 1) next
  plt <- radial.model[["plots"]][["ratecount"]][[p]]
  plt <- plt + scale_color_manual(
    values = v2
  ) 
  print(plt)
  new_ratecount[[p - 1]] <- plt
  print(p)
}

all_plots <- c(new_ratecount, radial.model[["plots"]][["timeseries"]])
all_plots <- as.vector(rbind(new_ratecount, radial.model[["plots"]][["timeseries"]]))
do.call(grid.arrange, c(all_plots, ncol = 2))

# Simulated data #######################################################################################################

# For the 250-sim benchmarking run, just run setup, come down to here, and run all lines to the bottom of the 250-sim loop.
set.seed(9999)

preprocess_Allen <- FALSE 
if (preprocess_Allen) {
  library(data.table)
  # Get complete Allen mouse brain data
  data_path <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/_molecular_mechanisms_of_ACx_lateralization/development_work/data/Allen_data/Allen_data.csv"
  start_time <- Sys.time()
  big_data_frame <- fread(data_path,
                          nrows = 1e6 )
  duration <- Sys.time() - start_time
  cat("Time to read in first million rows of data: ", duration, units(duration), "\n")
  start_time <- Sys.time()
  big_data_frame <- fread(data_path)
  duration <- Sys.time() - start_time
  cat("Time to read in data: ", duration, units(duration), "\n")
  View(big_data_frame[c(1:1000),])
  for (s in unique(big_data_frame$slice_num)) {
    mask <- big_data_frame$slice_num == s
    print(sum(mask)/1e6)
  }
  
  # Prune down to just a few interesting neuron markers
  neuron.gene.list <- c("Slc17a6", "Slc17a7", "Slc17a8", "Pvalb", "Tac2", "Vip")
  neuron.mask <- FALSE 
  for (g in neuron.gene.list) {
    neuron.mask <- neuron.mask | big_data_frame$trscrpt_gene_symb == g
  }
  print(sum(neuron.mask))
  count_data_neurons <- big_data_frame[neuron.mask,]
  
  # Collect more genes for enriched simulated reference atlas 
  set.seed(999) # ... all 20 remain with this seed
  enriched.gene.list <- sample(unique(big_data_frame$trscrpt_gene_symb), 20)
  enriched.gene.list <- enriched.gene.list[!enriched.gene.list %in% neuron.gene.list]
  
  # Make enriched data frame
  for (g in enriched.gene.list) {
    neuron.mask <- neuron.mask | big_data_frame$trscrpt_gene_symb == g
  }
  print(sum(neuron.mask))
  count_data_neurons_enriched <- big_data_frame[neuron.mask,]
  
  # Drop unneeded columns 
  count_data_neurons <- count_data_neurons[,c(1,2,3,6,8,9)]
  count_data_neurons_enriched <- count_data_neurons_enriched[,c(1,2,3,6,8,9)]
  
  # Rename columns 
  colnames(count_data_neurons) <- c("count", "cell_id", "slice_num", "gene", "coord_x", "coord_y")
  colnames(count_data_neurons_enriched) <- colnames(count_data_neurons)
  
  # Identify which vglut marker to keep
  for (g in neuron.gene.list) {
    mask <- count_data_neurons$gene == g
    cat("\ngene:", g, "count:", sum(count_data_neurons$count[mask]))
  }
  
  # Keep only Slc17a7
  count_data_neurons <- count_data_neurons[!c(count_data_neurons$gene == "Slc17a6" | count_data_neurons$gene == "Slc17a8"),]
  count_data_neurons_enriched <- count_data_neurons_enriched[!c(count_data_neurons_enriched$gene == "Slc17a6" | count_data_neurons_enriched$gene == "Slc17a8"),]
  
  # Identify slices to form the basis of the simulation
  start <- 40
  stop <- 50
  for (s in c(start:stop)) {
    df <- count_data_neurons[count_data_neurons$slice_num == s, c("gene", "coord_x", "coord_y")]
    plt <- ggplot(df, aes(x = coord_x, y = coord_y, color = gene)) +
      geom_point() + 
      ggtitle(as.character(s)) +
      theme_minimal()
    print(plt)
  }
  
  # Keep only these slices
  mask <- count_data_neurons$slice_num == 29 | count_data_neurons$slice_num == 33 | count_data_neurons$slice_num == 42
  count_data_neurons <- count_data_neurons[mask,]
  mask <- count_data_neurons_enriched$slice_num == 29 | count_data_neurons_enriched$slice_num == 33 | count_data_neurons_enriched$slice_num == 42
  count_data_neurons_enriched <- count_data_neurons_enriched[mask,]
  
  # Flip y-axis
  y_med <- median(count_data_neurons$coord_y)
  count_data_neurons$coord_y <- (count_data_neurons$coord_y - y_med) * -1 + y_med
  count_data_neurons_enriched$coord_y <- (count_data_neurons_enriched$coord_y - y_med) * -1 + y_med
  
  # Inspect data
  plt <- ggplot(count_data_neurons, aes(x = coord_x, y = coord_y, color = gene)) +
    geom_point(size = 0.5) + 
    facet_wrap(~ slice_num) +
    theme_minimal()
  print(plt)
  
  # Save count data
  write.csv(count_data_neurons, "Allen_data.csv", row.names = FALSE)
  
  # Create the reference count matrix 
  dt <- as.data.table(count_data_neurons_enriched)
  # ... sum counts per gene × cell_id
  dt <- dt[, .(count = sum(count)), by = .(gene, cell_id)]
  # pivot to wide format
  ref_counts <- dcast(
    dt,
    gene ~ cell_id,
    value.var = "count",
    fill = 0
  )
  # Subset to 10k random columns (cells) 
  count_idx <- which(colSums(ref_counts[,-1, with = FALSE]) >= 100) + 1
  set.seed(42123)
  ref_counts <- ref_counts[, c(1, sample(count_idx, 1e4, replace = FALSE)), with = FALSE]
  # convert to matrix with rownames as genes
  ref_counts <- as.matrix(ref_counts[, -1, with = FALSE])
  rownames(ref_counts) <- dt[, sort(unique(gene))]
  
  # Save 
  write.csv(ref_counts, "Allen_reference_counts_enriched_10kcells.csv", row.names = TRUE)
} else {
  # Import
  count_data_neurons <- read.csv("Allen_data.csv")
}

# Remove cell_id column from count_data_neurons 
count_data_neurons <- count_data_neurons[,c("count", "slice_num", "gene", "coord_x", "coord_y")]

# Cut down to just a single 2x2 patch from one slice
count_data_neurons_patch <- count_data_neurons[
  count_data_neurons$slice_num == 33 &
    count_data_neurons$coord_y >= 2 & count_data_neurons$coord_y <= 4 &
    count_data_neurons$coord_x >= 1 & count_data_neurons$coord_x <= 3,]

# DESeq2 benchmarking functions ########################################################################################

# Load library 
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")}
if (!require("apeglm", quietly = TRUE)) {BiocManager::install("apeglm")}
library(DESeq2) 

# Helper function to convert sim data into format expected by DESeq2
make_DESeq2_data <- function(
    sim_data
  ) {
    # Set reference level
    sim_data$fixedeffect <- relevel(as.factor(sim_data$fixedeffect), ref = "ref")
    # Get genes 
    genes <- unique(sim_data$gene)
    n_genes <- length(genes)
    # Grab variable levels
    fixedeffects <- unique(sim_data$fixedeffect)
    replicates <- unique(sim_data$replicate)
    sample_names <- c(
      paste0(as.character(replicates), "_", as.character(fixedeffects)[1]),
      paste0(as.character(replicates), "_", as.character(fixedeffects)[2])
    )
    # Make coldata
    coldata <- data.frame(
      fixedeffect = rep(fixedeffects, each = length(replicates)), 
      replicate = rep(replicates, length(fixedeffects))
    )
    rownames(coldata) <- sample_names
    # Make count matrix 
    cts <- array(NA, dim = c(n_genes, nrow(coldata)))
    rownames(cts) <- genes
    colnames(cts) <- rownames(coldata)
    for (i in genes) {
      for (j in rownames(coldata)) {
        ct_mask <- sim_data$replicate == coldata[j,"replicate"] & 
          sim_data$fixedeffect == coldata[j,"fixedeffect"] & 
          sim_data$gene == i
        cts[i,j] <- sum(sim_data[ct_mask, "count"])
      }
    }
    return(list(cts = cts, coldata = coldata))
  }

# Function to model simulation with DESeq2
model_attractor_simulation_DESeq2 <- function(
    sim,
    sim_num
  ) {
    
    # Make data for DESeq2
    ddata <- make_DESeq2_data(sim$data)
    dds <- DESeqDataSetFromMatrix(countData = ddata$cts, colData = ddata$coldata, design = ~ fixedeffect)
    
    # Fit model
    # ... run differential expression analysis on fixedeffect
    dds <- DESeq(estimateSizeFactors(dds), fitType = "mean", minReplicatesForReplace = 3)
    # ... apply shrinkage (ridge regression)
    resLFC_fixedeffect <- lfcShrink(dds, coef = "fixedeffect_trt_vs_ref", type = "apeglm")
    
    # Extract model results for comparing to ground truth
    MR <- as.data.frame(resLFC_fixedeffect@listData)[,c("log2FoldChange", "padj")]
    
    # Extract ground-truth from simulation
    GT <- attractor_simulation_ground_truth(sim)
    
    # Compile results in named vector and return
    return(
      data.frame(
        est = c(MR$log2FoldChange, MR$padj),
        true = c(GT$fse_true, GT$FSEs),
        param = c(
          rep("rate_effect", length(MR$log2FoldChange)), 
          rep("FSE", length(MR$padj))
        ),
        id = c(
          rownames(MR), 
          rownames(MR)
        ),
        method = "DESeq2",
        sim = sim_num
      )
    )
   
  }

# ELLA benchmarking functions ##########################################################################################

# Set up virtual python environment for reticulate: 
#   python3 -m venv ~/.virtualenvs/r-reticulate
#   source ~/.virtualenvs/r-reticulate/bin/activate
#   pip install --upgrade pip

library(reticulate)
if (Sys.info()["sysname"] == "Linux") {
  use_virtualenv("/home/oviedoworkstation/.virtualenvs/r-reticulate", required = TRUE)
} else {
  use_virtualenv("/Users/michaelbarkasi/.virtualenvs/r-reticulate", required = TRUE)
}

py_install(c(
  "numpy",
  "pandas",
  "matplotlib",
  "scipy",
  "scikit-learn",
  "torch",
  "tqdm",
  "ipdb",
  "rpy2"
), pip = TRUE)

# Install ELLA:
# git clone --branch ella1 https://github.com/jadexq/ELLA.git
# cd ELLA
# pip install .

# For debugging
#options(reticulate.python.stdout = TRUE)
#options(reticulate.python.stderr = TRUE)

# import ELLA module
ELLA_mod <- import_from_path("ELLA", path = "./ELLA/ELLA")

# Make function to convert simulation data into ELLA format
convert_sim_data_to_ELLA <- function(
    sim_data,
    radial_dim = "bin_x",
    theta = "bin_y"
  ) {
    
    # Initialize top-level list structure
    ella_data <- list()
    length(ella_data) <- 7
    names(ella_data) <- c(
      "types",        # character vector of cell types
      "cells",        # named list, each element is a character vector of cell IDs for each type (with the type as the list element name)
      "cells_all",    # character vector of all cell IDs
      "genes",        # named list, each element is a character vector of gene names for each type (with the type as the list element name)
      "cell_seg",     # data frame with columns x, y, and cell with points giving the cell boundary segmentation
      "nucleus_seg",  # data frame with columns x, y, and cell with points giving the nucleus boundary segmentation
      "expr"          # data frame with columns "type", "cell", "gene", "x", "y", "umi", "centerX", "centerY", "sc_total"
    )
    
    # Set types
    # ... use "ref" from fixedeffect as the cell type, throw out "trt" (because no ground truth about SVG for trt)
    ella_data$types <- list("ref")
    
    # Throw out trt data
    sim_data <- sim_data[sim_data$fixedeffect == "ref",]
    
    # Set cells and cells_all
    cell_names <- sort(unique(sim_data$replicate))
    ella_data$cells <- list(ref = cell_names)
    ella_data$cells_all <- cell_names
    
    # Set genes
    ella_data$genes <- list(ref = sort(unique(sim_data$gene)))
    
    # Make cell segmentation data 
    x_range <- range(sim_data$bin_x)
    y_range <- range(sim_data$bin_y)
    n_segs <- 2 * (x_range[2] - x_range[1] + 1) + 2 * (y_range[2] - y_range[1] + 1)
    n_cells <- length(cell_names)
    ella_data$cell_seg <- data.frame(
      x = rep(c(
        rep(x_range[1], y_range[2] - y_range[1] + 1),
        rep(x_range[2], y_range[2] - y_range[1] + 1),
        seq(x_range[1], x_range[2]),
        seq(x_range[1], x_range[2])
      ), n_cells),
      y = rep(c(
        seq(y_range[1], y_range[2]),
        seq(y_range[1], y_range[2]),
        rep(y_range[1], x_range[2] - x_range[1] + 1),
        rep(y_range[2], x_range[2] - x_range[1] + 1)
      ), n_cells),
      cell = rep(cell_names, each = n_segs)
    )
    
    # wrap patch around nuclear center 
    r <- sim_data[,radial_dim]
    theta <- sim_data[,theta]
    r <- r - min(r) + 1
    theta <- theta - min(theta) + 1
    X <- r * cos((theta/max(theta)) * 2 * pi)
    Y <- r * sin((theta/max(theta)) * 2 * pi)
    
    # wrap cell segmentation around nuclear center
    r_seg_dim <- substr(radial_dim, nchar(radial_dim), nchar(radial_dim))
    theta_seg_dim <- "y"
    if (r_seg_dim == "y") theta_seg_dim <- "x"
    theta_seg <- ella_data$cell_seg[,theta_seg_dim]
    theta_seg <- theta_seg - min(theta_seg) + 1
    r_seg <- ella_data$cell_seg[,r_seg_dim]
    r_seg <- r_seg - min(r_seg) + 1
    ella_data$cell_seg[,r_seg_dim] <- r_seg * cos((theta_seg/max(theta_seg)) * 2 * pi)
    ella_data$cell_seg[,theta_seg_dim] <- r_seg * sin((theta_seg/max(theta_seg)) * 2 * pi)
    
    # Compensate for radial dilution
    max_count <- max(sim_data$count)
    count <- sim_data$count * r
    count <- (count / max(count)) * max_count
    count <- as.integer(round(count, 0))
    
    # Get library counts 
    sc_total <- rep(0, nrow(sim_data))
    for (c in cell_names) {
      mask <- sim_data$replicate == c
      sc_total[mask] <- sum(sim_data$count[mask])
    }
    
    # Make expr data 
    ella_data$expr <- data.frame(
      type = sim_data$fixedeffect,
      cell = sim_data$replicate,
      gene = sim_data$gene,
      x = X,
      y = Y,
      umi = count,
      centerX = 0,
      centerY = 0,
      sc_total = sc_total
    )
    
    return(ella_data)
    
  }

# Function to run ELLA on simulated data
run_ELLA <- function(
    sim_data
  ) { 
    
    # construct ELLA object
    ELLA_class <- ELLA_mod$ELLA
    ella_sim <- ELLA_class(
      dataset = "sim_benchmark",
      adam_learning_rate_min = 1e-2,
      max_iter = 1000L, 
      L1_lam = 0.2
    )
    
    # convert simulated data to ELLA format
    ella_data <- r_to_py(convert_sim_data_to_ELLA(sim_data))
    # load data into ELLA object
    ella_sim$load_data(data_dict = ella_data)
    # register cells
    ella_sim$register_cells()
    # prepare data
    ella_sim$nhpp_prepare()
    # fit nhpp model
    ella_sim$nhpp_fit()
    # expression intensity estimation
    ella_sim$weighted_density_est()
    # likelihood ratio test
    ella_sim$compute_pv()
    
    return(ella_sim)
    
  }
#ella_sim <- run_ELLA(sim_data$data)

# Extract SVG p-values from ELLA results
extract_svg <- function(ella_sim) {
  pv_svg <- unlist(ella_sim$pv_fdr_tl[["ref"]])
  names(pv_svg) <- ella_sim$gene_list_dict[["ref"]]
  pv_svg
}
#extract_svg(ella_sim)

# Full model pipeline 
model_attractor_simulation_ELLA <- function(
    sim,
    sim_num 
  ) {
    
    # Run ELLA
    model <- run_ELLA(sim$data)
    
    # Extract model results for comparing to ground truth
    svg_pvalues <- unlist(model$pv_cauchy_tl[["ref"]])
    names(svg_pvalues) <- model$gene_list_dict[["ref"]]
    svg_pvalues
    
    # Extract ground-truth
    GT <- attractor_simulation_ground_truth(sim)
    SVGs <- GT$SVGs
    
    # Compile results in named vector 
    results <- data.frame(
      est = c(svg_pvalues),
      true = c(SVGs),
      param = c(
        rep("SVG", length(svg_pvalues))
      ),
      id = c(
        names(svg_pvalues)
      ),
      method = "ELLA",
      sim = sim_num
    )
    
    return(results)
    
  }

# Compare spatial distribution of gene along original x axis and transformed radial axis
plot_count_density <- function(
    df_rad, 
    df_car, 
    nbins = 25, 
    replicate = "replicate1"
  ) {
    
    ## ---------- filter ----------
    df_rad <- df_rad[df_rad$cell == replicate,]
    df_rad$count <- df_rad$umi
    df_car <- df_car[df_car$replicate == replicate,]
    df_car$x <- df_car$bin_x
    df_car$y <- df_car$bin_y
    
    ## ---------- scatter plots ----------
    p1 <- ggplot(df_car, aes(x = bin_x, y = bin_y, color = gene, size = count)) +
      geom_point() +
      ylab("y") +
      xlab("x") +
      theme_minimal()
    p2 <- ggplot(df_rad, aes(x = x, y = y, color = gene, size = count)) +
      geom_point() +
      theme_minimal() 
    
   
    ## ---------- density plots ----------
    eps <- 4e-3 # For matching the log-scale transforms
    
    ## ---------- density vs x ----------
    xbreaks <- seq(min(df_car$x), max(df_car$x), length.out = nbins + 1)
    xmid    <- 0.5 * (xbreaks[-1] + xbreaks[-length(xbreaks)])
    xwidth  <- diff(xbreaks)
    
    df_car$xbin <- cut(df_car$x, breaks = xbreaks, include.lowest = TRUE)
    
    xcount <- tapply(
      df_car$count,
      list(df_car$gene, df_car$xbin),
      sum
    )
    
    df_x <- data.frame(
      gene    = rep(rownames(xcount), times = ncol(xcount)),
      xmid    = rep(xmid, each = nrow(xcount)),
      density = as.vector(xcount) / rep(xwidth, each = nrow(xcount))
    )
    
    p_x <- ggplot(df_x, aes(x = xmid, y = density + eps, color = gene)) +
      geom_line() +
      scale_y_log10() +
      labs(x = "x", y = "count density (log10)") +
      theme_minimal()
    
    ## ---------- density vs radius ----------
    r <- sqrt(df_rad$x^2 + df_rad$y^2)
    
    rbreaks <- seq(min(r), max(r), length.out = nbins + 1)
    rmid    <- 0.5 * (rbreaks[-1] + rbreaks[-length(rbreaks)])
    rwidth  <- diff(rbreaks)
    
    df_rad$rbin <- cut(r, breaks = rbreaks, include.lowest = TRUE)
    
    rcount <- tapply(
      df_rad$count,
      list(df_rad$gene, df_rad$rbin),
      sum
    )
    
    df_r <- data.frame(
      gene = rep(rownames(rcount), times = ncol(rcount)),
      rmid = rep(rmid, each = nrow(rcount)),
      areal_density =
        as.vector(rcount) /
        (2 * pi * rep(rmid, each = nrow(rcount)) * rep(rwidth, each = nrow(rcount)))
    )
    
    p_r <- ggplot(df_r, aes(x = rmid, y = areal_density + eps, color = gene)) +
      geom_line() +
      scale_y_log10() +
      labs(x = "radius", y = "areal density (log10)") +
      theme_minimal()
    
    ## ---------- arrange ----------
    
    title <- textGrob(
      paste0("Density comparison of radial transform, ", replicate),
      gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    title_scatter <- textGrob(
      paste0("Scatter plot comparison of radial transform, ", replicate),
      gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    panels <- arrangeGrob(
      p_x,
      p_r,
      ncol = 2
    )
    
    panels_scatter <- arrangeGrob(
      p1,
      p2,
      ncol = 2
    )
    
    grid.arrange(
      title_scatter,
      panels_scatter,
      title,
      panels,
      ncol = 1,
      heights = c(0.1, 1, 0.1, 1)
    )
    
  }
#plot_count_density(ella_sim$data_df, sim_data$data)

# Show how ELLA fits the simulated data
plot_ella_fit <- function(
    ella_sim,
    sim_data,
    scalar = 10 # ad hoc, is an approximate guess
  ) {
    
    # Convert spatial coordinate to radial distance
    rad <- round(ella_sim$df_registered$d_c_s,2)*100
    rads <- sort(unique(rad))
    
    # Initialize data frame to hold results 
    n_rows <- length(rads) * length(ella_sim$gene_list_dict[["ref"]])
    counts <- data.frame(gene = character(n_rows), rad = numeric(n_rows), count = integer(n_rows), countc = integer(n_rows))
    
    # Grab simulation data with original Cartesian coordinates
    datac <- sim_data[sim_data$fixedeffect == "ref",]
    
    # Compute counts per radial bin for each gene
    i <- 0
    for (r in rads) {
      print(r)
      mask_r <- rad == r
      mask_rc <- datac$bin_x == r
      for (g in ella_sim$gene_list_dict[["ref"]]) {
        mask <- ella_sim$df_registered$gene == g & mask_r
        maskc <- datac$gene == g & mask_rc
        i <- i + 1
        counts$gene[i] <- g
        counts$rad[i] <- r
        n <- length(unique(ella_sim$df_registered$cell[mask]))
        # ... count after the radial transform
        counts$count[i] <- sum(ella_sim$df_registered$umi[mask])/n
        # ... count before the radial transform 
        counts$countc[i] <- sum(datac$count[maskc])/n
      }
    }
    
    # Grab ELLA predictions
    predictions <- ella_sim$weighted_lam_est[["ref"]]
    names(predictions) <- ella_sim$gene_list_dict[["ref"]]
    pred_wide <- as.data.frame(predictions)
    r <- c(1:nrow(pred_wide))
    pred_wide$r <- r
    pred <- data.frame()
    for (g in names(predictions)) {
      df_g <- data.frame(
        gene = g,
        r = r,
        lam = predictions[[g]]
      )
      pred <- rbind(pred, df_g)
    }
    
    # plot 
    # ... As the ELLA fit accounts for radial dilution, we need to use the undiluted counts pre-transform for comparison
    ggplot(pred, aes(x = r, y = log(lam+1), color = gene)) +
      geom_line() +
      labs(
        title = "ELLA fit vs. counts pre-radial transform",
        x = "radius",
        y = "expression level (log scale)"
      ) +
      geom_point(data = counts, aes(x = rad, y = log(countc/scaler+1), color = gene), alpha = 0.5) +
      theme_minimal()
    
  }
#plot_ella_fit(ella_sim, sim_data$data, scalar = 10)

# Run benchmarking #####################################################################################################

results <- run_attractor_sim_benchmarks(
  seed_data = count_data_neurons_patch,
  n_sims = 250,
  modeling_functions = list(
    wisp = model_attractor_simulation_wisp, 
    ELLA = model_attractor_simulation_ELLA,
    DESeq2 = model_attractor_simulation_DESeq2
    ),
  modeling_function_args = list(wisp = list(bs_num = 1e3, max_fork = bs_chunksize))
)

# Save results
write.csv(results, "benchmark_results.csv", row.names = FALSE)

# Load results 
results <- read.csv("benchmark_results.csv")

res <- results$results
meths <- unique(res$method)
param_mask <- res$param == "FSE" | res$param == "SVG"
for (m in meths) {
  mask <- res$method == m & param_mask
  if (m == "wisp") plot(res$est[mask], main = m)
  else plot(log(-log(res$est[mask])+1), main = m)
  cat("\n", m, ":", sum(res$est[mask] == 0)/sum(mask))
}

analyze_results <- function(
    results,
    sig_thresh = list(wisp = 0.05, ELLA = 0.05, DESeq2 = 0.001)
  ) {
    
    methods <- sort(unique(results$method))
    metrics <- c(
      "cor_rate", 
      "cor_ran", 
      "FPR_svg", 
      "FDR_svg", 
      "power_svg", 
      "FPR_fse", 
      "FDR_fse", 
      "power_fse"
    )
    
    summary <- as.data.frame(matrix(NA, nrow = length(metrics), ncol = length(methods)))
    colnames(summary) <- methods
    rownames(summary) <- metrics
   
    for (m in methods) {
      
      sigt <- 0.05
      if (m %in% names(sig_thresh)) sigt <- sig_thresh[[m]]
      
      m_mask <- results$method == m
      
      # Compute correlation for rate effect parameter, true vs est
      if (any("rate_effect" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "rate_effect"
        summary["cor_rate", m] <- cor(results[r_m_mask,"true"], results[r_m_mask,"est"], method = "pearson")
      }
      # Compute correlation for random effect parameter, true vs est
      if (any("random_effect" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "random_effect"
        summary["cor_ran", m] <- cor(results[r_m_mask,"true"], results[r_m_mask,"est"], method = "pearson")
      }
      # Compute FPR, FDR, power for SVG parameter
      if (any("SVG" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "SVG"
        FP <- sum(results[r_m_mask & results$true == FALSE, "est"] < sigt)
        TN <- sum(results[r_m_mask & results$true == FALSE, "est"] >= sigt)
        TP <- sum(results[r_m_mask & results$true == TRUE, "est"] < sigt)
        FN <- sum(results[r_m_mask & results$true == TRUE, "est"] >= sigt)
        summary["FPR_svg", m] <- FP / (FP + TN)
        summary["FDR_svg", m] <- FP / (FP + TP)
        summary["power_svg", m] <- TP / (TP + FN)
      }
      # Compute FPR, FDR, power for FSE parameter
      if (any("FSE" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "FSE"
        FP <- sum(results[r_m_mask & results$true == FALSE, "est"] < sigt)
        TN <- sum(results[r_m_mask & results$true == FALSE, "est"] >= sigt)
        TP <- sum(results[r_m_mask & results$true == TRUE, "est"] < sigt)
        FN <- sum(results[r_m_mask & results$true == TRUE, "est"] >= sigt)
        summary["FPR_fse", m] <- FP / (FP + TN)
        summary["FDR_fse", m] <- FP / (FP + TP)
        summary["power_fse", m] <- TP / (TP + FN)
      }
      
    }
    return(summary)
  }
results_summary <- analyze_results(results$results)
View(results_summary)


# Analyze results
results_rr <- results[grepl("_effect", results$param),]
ggplot(results_rr, aes(x = true, y = est, color = param)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "WSPmm parameter estimation on simulated data",
    x = "Ground truth",
    y = "Estimated"
  ) +
  scale_color_manual(
    values = c("rate_effect" = "blue", "random_effect" = "red"),
    name = "Parameter"
  )



# type I errors: FPR (FP / (FP + TN)) and FDR (FP / (FP + TP)) ... prob of false positive (i.e., of incorrectly rejecting null h, conditional on the null (FPR) or on a positive results (FDR)
# type II errors: power (TP / (TP + FN)) ... prob of detecting thing (i.e., of correctly rejecting null h), conditional on the alternative h
# 
# Metrics: cor_rate, cor_ran, FPR_svg, FDR_svg, power_svg, FPR_fse, FDR_fse, power_fse

# wisp: cor_rate, cor_ran, FPR_fse, FDR_fse, power_fse
# ELLA: FPR_svg, FDR_svg, power_svg
# C-SIDE
# SpaNorm
# Spark
# DESeq2

# end sink
sink(file = NULL)










