
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
bs_chunksize <- 20

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
    # Data to model
    count.data = count.data.WSPmm,
    # Variable labels
    variables = data.variables,
    # Settings used on R side
    use.median = FALSE,
    MCMC.settings = MCMC.settings,
    bootstraps.num = 1e3,
    converged.resamples.only = TRUE,
    max.fork = bs_chunksize,
    dim.bounds = colMeans(layer.boundary.bins),
    verbose = TRUE,
    print.plots = TRUE,
    # Setting to pass to C++ model
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

library(data.table)

preprocess_Allen <- FALSE 
if (preprocess_Allen) {
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
  ref_counts <- read.csv("Allen_reference_counts_enriched_10kcells.csv", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  ref_counts <- data.table::fread("Allen_reference_counts_enriched_10kcells.csv", data.table = FALSE)
  row_names <- ref_counts[,1]
  ref_counts <- as.matrix(ref_counts[,-1])
  rownames(ref_counts) <- row_names
}

# Remove cell_id column from count_data_neurons 
count_data_neurons <- count_data_neurons[,c("count", "slice_num", "gene", "coord_x", "coord_y")]

# Cut down to just a single 2x2 patch from one slice
count_data_neurons_patch <- count_data_neurons[
  count_data_neurons$slice_num == 33 &
    count_data_neurons$coord_y >= 2 & count_data_neurons$coord_y <=4 &
    count_data_neurons$coord_x >= 1 & count_data_neurons$coord_x <=3,]

# Helper function to make replicates
make_replicate <- function(
    data, 
    rate_scalar = 0.5, # number between 0 and 1 ... should this vary by gene??
    affine_transform = NULL,
    spatial_scalar = 0.05
  ) {
    data_shifted <- data
    # Make affine transformation to shift cells around
    if (is.null(affine_transform)) {
      Ax <- matrix(c(1,0,rnorm(1,0,spatial_scalar),1),2,2)  # shear x
      Ay <- matrix(c(1,rnorm(1,0,spatial_scalar),0,1),2,2)  # shear y
      As <- matrix(c(rnorm(1,0,spatial_scalar)+1, 0, 0, rnorm(1,0,spatial_scalar)+1),2,2) # scale
      A <- As %*% Ay %*% Ax
    } else {
      A <- affine_transform
    }
    # Make scales to keep cells roughly in bounds
    x_mid <- diff(range(data$coord_x))/2 + min(data$coord_x)
    distx <- (data$coord_x - x_mid)^2
    distx <- distx/max(distx)
    y_mid <- diff(range(data$coord_y))/2 + min(data$coord_y)
    disty <- (data$coord_y - y_mid)^2
    disty <- disty/max(disty)
    # Applied scaled transformation
    data_shifted$coord_x = distx*data$coord_x + (1-distx)* (A[1,1]*data$coord_x + A[1,2]*data$coord_y)
    data_shifted$coord_y = disty*data$coord_y + (1-disty)* (A[2,1]*data$coord_x + A[2,2]*data$coord_y)
    # Poisson resampling of genes
    data_shifted$count <- rpois(nrow(data), data_shifted$count * (2*rate_scalar))
    return(data_shifted)
  }

# Helper function to simulate non-SVGs
smooth_data <- function(
    data
  ) {
    idx_shuffled <- sample(nrow(data), nrow(data), replace = FALSE)
    data$coord_x[idx_shuffled] <- data$coord_x
    data$coord_y[idx_shuffled] <- data$coord_y
    return(data)
  }

# Helper function to simulate SV from an attractor point and spatial effect
induce_SV <- function(
    data_mixed, 
    gene, 
    attractor,    # number between 0 and 1
    effect = 0.5  # number between 0 and 1, to simulate fixed effects on rate
  ) {
    # Get attractor coordinates
    n_rows <- nrow(data_mixed)
    attractor_idx <- as.integer(n_rows * attractor)
    attractor_idx <- min(max(1, attractor_idx), n_rows)
    coord_attractor <- c(data_mixed$coord_x[attractor_idx], data_mixed$coord_y[attractor_idx])
    # Find distances from cells (i.e., transcripts) to attractor points
    gene_mask <- data_mixed$gene == gene
    x_diff <- data_mixed$coord_x[gene_mask] - coord_attractor[1]
    y_diff <- data_mixed$coord_y[gene_mask] - coord_attractor[2]
    d <- sqrt(x_diff^2 + y_diff^2)
    # Normalize distance
    d_norm <- 1 - d/max(d)
    # Add noise to normalization by using it as mean for beta distribution
    shape1 <- 1
    shape2 <- (shape1 - (shape1 * d_norm)) / d_norm
    n_points <- sum(gene_mask)
    d_norm_noise <- rbeta(n_points, shape1, shape2)
    # Scale differences with noisy normalized distance
    x_diff_norm <- x_diff * d_norm_noise
    y_diff_norm <- y_diff * d_norm_noise
    # Apply attraction
    data_attracted <- data_mixed
    data_attracted$coord_x[gene_mask] <- data_mixed$coord_x[gene_mask] - x_diff_norm 
    data_attracted$coord_y[gene_mask] <- data_mixed$coord_y[gene_mask] - y_diff_norm 
    # Apply effect
    data_attracted$count[gene_mask] <- rpois(sum(gene_mask), data_attracted$count[gene_mask] * (2*effect)^2)
    return(data_attracted)
  }

# Helper function to bin data
bin_data <- function(
    data,
    n_bins = 100,
    sum_counts = FALSE
  ) { 
    
    # Bin coordinates
    data$bin_x <- cut(
      data$coord_x,
      breaks = seq(min(data$coord_x), max(data$coord_x), length.out = n_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
    data$bin_y <- cut(
      data$coord_y,
      breaks = seq(min(data$coord_y), max(data$coord_y), length.out = n_bins + 1),
      include.lowest = TRUE,
      labels = FALSE
    )
    
    if (sum_counts) {
      # Get number of genes 
      genes <- unique(data$gene)
      n_genes <- length(genes)
      
      # Sum counts in bins by genes 
      countx <- rep(0, n_bins * n_genes)
      county <- rep(0, n_bins * n_genes)
      for (gi in c(1:n_genes)) {
        g <- genes[gi]
        maskg <- data$gene == g
        for (i in c(1:n_bins)) {
          maskx <- data$bin_x == i & maskg
          masky <- data$bin_y == i & maskg
          idx <- (gi - 1) * n_bins + i
          countx[idx] <- sum(data$count[maskx], na.rm = TRUE)
          county[idx] <- sum(data$count[masky], na.rm = TRUE)
        }
      }
      
      return(
        data.frame(
          bin = rep(c(1:n_bins), n_genes),
          countx = countx,
          county = county,
          gene = rep(genes, each = n_bins)
        )
      )
    } else {
      return(data)
    }
    
  }

# Function to simulate data 
simulate_data <- function(
    seed_data, 
    n_bins = 100,
    n_replicates = 4, # number of replicates per treatment condition 
    replicate_spatial_scalar = 0.05, 
    sum_counts = FALSE,
    print_plots = FALSE
  ) {
    
    # Print seed data 
    if (print_plots) {
      plt <- ggplot(seed_data, aes(x = coord_x, y = coord_y, color = gene, size = log(count + 1))) +
        geom_point() + 
        ggtitle("Seed data") +
        theme_minimal()
      print(plt)
    }
    
    # Get number and list of genes
    genes <- sort(unique(seed_data$gene))
    n_genes <- length(genes)
    
    # Randomly select genes to be spatially variable 
    # ... will select at least one
    SVGs <- sort(sample(n_genes, sample(n_genes, 1)))
    n_SVGs <- length(SVGs)
    
    # Randomly select SVGs to have FSEs
    # ... may select none
    FSEs <- c()
    n_FSEs <- sample(c(0:n_SVGs), 1)
    if (n_FSEs > 0) FSEs <- sort(sample(SVGs, n_FSEs))
    
    # Smooth seed data (no spatial variation)
    seed_data_smoothed <- smooth_data(seed_data)
    
    # Initialize variables
    sim_data_ref <- seed_data_smoothed
    sim_data_trt <- seed_data_smoothed
    
    # Select attractors
    attractor <- runif(n_genes)
    
    # Select effects
    effect <- rep(0.5, n_genes)
    if (n_FSEs > 0) effect[FSEs] <- runif(n_FSEs)
    
    # Smooth, apply attractors, apply effects
    for (g in c(1:n_genes)) {
      # Print smoothed gene
      if (print_plots) {
        mask <- sim_data_ref$gene == genes[g]
        sim_data_plt <- rbind(
          sim_data_ref[mask,],
          sim_data_trt[mask,]
        )
        sim_data_plt$fixedeffect <- c(rep("ref", sum(mask)), rep("trt", sum(mask)))
        plt <- ggplot(sim_data_plt, aes(x = coord_x, y = coord_y, color = fixedeffect, size = log(count + 1))) +
          geom_point() + 
          ggtitle(paste0(genes[g], ", smoothed")) +
          theme_minimal()
        print(plt)
      }
      if (any(g == SVGs)) {
        # Use attractor to induce spatial variability
        sim_data_ref <- induce_SV(sim_data_ref, genes[g], attractor[g])
        # Apply effects
        sim_data_trt <- induce_SV(sim_data_trt, genes[g], attractor[g], effect[g])
        # Print spatially variable gene
        if (print_plots) {
          mask <- sim_data_ref$gene == genes[g]
          sim_data_plt <- rbind(
            sim_data_ref[mask,],
            sim_data_trt[mask,]
          )
          sim_data_plt$fixedeffect <- c(rep("ref", sum(mask)), rep("trt", sum(mask)))
          plt <- ggplot(sim_data_plt, aes(x = coord_x, y = coord_y, color = fixedeffect, size = log(count + 1))) +
            geom_point() + 
            ggtitle(paste0(genes[g], ", spatially variable")) +
            theme_minimal()
          print(plt)
        }
      }
    }
    
    # Select replicate rate scalars 
    replicate_rate_scalars <- runif(n_replicates)
    
    # Make replicates and bin data
    rep_names <- paste0("replicate", c(1:n_replicates))
    nrow_rep <- nrow(sim_data_ref)
    ncol_rep <- ncol(sim_data_ref) + 2
    if (sum_counts) {
      nrow_rep <- n_bins * n_genes
      ncol_rep <- 4
    }
    sim_data_ref_reps <- as.data.frame(matrix(NA, nrow = n_replicates * nrow_rep, ncol = ncol_rep))
    sim_data_trt_reps <- as.data.frame(matrix(NA, nrow = n_replicates * nrow_rep, ncol = ncol_rep))
    for (r in c(1:n_replicates)) {
      # Make affine transform for replicate 
      Ax <- matrix(c(1,0,rnorm(1,0,replicate_spatial_scalar),1),2,2)  # shear x
      Ay <- matrix(c(1,rnorm(1,0,replicate_spatial_scalar),0,1),2,2)  # shear y
      As <- matrix(c(rnorm(1,0,replicate_spatial_scalar)+1, 0, 0, rnorm(1,0,replicate_spatial_scalar)+1),2,2) # scale
      A <- As %*% Ay %*% Ax
      # Set index
      idx_start <- (r - 1) * nrow_rep + 1
      idx_end <- r * nrow_rep
      idx <- c(idx_start:idx_end)
      # Make replicate and bin data
      sim_data_ref_reps[idx,] <- bin_data(make_replicate(data = sim_data_ref, rate_scalar = replicate_rate_scalars[r], affine_transform = A), n_bins, sum_counts)
      sim_data_trt_reps[idx,] <- bin_data(make_replicate(data = sim_data_trt, rate_scalar = replicate_rate_scalars[r], affine_transform = A), n_bins, sum_counts)
    }
    if (sum_counts) {
      colnames(sim_data_ref_reps) <- c("bin", "countx", "county", "gene")
      colnames(sim_data_trt_reps) <- c("bin", "countx", "county", "gene")
    } else {
      colnames(sim_data_ref_reps) <- c(colnames(sim_data_ref), "bin_x", "bin_y") 
      colnames(sim_data_trt_reps) <- c(colnames(sim_data_trt), "bin_x", "bin_y") 
    }
    sim_data_ref_reps$replicate <- rep(rep_names, each = nrow_rep)
    sim_data_trt_reps$replicate <- rep(rep_names, each = nrow_rep)
    sim_data_ref_reps$fixedeffect <- "ref"
    sim_data_trt_reps$fixedeffect <- "trt"
    
    # Combine binned data
    sim_data <- rbind(sim_data_ref_reps, sim_data_trt_reps)
    
    effect <- effect - 0.5
    names(effect) <- genes
    replicate_rate_scalars <- replicate_rate_scalars - 0.5
    names(replicate_rate_scalars) <- rep_names
    
    sim <- list(
      genes = genes,
      SVGs_idx = SVGs,
      SVGs = genes[SVGs],
      FSEs_idx = FSEs, 
      FSEs = genes[FSEs],
      attractor = attractor, 
      effect = effect, 
      replicate_rate_scalars = replicate_rate_scalars,
      data = sim_data
    )
    
    return(sim)
    
  }

# Function to extract ground-truth from sim data 
ground_truth <- function(sim) {
    fse_true <- sim$effect
    ran_true <- sim$replicate_rate_scalars
    SVGs <- rep(FALSE, length(sim$genes))
    names(SVGs) <- sim$genes
    SVGs[sim$SVGs] <- TRUE
    FSEs <- rep(FALSE, length(sim$genes))
    names(FSEs) <- sim$genes
    FSEs[sim$FSEs] <- TRUE 
    return(
      list(
        fse_true = fse_true,
        ran_true = ran_true,
        SVGs = SVGs,
        FSEs = FSEs
      )
    )
  }

# wisp benchmarking functions ##########################################################################################

# Rate effect is consistent across all of space, so, average over all blocks for each gene
extract_rate_effect_wisp <- function(model) {
  genes <- model[["grouping.variables"]][["species.lvls"]]
  rate_effect <- rep(NA, length(genes))
  names(rate_effect) <- genes
  for (g in genes) {
    mask <- grepl(paste0("beta_Rt_context_", g, "_trt_X"), names(model[["fitted.parameters"]]))
    rate_effect[g] <- mean(model[["fitted.parameters"]][mask])
  }
  return(rate_effect)
}

# Function to extract mean p-value across blocks for each gene 
extract_pvalue_wisp <- function(model) {
  genes <- model[["grouping.variables"]][["species.lvls"]]
  p_values <- rep(NA, length(genes))
  names(p_values) <- genes
  for (g in genes) {
    mask <- grepl(paste0("beta_Rt_context_", g, "_trt_X"), model[["stats"]][["parameters"]]$parameter)
    p_values[g] <- mean(model[["stats"]][["parameters"]]$p.value.adj[mask])
  }
  return(p_values)
}

# Replicate random noise is consistent across genes, so, average over all genes for each replicate
extract_random_effect_wisp <- function(model) {
  replicates <- model[["grouping.variables"]][["ran.lvls"]]
  random_effect <- rep(NA, length(replicates))
  names(random_effect) <- replicates
  for (r in replicates) {
    mask <- grepl(paste0("wfactor_rate_", r), names(model[["fitted.parameters"]]))
    random_effect[r] <- mean(model[["fitted.parameters"]][mask])
  }
  random_effect <- random_effect[-c(1)]
  return(random_effect)
} 

# Function to model simulation with wisp
model_sim_wisp <- function(
    sim, 
    sim_num, 
    fit_only = FALSE
  ) {
    
    # Transform data into count frame for wisp
    keep_cols <- !(colnames(sim$data) %in% c("slice_num", "coord_x", "coord_y", "bin_y"))
    sim_count <- sim$data[,keep_cols]
    
    # Set data variables
    data.variables <- list(
      bin = "bin_x", 
      count = "count",
      species = "gene",
      ran = "replicate"
    )
    
    # Fit model
    model <- wisp(
      count.data = sim_count,
      variables = data.variables,
      fit_only = fit_only,
      MCMC.settings = list(MCMC.steps = 0),
      bootstraps.num = 1e3,
      max.fork = bs_chunksize,
      plot.settings = list(print.plots = FALSE),
      verbose = FALSE,
      ran.seed = sim_num
    )
    
    # Extract model results for comparing to ground truth
    fse_est <- extract_rate_effect_wisp(model)
    ran_est <- extract_random_effect_wisp(model)
    fse_pvalues <- rep(NA, length(fse_est))
    names(fse_pvalues) <- names(fse_est)
    if (!fit_only) {
      fse_pvalues <- extract_pvalue_wisp(model)
    }
    
    # Extract ground-truth
    GT <- ground_truth(sim)
    fse_true <- GT$fse_true
    ran_true <- GT$ran_true
    FSEs <- GT$FSEs
    
    # Compile results in named vector 
    results <- data.frame(
      est = c(fse_est, ran_est, fse_pvalues),
      true = c(fse_true, ran_true, FSEs),
      param = c(
        rep("rate_effect", length(fse_est)), 
        rep("random_effect", length(ran_est)), 
        rep("FSE", length(fse_pvalues))
      ),
      id = c(
        names(fse_est), 
        names(ran_est), 
        names(fse_pvalues)
      )
    )
    
    # Add method and sim number
    results$method <- "wisp"
    results$sim <- sim_num
    
    return(results)
    
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
# /Users/michaelbarkasi/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade pip
# cd /Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/_molecular_mechanisms_of_ACx_lateralization/paper_WSPmm_method/wspmm_methods/ELLA
# /Users/michaelbarkasi/.virtualenvs/r-reticulate/bin/python -m pip install .

# For debugging
#options(reticulate.python.stdout = TRUE)
#options(reticulate.python.stderr = TRUE)

# import module
ELLA_mod <- import_from_path("ELLA", path = "./ELLA/ELLA")

# Make function to convert simulation data into ELLA format
convert_sim_data_to_ELLA <- function(
    sim_data,
    radial_dim = "bin_x",
    theta = "bin_y"
  ) {
    
    # ELLA: limit to "ref" fixedeffect, throw out "trt". replicate column becomes "cell". fixedeffect column becomes "type". 
    # bin_x and bin_y become x and y. centerX and centerY are the mean of bin_x and bin_y, for all rows. 
    # count becomes "umi". Compute sc_total. Throw out slice_num, coord_x, and coord_y.
    
    # Initialize top-level list structure
    ella_data <- list()
    length(ella_data) <- 7
    names(ella_data) <- c(
      "types",        # character vector of cell types
      "cells",        # named list, each element is a character vector of cell IDs for each type (with the type as the list element name)
      "cells_all",    # character vector of all cell IDsh
      "genes",        # named list, each element is a character vector of gene names for each type (with the type as the list element name)
      "cell_seg",     # data frame with columns x, y, and cell with points giving the cell boundary segmentation
      "nucleus_seg",  # data frame with columns x, y, and cell with points giving the nucleus boundary segmentation
      "expr"          # data frome with columns "type", "cell", "gene", "x", "y", "umi", "centerX", "centerY"  "sc_total"
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
      max_iter = 1000L
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
model_sim_ELLA <- function(
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
    GT <- ground_truth(sim)
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
      )
    )
    
    # Add method and sim number
    results$method <- "ELLA"
    results$sim <- sim_num
    
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

# CSIDE benchmarking functions #########################################################################################

# Load library 
# options(timeout = 600000000) ### set this to avoid timeout error
# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)

# Function to convert data to CSIDE format and run RCTD
convert_sim_data_to_CSIDE <- function(
    sim,
    ref_counts,
    max_cores = 2
  ) {
    
    # Need:
    #  DATA:
    #   count matrix (one per "puck"): rows as genes, columns as spots (i.e., spatial locations, or "pixels")
    #   coords matrix: rows as spots, columns as x and y coordinates
    #   nUMI vector: total counts per spot
    #  REFERENCE:
    #   ref_counts: matrix with rows as genes, columns as cells, pseudo-bulk library
    #   cell_types: named vector with cell type for each cell in ref_counts
    #   nUMI vector: total counts per cell in ref_counts
    
    # Grab data and split replicates by fixed effect
    sim_data <- sim$data
    sim_data$replicate <- paste(sim_data$replicate, sim_data$fixedeffect, sep = "_")
    replicates <- sort(unique(sim_data$replicate))
    grounp_ids <- as.integer(grepl("_trt", replicates)) + 1 # 1 for ref, 2 for trt
    
    # full grid of pixels
    rx <- range(sim_data$bin_x)
    ry <- range(sim_data$bin_y)
    pixels <- CJ(
      bin_x = rx[1]:rx[2],
      bin_y = ry[1]:ry[2]
    )
    pixels <- pixels[, pixel := paste(bin_x, bin_y, sep = "_")]
    pixels <- as.matrix(pixels)
    
    # Convert pixel matrix into coords matrix 
    coords <- pixels 
    rownames(coords) <- coords[, "pixel"]
    coords <- as.data.frame(coords[, c("bin_x", "bin_y")])
    colnames(coords) <- c("x", "y")
    coords$x <- as.numeric(coords$x)
    coords$y <- as.numeric(coords$y)
    
    # Create cell types for REFERENCE
    cell_types <- data.frame(
      cell = colnames(ref_counts),
      type = sample(c("celltype1", "celltype2"), ncol(ref_counts), replace = TRUE)
    )
    cell_types <- setNames(cell_types[[2]], cell_types[[1]])
    cell_types <- as.factor(cell_types) # convert to factor data type
    
    # Randomly create ten filler genes that are differentially expressed between cell types
    n_filler_genes <- 20
    n_original_genes <- length(unique(sim_data$gene))
    ct_mask <- cell_types == "celltype1"
    n_ct1_ref <- sum(ct_mask)
    n_ct2_ref <- sum(!ct_mask)
    this_slice_num <- sim_data$slice_num[1]
    filler_df <- data.frame()
    original_gene_names <- rownames(ref_counts)
    for (fg in c(1:n_filler_genes)) {
      
      # Select a gene from ref_counts
      fg_name <- sample(original_gene_names, 1)
      ct1_counts <- ref_counts[fg_name, ct_mask]
      ct1_counts <- ct1_counts * runif(1) * 2
      ct2_counts <- ct1_counts * (2 * runif(1))^2
      fg_name <- paste0("_filler_", fg, "_", fg_name)
      
      # Update ref_counts with filler gene
      new_row <- rep(0, ncol(ref_counts))
      new_row[ct_mask] <- rpois(n_ct1_ref, ct1_counts + 1)
      new_row[!ct_mask] <- rpois(n_ct2_ref, ct2_counts + 1)
      ref_counts <- rbind(ref_counts, new_row)
      rownames(ref_counts)[nrow(ref_counts)] <- fg_name
      
      for (r in replicates) {
        replicate_num <- as.integer(sub(".*?(\\d+).*", "\\1", r))
        rep_rate_scalar <- sim$replicate_rate_scalars[replicate_num] + 0.5
        r_mask <- sim_data$replicate == r
        n_spots <- sum(r_mask)/n_original_genes
        r_idx <- sample(which(r_mask), n_spots, replace = TRUE) 
        r_df <- sim_data[r_idx,]
        ct_mask_r_df <- r_df$bin_x >= max(r_df$bin_x)/2 # arbitrarily assign left half to celltype1, right half to celltype2
        n_ct1 <- sum(ct_mask_r_df)
        n_ct2 <- n_spots - n_ct1
        r_df$count[ct_mask_r_df] <- rpois(n_ct1, sample(ct1_counts * 2 * rep_rate_scalar, n_ct1, replace = TRUE) + 1)
        r_df$count[!ct_mask_r_df] <- rpois(n_ct2, sample(ct2_counts * 2 * rep_rate_scalar, n_ct2, replace = TRUE) + 1)
        r_df$gene <- fg_name
        filler_df <- rbind(filler_df, r_df)
      }
      
    }
    sim_data <- rbind(sim_data, filler_df)
    
    # Make reference library counts 
    nUMI_ref <- colSums(ref_counts) 
    names(nUMI_ref) <- colnames(ref_counts)
    
    # Make reference
    reference <- Reference(ref_counts, cell_types, nUMI_ref)
    
    # Make data for each replicate
    n_replicates <- length(replicates)
    pucks <- list()
    length(pucks) <- n_replicates
    names(pucks) <- replicates
    barcodes <- list()
    length(barcodes) <- n_replicates
    names(barcodes) <- replicates
    for (r in replicates) {
      
      # Prune down data to just this replicate and necessary columns
      mask_r <- sim_data$replicate == r
      data_thin <- sim_data[mask_r,c("gene", "count", "bin_x", "bin_y")]
      # Count matrix 
      dt <- as.data.table(data_thin)
      # aggregate counts
      dt[, pixel := paste(bin_x, bin_y, sep = "_")]
      dt <- dt[, .(count = sum(count)), by = .(gene, pixel)]
      # ensure all gene–pixel combinations exist
      dt_full <- merge(
        CJ(gene = unique(dt$gene), pixel = pixels[,"pixel"]),
        dt,
        by = c("gene", "pixel"),
        all.x = TRUE
      )
      dt_full[is.na(count), count := 0]
      # cast to matrix
      pixel_mat <- dcast(
        dt_full,
        gene ~ pixel,
        value.var = "count",
        fill = 0
      )
      pixel_mat <- as.matrix(pixel_mat[, -1, with = FALSE])
      rownames(pixel_mat) <- dt_full[, unique(gene)]
      
      # Get library counts 
      nUMI_spatial <- colSums(pixel_mat)
      
      # Make and save puck
      pucks[[r]] <- SpatialRNA(coords, pixel_mat, nUMI_spatial)
      barcodes[[r]] <- colnames(pucks[[r]]@counts)
      
    }
    
    barcodes <- barcodes[[1]]
    
    # Create RCTD object
    low_panel_cutoff <- 0 # ... set to zero to account for small panel
    myRCTD.reps <- create.RCTD.replicates(
      spatialRNA.replicates = pucks,
      reference = reference,
      replicate_names = replicates,
      group_ids = grounp_ids,
      max_cores = max_cores,
      gene_cutoff = low_panel_cutoff, 
      fc_cutoff = low_panel_cutoff,
      gene_cutoff_reg = low_panel_cutoff, 
      fc_cutoff_reg = low_panel_cutoff,
      # UMI_min = low_panel_cutoff,
      # UMI_min_sigma = low_panel_cutoff,
      CELL_MIN_INSTANCE = 1,
      CONFIDENCE_THRESHOLD = 1
    )
    myRCTD.reps <- run.RCTD.replicates(myRCTD.reps)
    
    return(
      list(
        RCTD = myRCTD.reps,
        barcodes = barcodes
      )
    )
    
  }

# Function to run CSIDE
run_CSIDE <- function(
    sim,
    ref_counts,
    max_cores = 2
  ) {
    
    # Convert data
    RCTD <- tryCatch(
      convert_sim_data_to_CSIDE(
        sim = sim,
        ref_counts = ref_counts,
        max_cores = max_cores
      ),
      error = function(e) {
        cat("\nError occurred with CSIDE data conversion or RCTD. Retrying with different random draws...")
        convert_sim_data_to_CSIDE(
          sim = sim,
          ref_counts = ref_counts,
          max_cores = max_cores
        )
      }
    )
    
    # Run CSIDE
    myRCTD.reps <- run.CSIDE.replicates(
      RCTD.replicates = RCTD$RCTD,
      cell_types = c("celltype1", "celltype2"),
      cell_type_threshold = 5, 
      weight_threshold = 0,
      fdr = 0.05,
      population_de = TRUE,
      de_mode = "nonparam",
      barcodes = RCTD$barcodes,
      log_fc_thresh = 0
    )
    
    # Run population inference with groups
    myRCTD.pop <- CSIDE.population.inference(
      myRCTD.reps, 
      fdr = 0.05, 
      use.groups = TRUE,
      MIN.CONV.REPLICATES = 1,
      MIN.CONV.GROUPS = 1,
      CT.PROP = 0,
      log_fc_thresh = 0
      )
    
    # Return completed model
    return(
      list(
        myRCTD.reps = myRCTD.reps,
        myRCTD.pop = myRCTD.pop
      )
    )
    
  }

test <- run_CSIDE(
  sim = sim_data,
  ref_counts = ref_counts,
  max_cores = bs_chunksize
)

# Full model pipeline
model_sim_CSIDE <- function(
    sim,
    sim_num,
    ref_counts,
    max_cores
  ) {
    
    # Run CSIDE
    cside_model <- run_CSIDE(
      sim = sim,
      ref_counts = ref_counts,
      max_cores = max_cores
    )
    
    # Extract model results for comparing to ground truth
    de_results <- cside_model@population_de_results
    # ... use just celltype1, discard celltype2 results (both cell types were simulated identically)
    de_results <- de_results$celltype1
    # ... filter down to just the genes from the simulation (discard filler genes)
    de_results <- de_results[rownames(de_results) %in% unique(sim$genes), ]
    # ... extract estimated replicate variance per gene 
    # ... extract log fold change per gene (FSE)
    fse_est <- de_results$log_fc_est
    ran_est
    fse_pvalues
    svg_pvalues
    
    # Extract and infer estimated random effects
    n_replicates <- length(cside_model@RCTD.reps)
    for (r in c(1:n_replicates)) {
      # Grab estimated expression per gene for this replicate for celltype1
      gene_exp_r <- cside_model@RCTD.reps[[r]]@de_results[["gene_fits"]][["mean_val"]][,"celltype1"]
    }
    
    
    # Extract model results for comparing to ground truth
    #fse_est <- extract_rate_effect_wisp(model)
    #ran_est <- extract_random_effect_wisp(model)
    #fse_pvalues <- rep(NA, length(fse_est))
    #names(fse_pvalues) <- names(fse_est)
    #if (!fit_only) {
    #  fse_pvalues <- extract_pvalue_wisp(model)
    #}
    
    # Extract ground-truth
    GT <- ground_truth(sim)
    fse_true <- GT$fse_true
    ran_true <- GT$ran_true
    FSEs <- GT$FSEs
    
    # Compile results in named vector 
    results <- data.frame(
      est = c(fse_est, ran_est, fse_pvalues),
      true = c(fse_true, ran_true, FSEs),
      param = c(
        rep("rate_effect", length(fse_est)), 
        rep("random_effect", length(ran_est)), 
        rep("FSE", length(fse_pvalues))
      ),
      id = c(
        names(fse_est), 
        names(ran_est), 
        names(fse_pvalues)
      )
    )
    
    # Add method and sim number
    results$method <- "wisp"
    results$sim <- sim_num
    
    return(results)
  
  }



# Run benchmarking #####################################################################################################

n_sims <- 100
results <- data.frame()
for (s in c(1:n_sims)) {
  
  d_sim <- rep(NA, n_sims)
  d_wisp <- rep(NA, n_sims)
  d_ella <- rep(NA, n_sims)
  
  cat("\n\n----------------------")
  cat("\nRunning simulation ", s, "/", n_sims)
  
  # Simulate data and extract count data for modeling
  t_stim_start <- Sys.time()
  sim_data <- simulate_data(count_data_neurons_patch)
  dsim <- Sys.time() - t_stim_start
  units(dsim) <- "secs"
  d_sim[s] <- dsim
  cat("\n\tdata sim time:", round(d_sim[s],3), "s") 
  
  # Model simulated data with wisp
  t_wisp_start <- Sys.time()
  wisp_results <- model_sim_wisp(sim_data, sim_num = s)
  dwisp <- Sys.time() - t_wisp_start
  units(dwisp) <- "secs"
  d_wisp[s] <- dwisp
  results <- rbind(results, wisp_results)
  cat("\n\twisp time:", round(d_wisp[s],3), "s") 
  
  # Model simulated data with ELLA
  cat("\n\n") # ... skip line for readability
  t_ella_start <- Sys.time()
  ella_results <- model_sim_ELLA(sim_data, sim_num = s)
  della <- Sys.time() - t_ella_start
  units(della) <- "secs"
  d_ella[s] <- della
  results <- rbind(results, ella_results)
  cat("\n\tELLA time:", round(d_ella[s],3), "s")
  
}

# Save results
write.csv(results, "benchmark_results.csv", row.names = FALSE)

# Load results 
results <- read.csv("benchmark_results.csv")

analyze_results <- function(
    results
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
        FP <- sum(results[r_m_mask & results$true == FALSE, "est"] < 0.05)
        TN <- sum(results[r_m_mask & results$true == FALSE, "est"] >= 0.05)
        TP <- sum(results[r_m_mask & results$true == TRUE, "est"] < 0.05)
        FN <- sum(results[r_m_mask & results$true == TRUE, "est"] >= 0.05)
        summary["FPR_svg", m] <- FP / (FP + TN)
        summary["FDR_svg", m] <- FP / (FP + TP)
        summary["power_svg", m] <- TP / (TP + FN)
      }
      # Compute FPR, FDR, power for FSE parameter
      if (any("FSE" == results$param[m_mask])) {
        r_m_mask <- m_mask & results$param == "FSE"
        FP <- sum(results[r_m_mask & results$true == FALSE, "est"] < 0.05)
        TN <- sum(results[r_m_mask & results$true == FALSE, "est"] >= 0.05)
        TP <- sum(results[r_m_mask & results$true == TRUE, "est"] < 0.05)
        FN <- sum(results[r_m_mask & results$true == TRUE, "est"] >= 0.05)
        summary["FPR_fse", m] <- FP / (FP + TN)
        summary["FDR_fse", m] <- FP / (FP + TP)
        summary["power_fse", m] <- TP / (TP + FN)
      }
      
    }
    return(summary)
  }
results_summary <- analyze_results(results)
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










