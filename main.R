
# Code for the paper "Logistic regression for estimating functional effects with spatial transcriptomics"
# Preprint: https://doi.org/10.1101/2025.06.11.659209

# Setup ################################################################################################################

# Clear global environment and set project folder
rm(list = ls())
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

# For for WSPmm and snk ("sink") printing
# ... version 1.0 used for the preprint (draft 1)
# ... version 2.0 used for the resubmission to NAR (draft 2)
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
snk.report("Code for the paper 'Logistic regression for estimating functional effects with spatial transcriptomics'", end_breaks = 1)

# Set file path, and bootstrap chunk size
data_path <- paste0(projects_folder, "_molecular_mechanisms_of_ACx_lateralization/data_SSp/")
if (Sys.info()[["sysname"]] == "Linux") { # ... for workstation: 
  data_path <- "/home/oviedoworkstation/Desktop/SharedStorageMount/SSD_Main/MERFISH/Mike/data_SSp/"
} 
bs_chunksize <- 20

# Preprocessing MERFISH data ###########################################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Preprocessing MERFISH data for RORB analysis", end_breaks = 1)

# Define list of genes to analyze
laminar.gene.list <- c("Bcl11b", "Fezf2", "Satb2", "Nxph3", "Cux2", "Rorb")  

# Load and parse raw visgen data from files
count_data <- make_count_data(
    data_path,
    remove_L1 = TRUE,        # If TRUE, removes any points labeled as layer 1
    ROIname = "Primary somatosensory area",
    raw = TRUE               # If FALSE, uses normalized data (that is not desired for this paper)
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

# Label layer boundary bins 
layer_boundaries <- colMeans(layer.boundary.bins)[-c(1, ncol(layer.boundary.bins))]
names(layer_boundaries) <- c("L2/3", "L4", "L5", "L6")

# Simplify cell types and limit to just neurons 
snk.report("Simplifying to just Glut and GABA cells", end_breaks = 1)
# ... Note: types come from naive application of MapMyCells to counts from full gene panel
count_data$celltype_MMC <- as.character(count_data$celltype_MMC)
# ... Make masks
Glut_mask <- grepl("Glut", count_data$celltype_MMC)
GABA_mask <- grepl("GABA", count_data$celltype_MMC)
# ... Apply new annotations
count_data$celltype_MMC[Glut_mask] <- "Glut"
count_data$celltype_MMC[GABA_mask] <- "GABA"
# ... Remove other cell types
count_data <- count_data[Glut_mask | GABA_mask, ]
# ... Rename column 
colnames(count_data)[which(colnames(count_data) == "celltype_MMC")] <- "celltype"

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
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Running wisp analysis of RORB MERFISH data", end_breaks = 1)

# Define fixed effects to test
fixed.effect.names <- c("hemisphere", "age")

# Create count data for WSPmm object, from preprocessed count_data, using laminar axis (y)
count.data.WSPmm <- create.count.data.WSPmm(
    df.merfish = count_data,
    bin.dim = "y_bins",
    gene.list = laminar.gene.list,
    fixed.effect.names = fixed.effect.names, 
    context = "celltype"
  )

# Make genes upper case 
count.data.WSPmm$gene <- toupper(count.data.WSPmm$gene)

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
    context = "celltype", 
    species = "gene",
    ran = "mouse",
    timeseries = "age",
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
    MCMC.step.size = 0.1,
    MCMC.prior = 1.0, 
    MCMC.neighbor.filter = 2
  )

# Fit model
laminar.model <- wisp(
    count.data = count.data.WSPmm,
    variables = data.variables,
    use.median = FALSE,
    bootstraps.num = 1e4,
    converged.resamples.only = TRUE,
    max.fork = bs_chunksize,
    verbose = TRUE,
    plot.settings = list(
      print.plots = TRUE, 
      dim.bounds = layer_boundaries
    ),
    MCMC.settings = MCMC.settings,
    model.settings = model.settings
  )

# Save
# ... load with: laminar.model <- readRDS("saved_laminar_model.rds")
n <- length(laminar.model)
k <- 14
saveRDS(laminar.model[1:k], "saved_laminar_model_part1.rds")
saveRDS(laminar.model[(k+1):n], "saved_laminar_model_part2.rds")

# Make and export figures for MERFISH data #############################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Making and exporting figures for MERFISH data", end_breaks = 2)

load_saved_laminar_model <- FALSE 
if (load_saved_laminar_model) {
  laminar.model1 <- readRDS("saved_laminar_model_NARresub_feb14_part1.rds")
  laminar.model2 <- readRDS("saved_laminar_model_NARresub_feb14_part2.rds")
  laminar.model <- c(laminar.model1, laminar.model2)
  plots <- list(
    ratecount = plot.ratecount(laminar.model)
  )
  laminar.model[["plots"]] <- plots
}


plt_laminar_ratecount <- laminar.model$plots$ratecount[["plot_all"]] +
  theme(legend.position = "bottom") + 
  ggtitle("Cortical layer model: Rate by time") +
  labs(x = "Laminar position")

# Export rate-count plot
ggsave(
  "fig_laminar_ratecount.png", 
  plot = plt_laminar_ratecount, 
  width = 11, height = 8, dpi = 900
  )

# Make and export residuals figure (histogram and qq plot)
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

# Make and export MCMC vs. bootstrap comparison figure (walks, autocorrelation, normality)
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

# Make results table for stats 
snk.report("Making stat results table for MERFISH data", end_breaks = 2)

kableExtra::kbl(
  laminar.model[["stats"]][["parameters"]], 
  format = "latex", 
  booktabs = TRUE, 
  escape = FALSE, 
  caption = "Parameter estimates for MERFISH cortex model.\\label{table:FEestimatesCortex}", 
  linesep = "")

# Fit WSPmm model to radial data #######################################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Running wisp analysis of radial liver data", end_breaks = 1)

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
  snk.report("Importing pre-preprocessed Droin MatLab data", end_breaks = 1)
  # Import
  count_data <- read.csv("Droin_radial_count_data_sim.csv")
}

# Make gene column uppercase
count_data$gene <- toupper(count_data$gene)

# Define fixed effects to test
fixed.effect.names <- c("ZT")

# Define variables in the dataframe for the model
data.variables <- list(
  count = "count",
  bin = "bin", 
  context = "hepatocytes", 
  species = "gene",
  ran = "mouse",
  timeseries = "ZT",
  fixedeffects = fixed.effect.names
)

# Model settings 
model.settings <- list(
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

# Settings for plotting
plot.settings <- list(
  print.plots = FALSE,
  title_size = 20,
  count_size = 2.5,
  count_jitter = 0.0,
  count.alpha.none = 0.8,
  count.alpha.ran = 0.5,
  pred.alpha.ran = 0.0
)

# Fit model
radial.model <- wisp(
  count.data = count_data,
  variables = data.variables,
  bootstraps.num = 1e4,
  max.fork = bs_chunksize,
  model.settings = model.settings,
  MCMC.settings = list(MCMC.steps = 0),
  plot.settings = plot.settings
)

# Adjust and export figures for liver data
snk.report("Exporting figures for liver data", end_breaks = 2)

plt_radial_ratecount <- radial.model$plots$ratecount[["plot_all"]] +
  facet_wrap(~ species + context, scales = "free_y", ncol = 2) + 
  theme(legend.position = "bottom") + 
  ggtitle("Radial liver model: Rate by zone") +
  labs(x = "Lobule zone", fill = "ZT", color = "ZT") 
plt_radial_timeseries <- radial.model$plots$timeseries[["plot_all"]] + 
  facet_wrap(~ species + context, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = c(0, 6, 12, 18)) +
  labs(x = "Zeitgeber Time (ZT)") +
  theme(legend.position = "bottom") + 
  ggtitle("Radial liver model: Rate by time")

ggsave(
  "fig_radial_ratecount.png", 
  plot = plt_radial_ratecount, 
  width = 5.5, height = 8, dpi = 900
  )
ggsave(
  "fig_radial_timeseries.png", 
  plot = plt_radial_timeseries, 
  width = 5.5, height = 8, dpi = 900
  )

# Make results table for stats 
snk.report("Making stat results table for MERFISH data", end_breaks = 2)

kableExtra::kbl(
  radial.model[["stats"]][["parameters"]], 
  format = "latex", 
  booktabs = TRUE, 
  escape = FALSE, 
  caption = "Parameter estimates for radial liver model.\\label{table:FEestimatesLiver}", 
  linesep = "")

# Save
radial.model$Cpp_model <- NULL
radial.model$plots <- NULL
saveRDS(radial.model, "saved_radial_model.rds")

# Simulated data #######################################################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Making seed data for benchmark simulations", end_breaks = 1)

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
  snk.report("Importing pre-preprocessed Allen data", end_breaks = 1)
  # Import
  count_data_neurons <- read.csv("Allen_data.csv")
}

# Remove cell_id column from count_data_neurons 
count_data_neurons <- count_data_neurons[,c("count", "slice_num", "gene", "coord_x", "coord_y")]
count_data_neurons$gene <- toupper(count_data_neurons$gene)

# Cut down to just a single 2x2 patch from one slice
count_data_neurons_patch <- count_data_neurons[
  count_data_neurons$slice_num == 33 &
    count_data_neurons$coord_y >= 2 & count_data_neurons$coord_y <= 4 &
    count_data_neurons$coord_x >= 1 & count_data_neurons$coord_x <= 3,]

# Print demo data 
snk.print_table("Benchmark simulation seed patch data", count_data_neurons_patch)

# DESeq2 benchmarking functions ########################################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Defining DESeq2 benchmarking functions", end_breaks = 1)

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
    sim_num,
    use.adaptive.shrinkage = TRUE,
    multiplier = 1,     # multiplier for adaptive shrinkage; default = 1
    prior.scale = 1,    # analog of a ridge/lasso penalty magnitude; default = 1
    prior.df = 1        # smaller = strong strinkage of small effects, weak of large ones; default = 1
  ) {
    
    # Make data for DESeq2
    ddata <- make_DESeq2_data(sim$data)
    dds <- DESeqDataSetFromMatrix(countData = ddata$cts, colData = ddata$coldata, design = ~ fixedeffect)
    
    # Fit model
    # ... run differential expression analysis on fixedeffect
    dds <- DESeq(estimateSizeFactors(dds), fitType = "mean", minReplicatesForReplace = 3)
    # ... apply shrinkage (ridge regression)
    if (use.adaptive.shrinkage) {
      resLFC_fixedeffect <- lfcShrink(
        dds, 
        coef = "fixedeffect_trt_vs_ref", 
        type = "apeglm",
        apeAdapt = use.adaptive.shrinkage,
        multiplier = multiplier
      )
    } else {
      resLFC_fixedeffect <- lfcShrink(
        dds, 
        coef = "fixedeffect_trt_vs_ref", 
        type = "apeglm", 
        apeAdapt = use.adaptive.shrinkage,
        prior.control = list(
          no.shrink = 1,
          prior.mean = 0,
          prior.scale = prior.scale,    
          prior.df = prior.df,       
          prior.no.shrink.mean = 0,
          prior.no.shrink.scale = 15
        )
      )
    }
    
    # For debugging
    # print(resLFC_fixedeffect@priorInfo)
    
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
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Defining ELLA benchmarking functions", end_breaks = 1)

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
make_ELLA_data <- function(
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
    sim_data,
    L1.lam = 0.0
  ) { 
    
    # construct ELLA object
    ELLA_class <- ELLA_mod$ELLA
    ella_sim <- ELLA_class(
      dataset = "sim_benchmark",
      adam_learning_rate_min = 1e-2,
      max_iter = 1000L, 
      L1_lam = L1.lam
    )
    
    # convert simulated data to ELLA format
    ella_data <- r_to_py(make_ELLA_data(sim_data))
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

# Extract SVG p-values from ELLA results
extract_svg <- function(ella_sim) {
  pv_svg <- unlist(ella_sim$pv_fdr_tl[["ref"]])
  names(pv_svg) <- ella_sim$gene_list_dict[["ref"]]
  pv_svg
}

# Full model pipeline 
model_attractor_simulation_ELLA <- function(
    sim,
    sim_num,
    L1.lam = 0.0
  ) {
    
    # Run ELLA
    model <- run_ELLA(sim$data, L1.lam)
    
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
    
    text_scalar <- 1.25
    
    # filter
    df_rad <- df_rad[df_rad$cell == replicate,]
    df_rad$count <- df_rad$umi
    df_car <- df_car[df_car$replicate == replicate,]
    df_car$x <- df_car$bin_x
    df_car$y <- df_car$bin_y
    
    # scatter plots 
    p1 <- ggplot(
        test_sim$data[test_sim$data$replicate == "replicate1" & test_sim$data$fixedeffect == "ref",], 
        aes(x = bin_x, y = bin_y, color = gene, size = count)
      ) +
      geom_point(alpha = 0.25) +
      ggtitle("Simulation, native Cartesian coordinates") +
      labs(x = "x bin", y = "y bin") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = laminar.model$plot.settings$title_size),
        axis.title = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        axis.text = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        legend.title = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        legend.text = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position = "bottom"
      ) +
      guides(size = "none")
    p2 <- ggplot(ella_sim$expr[ella_sim$expr$cell == "replicate1",], aes(x = x, y = y, color = gene, size = umi)) +
      geom_point(alpha = 0.5) +
      ggtitle("Simulation, transformed radial coordinates") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = laminar.model$plot.settings$title_size),
        axis.title = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        axis.text = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        legend.title = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        legend.text = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position = "bottom"
      ) +
      guides(size = "none")
    
    # density plots 
    eps <- 4e-3 # For matching the log-scale transforms
    
    # density vs x 
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
      labs(x = "x bin", y = "count density (log10)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = text_scalar * laminar.model$plot.settings$title_size),
        axis.title = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        axis.text = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        legend.title = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        legend.text = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position = "bottom"
      )
    
    ## density vs radius
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
      radial_density =
        as.vector(rcount) /
        (2 * pi * rep(rmid, each = nrow(rcount)) * rep(rwidth, each = nrow(rcount)))
    )
    p_r <- ggplot(df_r, aes(x = rmid, y = radial_density + eps, color = gene)) +
      geom_line() +
      scale_y_log10() +
      labs(x = "radius from origin", y = "radial density (log10)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = text_scalar * laminar.model$plot.settings$title_size),
        axis.title = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        axis.text = element_text(size = text_scalar * laminar.model$plot.settings$axis_size),
        legend.title = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        legend.text = element_text(size = text_scalar * laminar.model$plot.settings$legend_size),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        legend.position = "bottom"
      )
    
    # arrange
    title <- grobTree(
      rectGrob(gp = gpar(fill = "white", col = NA)),
      textGrob("Resulting simulation", gp = gpar(fontsize = 1.25 * text_scalar * laminar.model$plot.settings$title_size))
    )
    panels <- arrangeGrob(p_x, p_r, ncol = 2)
    panels_scatter <- arrangeGrob(p1, p2, ncol = 2)
    
    return(
      grid.arrange(
        title,
        panels_scatter,
        panels,
        ncol = 1,
        heights = c(0.1, 1, 1)
      )
    )
    
  }

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

# Make simulation demo plots ###########################################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Making and exporting demo plots of attractor simulations", end_breaks = 1)

# Plot seed data
plt_sim_seed <- ggplot(count_data_neurons[count_data_neurons$slice_num == 33,], aes(x = coord_x, y = coord_y, color = gene)) +
  geom_point(size = 0.1) + 
  ggtitle("Seed data: Coronal slice from Allen MERFISH data") +
  labs(x = "x coordinate", y = "y coordinate") +
  geom_rect(
    xmin = 1, xmax = 3,
    ymin = 2, ymax = 4,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = laminar.model$plot.settings$title_size),
    axis.title = element_text(size = laminar.model$plot.settings$axis_size),
    axis.text = element_text(size = laminar.model$plot.settings$axis_size),
    legend.title = element_text(size = laminar.model$plot.settings$legend_size),
    legend.text = element_text(size = laminar.model$plot.settings$legend_size),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    legend.position = "bottom"
  )

ggsave(
  "fig_sim_demo_seed.png", 
  plot = plt_sim_seed, 
  width = 8, height = 8, dpi = 900
)

# Make sample simulation 
# ... Ensure reproducibility
set.seed(12399)
# ... Run simulation
test_sim <- attractor_simulation(
  seed_data = count_data_neurons_patch, 
  n_bins = 100,
  n_replicates = 4,
  replicate_spatial_scalar = 0.05, 
  min_effect_size = 0.05,
  print_plots = TRUE
)

# Plot Tac2 simulation
recolored_TAC2 <- list()
for (p in c(1:length(test_sim$plots[["TAC2"]]))) {
  recolored_TAC2[[p]] <- test_sim$plots[["TAC2"]][[p]]  + 
    labs(
      x = "x coordinate",
      y = "y coordinate"
    ) 
  recolored_TAC2[[p]]$layers[[1]]$aes_params$colour <- scales::hue_pal()(4)[3]
}
plt_TAC2 <- do.call(grid.arrange, c(recolored_TAC2, ncol = 2))

ggsave(
  "fig_sim_demo_effects.png", 
  plot = plt_TAC2, 
  width = 8, height = 8, dpi = 900
)

# Plot simulated replicates
mask <- test_sim$data$gene == "TAC2" & test_sim$data$fixedeffect == "ref"
data <- test_sim$data[mask,]
plt_TAC2_reps <- ggplot(data, aes(x = coord_x, y = coord_y, size = log(count + 1))) +
  geom_point(alpha = 0.5, color = scales::hue_pal()(4)[3]) + 
  facet_wrap(~ replicate) +
  ggtitle("TAC2 replicates under reference condition") +
  labs(x = "x coordinate", y = "y coordinate") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = laminar.model$plot.settings$title_size),
    axis.title = element_text(size = laminar.model$plot.settings$axis_size),
    axis.text = element_text(size = laminar.model$plot.settings$axis_size),
    legend.title = element_text(size = laminar.model$plot.settings$legend_size),
    legend.text = element_text(size = laminar.model$plot.settings$legend_size),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    strip.text = element_text(size = laminar.model$plot.settings$title_size/1.5)
    )

ggsave(
  "fig_sim_demo_replicates.png", 
  plot = plt_TAC2_reps, 
  width = 8, height = 8, dpi = 900
)

# Plot transformation for ELLA
ella_sim <- make_ELLA_data(test_sim$data)
plt_ELLA_transform <- plot_count_density(ella_sim$expr, test_sim$data)
ggsave(
  "fig_sim_demo_ELLA_transform.png", 
  plot = plt_ELLA_transform, 
  width = 12, height = 12, dpi = 900
)

# Run benchmarking #####################################################################################################
snk.horizontal_rule(initial_breaks = 2, end_breaks = 0)
snk.report("Running and analyzing benchmarks", end_breaks = 1)

results <- run_attractor_sim_benchmarks(
  seed_data = count_data_neurons_patch,
  n_sims = 250,
  modeling_functions = list(
    wisp = model_attractor_simulation_wisp, 
    ELLA = model_attractor_simulation_ELLA,
    DESeq2 = model_attractor_simulation_DESeq2
    ),
  modeling_function_args = list(
    wisp = list(bs_num = 1e3, max_fork = 100),
    ELLA = list(L1.lam = 0.2),
    DESeq2 = list(use.adaptive.shrinkage = TRUE, multiplier = 1, prior.scale = 1, prior.df = 1)
  )
)

# Save results
write.csv(results, "benchmark_results_temp.csv", row.names = FALSE)

# Load results 
#results <- read.csv("benchmark_results.csv")

# Analyze results
results_summary <- analyze_attractor_sim_benchmarks(results)
snk.print_table("Benchmark Results", results_summary, head = FALSE, initial_breaks = 1)

# Print session info
snk.report("Done! Printing session info", end_breaks = 2)
devtools::session_info()

# end sink
sink(file = NULL)










