
# Setup ################################################################################################################
# Analysis of MERFISH data with Warped Sigmoid, Poisson-Process Mixed-Effects Model (WSPmm)

# Clear global environment
rm(list = ls())
# If using VS Code with httpgd, ensure clean start
httpgd::hgd_close() 
projects_folder <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/"

# For for WSPmm 
wispack_path <- paste0(projects_folder, "R_packages/wispack/wispack_1.0.tar.gz")
#install.packages(wispack_path, repos = NULL)
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
  wf_initial = 0.5,                                     # initial value for wfactor (0.5 results in faster and better fits than 0.25 or 0.05)
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
  bootstraps.num = 1e2,
  converged.resamples.only = FALSE,
  max.fork = bs_chunksize,
  batch.size = bs_chunksize,
  dim.bounds = colMeans(layer_boundary_bins),
  verbose = TRUE,
  print.child.summaries = TRUE, 
  model.settings = model.settings
)



##############################
##############################

# Grab resulting data frame
count <- merfish_laminar_model[["count.data.summed"]]

# Grab random effects levels 
mice <- unique(count$ran[count$ran != "none"])
n_mice <- length(mice)

# Get treatment levels for each mouse
trts <- merfish_laminar_model$treatment$names
trts_comps <- merfish_laminar_model$treatment$components
names(trts_comps) <- trts
n_trts <- length(trts)
trt_lvls <- list()
for (m in mice) {
  mask <- count$ran == m & !is.na(count$count)
  trt_lvls[[m]] <- unique(count$treatment[mask])
}

# Grab fixed effects levels
fix <- merfish_laminar_model$fix
names(fix$lvls) <- fix$name
names(fix$treat.lvl) <- fix$name
names(fix$ref.lvl) <- fix$name

# Make treatment component array
trt_comp_array <- array(
  as.character(NA),
  dim = c(n_trts, 2),
  dimnames = list(trts, fix$name)
)
for (t in trts) {
  if (t == "ref") {
    trt_comp_array[t,] <- fix$ref.lvl
  } else {
    for (fe in fix$name) {
      if (any(trts_comps[[t]] %in% fix$lvls[[fe]])) {
        trt_comp_array[t, fe] <- fix$treat.lvl[[fe]]
      } else {
        trt_comp_array[t, fe] <- fix$ref.lvl[fe]
      }
    }
  }
}

# Make data matrices 
genes <- unique(count$child)
n_genes <- length(genes)
bins <- unique(count$bin)
n_bins <- length(bins)
data0 <- array(
  NA, 
  dim = c(
    n_bins, n_genes,
    n_mice,
    length(trt_lvls)
    ),
  dimnames = list(bin = bins, child = genes, ran = mice, fixedeffect = trts)
  )
for (m in mice) {
  for (t in trt_lvls[[m]]) {
    for (g in genes) {
      mask <- count$ran == m & count$treatment == t & count$child == g
      if (sum(mask) == n_bins) data0[, g, m, t] <- count$pred.log[mask]
    }
  }
}
param_names <- names(merfish_laminar_model$fitted.parameters)
tpoint_mask <- grepl("tpoint", param_names) & grepl("baseline", param_names)
degs <- rep(0, n_genes)
tpoints <- list()
names(degs) <- genes
for (g in genes) {
  gmask <- grepl(g, param_names) & tpoint_mask
  print(param_names[gmask])
  cat("\n\n")
  degs[g] <- sum(gmask)
  if (any(gmask)) {
    tpoints[[g]] <- as.integer(c(merfish_laminar_model$fitted.parameters[gmask]))
  }
}
max_blocks <- max(degs) + 1
data <- array(
  NA, 
  dim = c(
    max_blocks, n_genes,
    n_mice,
    length(trt_lvls)
  ),
  dimnames = list(bin = 1:max_blocks, child = genes, ran = mice, fixedeffect = trts)
)
for (m in mice) {
  for (t in trt_lvls[[m]]) {
    if (length(data0[, , m, t]) > 0) {
      for (g in genes) {
        if (degs[g] == 0) {
          data[1, g, m, t] <- mean(data0[, g, m, t])
        } else {
          for (i in 1:(degs[g] + 1)) {
            if (i == 1) {
              row_batch <- 1:(tpoints[[g]][i] - 1) 
            } else if (i > degs[g]) {
              row_batch <- tpoints[[g]][i - 1]:n_bins
            } else {
              row_batch <- tpoints[[g]][i - 1]:(tpoints[[g]][i] - 1) 
            }
            data[i, g, m, t] <- mean(data0[row_batch, g, m, t])
          }
        }
      }
    }
  }
}


# Find which treatment levels appear in all mice (are "universal") 
univ_trt <- rep(TRUE, length(trts))
names(univ_trt) <- trts
for (t in trts) {
  for (m in mice) {
    t_here <- FALSE
    mtrts <- trt_lvls[[m]]
    for (mt in mtrts) {
      mt_comps <- trts_comps[[mt]]
      t_here <- t_here || t %in% mt_comps
    }
    univ_trt[t] <- univ_trt[t] && t_here
  }
}

# For each universal treatment, find all treatment levels containing it
univ_lvls <- array(FALSE, dim = c(n_trts, n_trts), dimnames = list(trts, trts))
for (t in trts) {
  if (univ_trt[t]) {
    for (t2 in trts) {
      if (t %in% trts_comps[[t2]]) {
        univ_lvls[t, t2] <- TRUE
      }
    }
  }
}
univ_lvls <- univ_lvls[-c(which(rowSums(univ_lvls) == 0)),]
if (is.null(dim(univ_lvls)) && length(univ_lvls) > 0) {
  dim(univ_lvls) = c(1, length(univ_lvls))
  colnames(univ_lvls) <- trts
  rownames(univ_lvls) <- trts[univ_trt]
}

# For each fixed effect, check if it is universal, and if so, which treatment levels are universal within it
n_eff <- length(fix$name)
univ_eff <- rep(FALSE, n_eff)
names(univ_eff) <- fix$name
for (fe in fix$name) {
  for (t in trts) {
    if (univ_trt[t] && t %in% fix$lvls[[fe]]) {
      univ_eff[fe] <- TRUE
      rownames(univ_lvls)[rownames(univ_lvls) == t] <- fe
    }
  }
}

# For each fixed effect
ran_variation <- list()
for (fe in fix$name) {
  running_tally <- array(
    as.character(NA),
    dim = c(2, 4)
  )
  norm_diffs <- c()
  # ... then for each treatment level
  for (t in trts) {
    t1 <- t
    t2 <- t # assuming universal effect
    if (!univ_eff[fe]) {
      # ... but if not, pick treatment which matches t, except on fe
      this_fe <- trt_comp_array[t, fe]
      fe_col <- trt_comp_array[, fe]
      fe_mismatches <- fe_col != this_fe 
      these_fe <- trt_comp_array[t, fix$name != fe]
      these_col <- trt_comp_array[, fix$name != fe]
      if (!is.null(dim(these_col))) these_matches <- apply(these_col, 1, function(x) all(x == these_fe))
      else these_matches <- these_col == these_fe
      t2 <- trts[these_matches & fe_mismatches]
    }
    
    # For each unique mouse pair
    for (m1 in mice) {
      for (m2 in mice) {
        if (m1 != m2) {
          this_row <- c(m1, m2, t1, t2)
          if (any(apply(running_tally, 1, function(x) all(this_row %in% x)))) {
            next
          } else {
            values1 <- data[, , m1, t1]
            values2 <- data[, , m2, t2]
            
            diffs <- abs(values1 - values2)
            if (sum(!is.na(c(diffs))) > 0) {
              
              diffs <- diffs/values1
              diffs <- na.omit(c(diffs))
              norm_diffs <- c(norm_diffs, diffs)
              running_tally <- rbind(running_tally, this_row)
            }
          }
          
        }
      }
    }
  }
  ran_variation[[fe]] <- norm_diffs
}
# ... if that level is from the universal effect

hemiT <- ran_variation[[1]]
ageT <- ran_variation[[2]]

df <- data.frame(
  value = c(hemiT, ageT),
  group = c(rep("Ran.only", length(hemiT)), rep("Age+Ran", length(ageT)))
)

# Plot
ggplot(df, aes(x = value, color = group)) +
  geom_density(size = 1) +
  labs(x = "Value", y = "Density", title = "Density Plot of Count Varation") +
  theme_minimal()

max(ageT)
max(hemiT)
mean(ageT)
mean(hemiT)
mean(ageT)/mean(hemiT)


length(ageT)
length(hemiT)
