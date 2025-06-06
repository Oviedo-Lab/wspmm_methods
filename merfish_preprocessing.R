# Preprocessing MERFISH data ###########################################################################################

# Packages for plotting
library(ggplot2)
library(grid)
library(gridExtra)

# hdf5r package for extracting data
library(hdf5r)

# Report 
snk.report...("Loading helper function definitions for preprocessing MERFISH data", initial_breaks = 2)
RemoveL1 <- TRUE

# Load raw data ########################################################################################################

# Helper function, make rotation matrix
rot_matrix <- function(
    theta               # angle to rotate coordinates, in radians
  ) {
    return(
      matrix(
        c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
        nrow = 2, 
        byrow = TRUE
        )
      )
  }

# Helper function, matrix transform
matrix_transform <- function(
    coord,              # Set of coordinates to rotate, columns as dimensions, rows as points
    transform           # A matrix, or angle in radians
  ) {
    if (is.null(dim(transform))) transform <- rot_matrix(transform)
    return(as.matrix(coord) %*% t(transform))
  }

# Helper function for parameterized affine transforms
# ... not currently used, leaving here for future utility
affine <- function(
    coord,              # Set of coordinates to rotate, columns as dimensions, rows as points
    scale_H = 1,
    scale_W = 1,
    theta = 0,          # angle to rotate coordinates, in radians
    shear_x = 0,
    shear_y = 0
  ) {
    scale_matrix <- matrix(c(scale_W, 0, 0, scale_H), nrow = 2)
    rot_matrix <- rot_matrix(theta)
    shear_matrix <- matrix(c(1, shear_y, shear_x, 1), nrow = 2)
    return(as.matrix(coord) %*% t(scale_matrix %*% rot_matrix %*% shear_matrix))
  }

# Helper function to load and parse data, HDF5
parse_hdf5 <- function(
    file_path,
    mouse_num,
    z_view_bottom = TRUE,
    raw = TRUE,             # Grab normalized transcript counts, or raw ones?
    ROI_only = TRUE,        # Grab all transcripts, or just those in ROI?
    remove_L1 = RemoveL1,   # Exclude transcripts in layer 1?
    ROIname = "Primary auditory area"
  ) {
    
    # Load data
    file <- H5File$new(file_path, mode = "r")
    if (!raw) warning("Using normalized counts for hdf5 parsing")
    
    # Cell type codes for each cell
    celltype_MMC <- file[["/obs/MapMyCells_class/codes"]][] + 1          # Starts at 0, want numbers to match category names
    cellsubclass_MMC <- file[["/obs/MapMyCells_subclass/codes"]][] + 1   # Starts at 0, want numbers to match category names
    cellsupertype_MMC <- file[["/obs/MapMyCells_supertype/codes"]][] + 1 # Starts at 0, want numbers to match category names
    
    # Cell type names and their codes
    celltype_names <- gsub("^...", "", file[["/obs/MapMyCells_class/categories"]][])
    cellsubclass_names <- gsub("^...", "", file[["/obs/MapMyCells_subclass/categories"]][])
    cellsupertype_names <- gsub("^...", "", file[["/obs/MapMyCells_supertype/categories"]][])
    celltype_names[celltype_names == " Transcripts (filtered)"] <- "TSCRPT_filtered"
    celltype_names[celltype_names == " Bootstrapping Probability (filtered)"] <- "filtered"
    cellsubclass_names[cellsubclass_names == " Bootstrapping Probability (filtered)"] <- "filtered"
    cellsupertype_names[cellsupertype_names == " Bootstrapping Probability (filtered)"] <- "filtered"
    
    # Gene names
    nonblanks <- file[["/var/Genes"]][]
    gene_names <- file[["/var/_index"]][nonblanks]
    
    # Raw transcript counts, rows are cells, columns genes
    if (raw) transcript_counts_raw <- t(file[["/raw/X"]][,])     # should have integer elements, not normalized
    else transcript_counts_raw <- t(file[["/X"]][,])             # should have integer elements, normalized for MapMyCells
    transcript_counts_raw <- transcript_counts_raw[, nonblanks]  # drop blanks
    colnames(transcript_counts_raw) <- gene_names               
    n_cells <- nrow(transcript_counts_raw)
    
    # Grab rows for cortical layers
    ROI_names <- list(
      LL1 = paste0("/obs/ROI__Left ", ROIname, ", layer 1"),
      LL23 = paste0("/obs/ROI__Left ", ROIname, ", layer 2/3"),
      LL4 = paste0("/obs/ROI__Left ", ROIname, ", layer 4"),
      LL5 = paste0("/obs/ROI__Left ", ROIname, ", layer 5"),
      LL6a = paste0("/obs/ROI__Left ", ROIname, ", layer 6a"),
      LL6b = paste0("/obs/ROI__Left ", ROIname, ", layer 6b"),
      RL1 = paste0("/obs/ROI__Right ", ROIname, ", layer 1"),
      RL23 = paste0("/obs/ROI__Right ", ROIname, ", layer 2/3"),
      RL4 = paste0("/obs/ROI__Right ", ROIname, ", layer 4"),
      RL5 = paste0("/obs/ROI__Right ", ROIname, ", layer 5"),
      RL6a = paste0("/obs/ROI__Right ", ROIname, ", layer 6a"),
      RL6b = paste0("/obs/ROI__Right ", ROIname, ", layer 6b")
    )
    L1_mask <- file[[ROI_names[["LL1"]]]][] | file[[ROI_names[["RL1"]]]][]
    L23_mask <- file[[ROI_names[["LL23"]]]][] | file[[ROI_names[["RL23"]]]][]
    L4_mask <- file[[ROI_names[["LL4"]]]][] | file[[ROI_names[["RL4"]]]][]
    L5_mask <- file[[ROI_names[["LL5"]]]][] | file[[ROI_names[["RL5"]]]][]
    L6a_mask <- file[[ROI_names[["LL6a"]]]][] | file[[ROI_names[["RL6a"]]]][]
    L6b_mask <- file[[ROI_names[["LL6b"]]]][] | file[[ROI_names[["RL6b"]]]][]
    layer <- rep("notROI", n_cells)
    layer[L1_mask] <- "L1"
    layer[L23_mask] <- "L23"
    layer[L4_mask] <- "L4"
    layer[L5_mask] <- "L5"
    layer[L6a_mask] <- "L6a"
    layer[L6b_mask] <- "L6b"
    layer <- as.factor(layer)
    
    # Grab rows for left and right hemisphere
    left_mask <- file[[ROI_names[["LL1"]]]][] | 
      file[[ROI_names[["LL23"]]]][] | 
      file[[ROI_names[["LL4"]]]][] | 
      file[[ROI_names[["LL5"]]]][] | 
      file[[ROI_names[["LL6a"]]]][] | 
      file[[ROI_names[["LL6b"]]]][]
    right_mask <- file[[ROI_names[["RL1"]]]][] | 
      file[[ROI_names[["RL23"]]]][] | 
      file[[ROI_names[["RL4"]]]][] | 
      file[[ROI_names[["RL5"]]]][] | 
      file[[ROI_names[["RL6a"]]]][] | 
      file[[ROI_names[["RL6b"]]]][]
    hemisphere <- rep("notROI", n_cells)
    hemisphere[left_mask] <- "left"
    hemisphere[right_mask] <- "right"
    hemisphere <- as.factor(hemisphere)
    
    # Grab metadata and form columns
    age <- as.factor(rep(as.integer(sub("P", "", file[["/uns/info/metadata/age"]][])), n_cells))
    sex <- as.factor(rep(file[["/uns/info/metadata/sex"]][], n_cells))
    strain <- as.factor(rep(file[["/uns/info/metadata/strain"]][], n_cells))
    experience <- as.factor(rep(file[["/uns/info/metadata/treatment"]][], n_cells))
    
    # Grab spatial coordinates and form columns 
    # ... units are in microns (um)
    # ... x_max <- file[["/obs/max_x"]][] 
    # ... x_min <- file[["/obs/min_x"]][]
    # ... y_max <- file[["/obs/max_y"]][] 
    # ... y_min <- file[["/obs/min_y"]][]
    x_coord <- file[["/obs/center_x"]][]
    y_coord <- file[["/obs/center_y"]][]
    
    # Assign a number to the mouse
    mouse <- as.factor(rep(mouse_num, n_cells))
    
    # Reorient spatial coordinates of whole slice 
    # ... center around middle of L5 left cortical area
    L5_left_mask <- L5_mask & left_mask
    x_coord <- x_coord - mean(x_coord[L5_left_mask])
    y_coord <- y_coord - mean(y_coord[L5_left_mask])
    # ... align x axis and line through L5 cortical regions
    L5_right_mask <- L5_mask & right_mask
    right_xlarger <- mean(x_coord[L5_right_mask]) > mean(x_coord[L5_left_mask])
    tilt_slope <- mean(y_coord[L5_right_mask]) / mean(x_coord[L5_right_mask])
    if (tilt_slope < 0) {
      y_tilt_radians <- -atan(-tilt_slope)
    } else {
      y_tilt_radians <- atan(tilt_slope)
    }
    aligned_coord <- matrix_transform(
      coord = cbind(x_coord, y_coord), 
      transform = y_tilt_radians
    )
    # ... save transformed coordinates 
    x_coord <- aligned_coord[, 1]
    y_coord <- aligned_coord[, 2]
    # ... align y axis and perpendicular bisection of cortical regions 
    x_coord <- x_coord - mean(x_coord[L5_mask])
    # ... make sure anterior is up
    if (z_view_bottom) {
      if (right_xlarger) anterior_up <- FALSE
      else anterior_up <- TRUE
    } else {
      if (right_xlarger) anterior_up <- TRUE
      else anterior_up <- FALSE
    }
    if (!anterior_up) y_coord <- -y_coord
    # ... make sure right is positive
    if (!right_xlarger) x_coord <- -x_coord
    # ... push into positive quadrant corner
    x_coord <- x_coord - min(x_coord)
    y_coord <- y_coord - min(y_coord)
    
    # Make data frame
    transcript_counts <- data.frame(
      mouse, 
      celltype_MMC, # cellsubclass_MMC, cellsupertype_MMC, 
      hemisphere, layer, 
      age, sex, 
      strain, experience, 
      x_coord, y_coord,
      transcript_counts_raw 
    )
    
    # Make plot of whole cortical slice 
    layer_colors <- c("notROI" = "gray", "L1" = "gray4", "L23" = "tomato1", "L4" = "orange", "L5" = "springgreen2", "L6a" = "steelblue1", "L6b" = "purple") 
    slice_plot <- ggplot(transcript_counts, aes(x = x_coord, y = y_coord, color = factor(layer))) +
      geom_point() +
      scale_color_manual(values = layer_colors) +
      labs(title = paste("Mouse", mouse_num), color = "Layer") +
      theme_minimal() + theme(legend.position = "none")
    
    # Prune to just ROI
    if (ROI_only) {
      transcript_counts <- transcript_counts[hemisphere != "notROI",]
    }
    
    # Replace cell type number with cell type name
    # ... transcript_counts$cellsubclass_MMC <- cellsubclass_names[transcript_counts$cellsubclass_MMC]
    # ... transcript_counts$cellsupertype_MMC <- cellsupertype_names[transcript_counts$cellsupertype_MMC]
    transcript_counts$celltype_MMC <- celltype_names[transcript_counts$celltype_MMC]
    
    # Collapse GABA cell types together 
    GABA_mask <- which(transcript_counts$celltype_MMC == "CTX-CGE GABA" | transcript_counts$celltype_MMC == "CTX-MGE GABA")
    transcript_counts$celltype_MMC[GABA_mask] <- "CTX-CGE/MGE GABA"
    
    # Convert to factor 
    transcript_counts$celltype_MMC <- as.factor(transcript_counts$celltype_MMC)
    
    # Remove Layer 1
    if (remove_L1) {
      L1_mask <- transcript_counts$layer == "L1"
      transcript_counts <- transcript_counts[!L1_mask,]
    }
    
    # Release file
    file$close_all()
    
    return(list(transcript_counts = transcript_counts, slice_plot = slice_plot))
    
  }

# Helper function to load and parse data, csv
parse_csv <- function(
    file_path,
    mouse_num,
    z_view_bottom = TRUE,
    remove_L1 = RemoveL1    # Exclude transcripts in layer 1?  
  ) {
    
    # Load data
    file <- data.table::fread(file_path, showProgress = FALSE)
    
    # Make transcript count column
    # ... rows in "file" are transcripts
    n_transcripts <- nrow(file)
    transcript_counts_raw <- rep(1, n_transcripts)
    
    # Make layer column 
    layer <- file$roi_labels
    layer[layer == "Primary auditory area, layer 1"] <- "L1"
    layer[layer == "Primary auditory area, layer 2/3"] <- "L23"
    layer[layer == "Primary auditory area, layer 4"] <- "L4"
    layer[layer == "Primary auditory area, layer 5"] <- "L5"
    layer[layer == "Primary auditory area, layer 6a"] <- "L6a"
    layer[layer == "Primary auditory area, layer 6b"] <- "L6b"
    layer <- as.factor(layer)
    
    # Make hemisphere column
    hemisphere <- as.factor(file$hemi_labels)
    
    # Make list of ages 
    age_list <- list(
      "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/data_ACx_transcripts//transcripts-CCFreg_ACxDev1_CBA_CaJ_FFPE_P12_z0-3.csv" = 12,
      "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/data_ACx_transcripts//transcripts-CCFreg_ACxDev1_CBA_CaJ_FFPE_P18.csv" = 18,
      "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/data_ACx_transcripts//transcripts-CCFreg_ACxDev1_CBA_CaJ_FFPE_P7.csv" = 7,
      "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/data_ACx_transcripts//transcripts-CCFreg_ACxDev1_CBA_CaJ_P12.csv" = 12, 
      "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/data_ACx_transcripts//transcripts-CCFreg_ACxDev1_CBA_CaJ_P18.csv" = 18
    )
    
    # Grab metadata and form columns
    age <- as.factor(rep(age_list[[file_path]], n_transcripts))
    sex <- as.factor(rep("male", n_transcripts))
    strain <- as.factor(rep("CBA-CAJ", n_transcripts))
    experience <- as.factor(rep("control", n_transcripts))
    
    # Grab spatial coordinates and form columns 
    # ... units are in microns (um)
    x_coord <- file$global_x
    y_coord <- file$global_y
    
    # Assign a number to the mouse
    mouse <- as.factor(rep(mouse_num, n_transcripts))
    
    # Reorient spatial coordinates of whole slice 
    L5_mask <- layer == "L5"
    left_mask <- hemisphere == "left"
    right_mask <- hemisphere == "right"
    # ... center around middle of L5 left cortical area
    L5_left_mask <- L5_mask & left_mask
    x_coord <- x_coord - mean(x_coord[L5_left_mask])
    y_coord <- y_coord - mean(y_coord[L5_left_mask])
    # ... align x axis and line through L5 cortical regions
    L5_right_mask <- L5_mask & right_mask
    right_xlarger <- mean(x_coord[L5_right_mask]) > mean(x_coord[L5_left_mask])
    tilt_slope <- mean(y_coord[L5_right_mask]) / mean(x_coord[L5_right_mask])
    if (tilt_slope < 0) {
      y_tilt_radians <- -atan(-tilt_slope)
    } else {
      y_tilt_radians <- atan(tilt_slope)
    }
    aligned_coord <- matrix_transform(
      coord = cbind(x_coord, y_coord), 
      transform = y_tilt_radians
    )
    # ... save transformed coordinates 
    x_coord <- aligned_coord[, 1]
    y_coord <- aligned_coord[, 2]
    # ... align y axis and perpendicular bisection of cortical regions 
    x_coord <- x_coord - mean(x_coord[L5_mask])
    # ... make sure anterior is up
    if (z_view_bottom) {
      if (right_xlarger) anterior_up <- FALSE
      else anterior_up <- TRUE
    } else {
      if (right_xlarger) anterior_up <- TRUE
      else anterior_up <- FALSE
    }
    if (!anterior_up) y_coord <- -y_coord
    # ... make sure right is positive
    if (!right_xlarger) x_coord <- -x_coord
    # ... push into positive quadrant corner
    x_coord <- x_coord - min(x_coord)
    y_coord <- y_coord - min(y_coord)
    
    # Make celltype columns
    cell_id <- file$cell_id
    celltype_MMC <- as.factor(rep("unknown", n_transcripts))
    celltype <- rep("intracellular", n_transcripts)
    celltype[file$cell_id == -1] <- "extracellular"
    celltype <- as.factor(celltype)
    celltype <- relevel(celltype, ref = "intracellular")
    
    # Make gene column 
    trscrpt_gene_symb <- file$gene
    
    # Make data frame
    transcript_counts <- data.frame(
      mouse, 
      celltype_MMC, celltype, cell_id, 
      hemisphere, layer, 
      age, sex, 
      strain, experience, 
      x_coord, y_coord,
      transcript_counts_raw, 
      trscrpt_gene_symb
    )
    
    # Make plot of whole cortical slice 
    transcript_counts_downsampled <- transcript_counts[sample(1:nrow(transcript_counts), 10000),]
    layer_colors <- c("notROI" = "gray", "L1" = "gray4", "L23" = "tomato1", "L4" = "orange", "L5" = "springgreen2", "L6a" = "steelblue1", "L6b" = "purple") 
    slice_plot <- ggplot(transcript_counts_downsampled, aes(x = x_coord, y = y_coord, color = factor(layer))) +
      geom_point() +
      scale_color_manual(values = layer_colors) +
      labs(title = paste("Mouse", mouse_num), color = "Layer") +
      theme_minimal() + theme(legend.position = "none")
    
    # Remove Layer 1
    if (remove_L1) {
      L1_mask <- transcript_counts$layer == "L1"
      transcript_counts <- transcript_counts[!L1_mask,]
    }
    
    return(list(transcript_counts = transcript_counts, slice_plot = slice_plot))
    
  }

# Load and parse data, HDF5
make_count_data <- function(
    data_path, 
    remove_L1 = RemoveL1,   # Exclude transcripts in layer 1?
    ROIname = "Primary auditory area",
    raw = TRUE,             # Grab normalized transcript counts, or raw ones?
    verbose = TRUE
  ) {
    
    # Get a list of all HDF5 files in the "data" folder
    files <- list.files(
      path = data_path, # Defined in the main script
      pattern = "\\_withCCF.hdf5$", 
      full.names = TRUE
    )
    names(files) <- paste("mouse", seq_along(files)) # These numbers will correspond to the ran levels assigned latter
    if (verbose) {
      snk.report("Loading raw data")
      snk.horizontal_rule(reps = snk.simple_break_reps)
      snk.report...(paste("Found", length(files), "HDF5 files."))
      snk.print_var_list("File names",files)
    }
    
    # Loop through each file and parse it
    cat("Loading and parsing file for mouse number: ")
    mean_rates <- c()
    cells_per_mouse <- c()
    slice_plots <- list()
    for (f in seq_along(files)) {
      if (f < length(files)) cat(f, ", ", sep="")
      else cat(f, "\n")
      assign(paste0("mouse", f), parse_hdf5(file_path = files[f], mouse_num = f, raw = raw, remove_L1 = remove_L1, ROIname = ROIname))
      slice_plots[[paste0("slice_plot", f)]] <- get(paste0("mouse", f))$slice_plot
      assign(paste0("mouse", f), get(paste0("mouse", f))$transcript_counts)
      ind_var_fields <- which(colnames(get(paste0("mouse", f))) %in% c(
        "mouse", 
        "celltype_MMC", # "cellsubclass_MMC", "cellsupertype_MMC", 
        "hemisphere", "layer", 
        "age", "sex", 
        "strain", "experience",
        "x_coord", "y_coord"))
      noise_columns <- grep("_noise$", colnames(get(paste0("mouse", f))), value = FALSE)
      mean_rates <- c(mean_rates, mean(as.matrix(get(paste0("mouse", f))[,-c(ind_var_fields, noise_columns)])))
      cells_per_mouse <- c(cells_per_mouse, nrow(get(paste0("mouse", f))))
    }
    
    # Check the mean of each run to see if there's some reason to suspect a systematic difference
    if (verbose) {
      snk.print_vec("Mean transcript counts per cell for each mouse", mean_rates)
      snk.print_vec("Mean of means", mean(mean_rates))
      snk.print_vec("Standard deviation of means", sd(mean_rates))
    }
    
    # Combine data from all mice
    count_data <- data.frame()
    for (f in 1:length(files)) {
      count_data <- rbind(count_data, get(paste0("mouse", f)))
      rm(list = paste0("mouse", f))
    }
    count_data$hemisphere <- droplevels(count_data$hemisphere)
    count_data$hemisphere <- relevel(count_data$hemisphere, ref = "left")
    
    # Print summary
    if (verbose) {
      snk.print_table(
        "Total cells by type (Map My Cells)", 
        table(count_data$celltype_MMC),
        head = FALSE,
        initial_breaks = 1
      )
      snk.print_vec("Number of mice", length(files))
      snk.print_vec("Cells per mouse", cells_per_mouse)
      snk.print_vec("Means cells per mouse", mean(cells_per_mouse))
      snk.print_vec("Total cells", nrow(count_data))
    }
    
    # Initialize new coordinate columns 
    count_data$x_trans <- rep(0,nrow(count_data))
    count_data$y_trans <- rep(0,nrow(count_data))
    count_data$x_bins_raw <- rep(0,nrow(count_data))
    count_data$y_bins_raw <- rep(0,nrow(count_data))
    count_data$x_bins <- rep(0,nrow(count_data))
    count_data$y_bins <- rep(0,nrow(count_data))
    
    return(list(count_data = count_data, slice_plots = slice_plots))
    
  }

# Load and parse data, csv
make_count_data_csv <- function(
    data_path, 
    remove_L1 = RemoveL1,   # Exclude transcripts in layer 1?
    verbose = TRUE
  ) {
    
    # Get a list of all csv files in the "data" folder
    files <- list.files(
      path = data_path, # Defined in the main script
      pattern = "\\.csv$", 
      full.names = TRUE
    )
    names(files) <- paste("mouse", seq_along(files)) # These numbers will correspond to the ran levels assigned latter
    if (verbose) {
      snk.report("Loading raw data")
      snk.horizontal_rule(reps = snk.simple_break_reps)
      snk.report...(paste("Found", length(files), "files."))
      snk.print_var_list("File names",files)
    }
    
    # Loop through each file and parse it
    slice_plots <- list()
    count_list <- vector("list", length(files))
    cat("Loading and parsing file for mouse number: ")
    for (f in seq_along(files)) {
      if (f < length(files)) cat(f, ", ", sep="") else cat(f, "\n")
      parsed <- parse_csv(file_path = files[f], mouse_num = f, remove_L1 = remove_L1)
      slice_plots[[paste0("slice_plot", f)]] <- parsed$slice_plot
      count_list[[f]] <- parsed$transcript_counts
    }
    
    # Combine data from all mice
    count_data <- do.call(rbind, count_list)
    count_data$hemisphere <- droplevels(count_data$hemisphere)
    count_data$hemisphere <- relevel(count_data$hemisphere, ref = "left")
    
    # Print summary
    if (verbose) {
      snk.print_vec("Number of mice", length(files))
      snk.print_vec("Total transcripts", nrow(count_data))
    }
    
    # Initialize new coordinate columns 
    count_data$x_trans <- rep(0,nrow(count_data))
    count_data$y_trans <- rep(0,nrow(count_data))
    count_data$x_bins_raw <- rep(0,nrow(count_data))
    count_data$y_bins_raw <- rep(0,nrow(count_data))
    count_data$x_bins <- rep(0,nrow(count_data))
    count_data$y_bins <- rep(0,nrow(count_data))
    
    return(list(count_data = count_data, slice_plots = slice_plots))
    
  }

# Functions to transform into laminar and columnar coordinates #########################################################

# Note: 
#  laminar axis goes to y, with cortical surface in the positive direction (zero is deepest layer, L6)
#  columnar axis goes to x, with anterior in the positive direction (zero is most posterior)

# Helper function for plotting
plot_results <- function(
    df, 
    coordinates,         # Set of coordinates to rotate, columns as dimensions, rows as points
    mouse_num,
    plot_title,
    separate_hemi = TRUE
  ) {
   
    # Check if mask has been applied 
    mask <- df$mouse == mouse_num 
    df <- df[mask,]
    if (nrow(df) != nrow(coordinates)) coordinates <- coordinates[mask,]
    
    # Force correct coordinate labels 
    x_coord <- "x_coord"
    y_coord <- "y_coord"
    colnames(coordinates) <- c(x_coord, y_coord)
    range_limit <- range(c(coordinates[, x_coord], coordinates[, y_coord]), na.rm = TRUE)
    
    # Grab masks
    mask_left <- df[,"hemisphere"] == "left"
    mask_right <- df[,"hemisphere"] == "right"
    
    # Get right hemisphere data
    layer_right <- df$layer[mask_right]
    df_right <- coordinates[mask_right,]
    if (separate_hemi) {
      df_right[,x_coord] <- df_right[,x_coord] - min(df_right[,x_coord])
      df_right[,y_coord] <- df_right[,y_coord] - min(df_right[,y_coord])
    }
    df_right <- cbind(df_right, layer_right)
    colnames(df_right)[ncol(df_right)] <- "layer"
    df_right$hemisphere <- "right"
    
    # Get left hemisphere data
    layer_left <- df$layer[mask_left]
    df_left <- coordinates[mask_left,]
    if (separate_hemi) {
      df_left[,x_coord] <- df_left[,x_coord] - min(df_left[,x_coord])
      df_left[,y_coord] <- df_left[,y_coord] - min(df_left[,y_coord])
    }
    df_left <- cbind(df_left, layer_left)
    colnames(df_left)[ncol(df_left)] <- "layer"
    df_left$hemisphere <- "left"
    
    # Combine data
    layers_both <- rbind(df_right, df_left)
    
    # Downsample 
    if (nrow(layers_both) > 10000) {
      layers_both <- layers_both[sample(1:nrow(layers_both), 10000),]
    }
    
    plot <- ggplot(layers_both, aes(x = x_coord, y = y_coord, color = layer)) +
      geom_point(na.rm = TRUE) +
      facet_wrap(~hemisphere) +
      labs(x = "X Coordinate", y = "Y Coordinate", color = "Layer", title = plot_title) +
      theme_minimal()
    
    if (!separate_hemi) {
      plot <- plot + 
        scale_x_continuous(limits = range_limit) + 
        scale_y_continuous(limits = range_limit)
    }
    
    return(plot)
    
  }

# Define helper function for transforming coordinates
coordinate_transform <- function(
    mouse_num, 
    df,
    nat_left = FALSE,    # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    nat_right = FALSE,   # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    verbose = TRUE
  ) {
    
    if (verbose) snk.report...("Grabbing coordinates and defining layers and hemispheres")
    
    # Grab indexes
    mask <- df$mouse == mouse_num 
    mask_left <- df[mask,"hemisphere"] == "left"
    mask_right <- df[mask,"hemisphere"] == "right"
    
    # Define columns for coordinate
    x_coord <- "x_coord"
    y_coord <- "y_coord"
    all_coord <- c(x_coord, y_coord)
    
    # Grab coordinates 
    coordinates <- df[mask, all_coord]
    
    # Define layer rows
    mask_left_L1 <- mask_left & df[mask, "layer"] == "L1"
    mask_right_L1 <- mask_right & df[mask, "layer"] == "L1"
    mask_left_L23 <- mask_left & df[mask, "layer"] == "L23"
    mask_right_L23 <- mask_right & df[mask, "layer"] == "L23"
    mask_left_L4 <- mask_left & df[mask, "layer"] == "L4"
    mask_right_L4 <- mask_right & df[mask, "layer"] == "L4"
    mask_left_L5 <- mask_left & df[mask, "layer"] == "L5"
    mask_right_L5 <- mask_right & df[mask, "layer"] == "L5"
    mask_left_L6a <- mask_left & df[mask, "layer"] == "L6a"
    mask_right_L6a <- mask_right & df[mask, "layer"] == "L6a"
    mask_left_L6b <- mask_left & df[mask, "layer"] == "L6b"
    mask_right_L6b <- mask_right & df[mask, "layer"] == "L6b"
    
    # Put layer rows in convenient lists
    # ... left hemisphere
    layer_rows_left <- list(
      L1 = mask_left_L1,
      L23 = mask_left_L23,
      L4 = mask_left_L4,
      L5 = mask_left_L5,
      L6a = mask_left_L6a,
      L6b = mask_left_L6b
    )
    # ... right hemisphere
    layer_rows_right <- list(
      L1 = mask_right_L1,
      L23 = mask_right_L23,
      L4 = mask_right_L4,
      L5 = mask_right_L5,
      L6a = mask_right_L6a,
      L6b = mask_right_L6b
    )
    
    # Plot original data
    range_limit <- range(c(df[, x_coord], df[, y_coord]), na.rm = TRUE)
    if (nrow(df) > 10000) {
      these_rows <- sample(1:nrow(df), 10000)
    } else {
      these_rows <- 1:nrow(df)
    }
    df2 <- df[these_rows, ]
    plot_untransformed <- ggplot(df2[df2$mouse == mouse_num,], aes(x = x_coord, y = y_coord, color = layer)) +
      geom_point(na.rm = TRUE) +
      labs(
        x = "X Coordinate", 
        y = "Y Coordinate", 
        color = "Layer", 
        title = paste("Untransformed Cortical layers, mouse", mouse_num)
      ) +
      theme_minimal() +
      scale_x_continuous(limits = range_limit) + 
      scale_y_continuous(limits = range_limit)
    plot_list <- list(plot_untransformed = plot_untransformed)
    
    # Step 1: Center each patch (left and right) around the mean point of L5
    if (verbose) snk.report...("Step 1, centering each patch around the mean point of L5")
    # Helper function for recentering data
      recenter_coordinates <- function(
      coord, 
      mask_hemisphere = NULL, 
      mask_layer = NULL
    ) {
      
      if (is.null(mask_layer)) mask_layer <- TRUE
      if (is.null(mask_hemisphere)) mask_hemisphere <- TRUE
      
      mean_x <- mean(coord[mask_layer, x_coord])
      mean_y <- mean(coord[mask_layer, y_coord])
      
      coord[mask_hemisphere,x_coord] <- coord[mask_hemisphere,x_coord] - mean_x
      coord[mask_hemisphere,y_coord] <- coord[mask_hemisphere,y_coord] - mean_y
      
      return(coord)
      
    }
    # ... perform recentering
    coordinates <- recenter_coordinates(coordinates, mask_right, mask_right_L5)
    coordinates <- recenter_coordinates(coordinates, mask_left, mask_left_L5)
    
    # Test by plotting (recentered)
    plot_recenter <- plot_results(
      df, coordinates, mouse_num,
      paste("Untransformed cortical layers (recentered), mouse", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_recenter = plot_recenter))
    
    # Helper function for leveling
    level_layer <- function(
      data, 
      hemisphere, 
      mask_hemisphere, 
      mask_layer, 
      flip_right = TRUE,
      natural_left = FALSE,
      natural_right = FALSE,
      verbose = FALSE
    ) {
      # Fit linear model to the mask_layer coordinates to get the angle of the slice
      model <- lm(y_coord ~ x_coord, data = data[mask_layer,])
      # Find slope and angle
      slope <- model$coefficients[x_coord]
      if (slope < 0) {
        angle <- -atan(slope)
        if (natural_left || natural_right) angle <- -angle
      } else {
        angle <- atan(slope)
      }
      if (verbose) cat("\nangle: ", angle*57.3)
      # Center layer on x axis 
      y_mean <- mean(data[mask_layer, y_coord])
      new_coord <- data[mask_hemisphere,]
      colnames(new_coord) <- colnames(data)
      new_coord[, y_coord] <- new_coord[, y_coord] - y_mean
      # Level by rotating
      new_coord <- matrix_transform(new_coord, angle)
      colnames(new_coord) <- colnames(data)
      # Undo centering
      new_coord[, y_coord] <- new_coord[, y_coord] + y_mean
      # Flip hemispheres if necessary
      if (hemisphere == "right" && flip_right && !natural_right) new_coord[, y_coord] <- -new_coord[, y_coord]
      if (hemisphere == "right" && flip_right && natural_right) new_coord[, x_coord] <- -new_coord[, x_coord]
      if (natural_left) new_coord[, y_coord] <- -new_coord[, y_coord]
      # Return just the transformed hemisphere
      return(new_coord)
    }
    
    # Step 2: Rotate each patch so that L4 aligns with the x-axis and L1 (or L23, if L1 removed) is on top
    if (verbose) snk.report...("Step 2, rotating each patch so that L4 aligns with the x-axis with anterior in positive y direction")
    coordinates[mask_right,] <- level_layer(coordinates, "right", mask_right, mask_right_L4, natural_right = nat_right, verbose = FALSE)
    coordinates[mask_left,] <- level_layer(coordinates, "left", mask_left, mask_left_L4, natural_left = nat_left, verbose = FALSE)
    
    # Test by plotting (L4 leveled)
    plot_level_L4 <- plot_results(
      df, coordinates, mouse_num,
      paste("Transformed cortical layers (L4 leveled), mouse", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_level_L4 = plot_level_L4))
    
    # Step 3: Model and flatten laminar curve based on L4
    if (verbose) snk.report...("Step 3, modeling and flattening laminar curve based on L4")
    
    # Helper function for flattening data
    flatten <- function(
      data, 
      mask_hemisphere, 
      mask_layer
    ) {
      
      # Down-sample data 
      just_these <- sample(1:sum(mask_layer), as.integer(sum(mask_layer)/2))
      data_sampled <- data[mask_layer,][just_these,]
      
      # Fit curve
      model <- nls(
        y_coord ~ -c * (x_coord - a)^2 + b,
        data = data_sampled,
        start = list(a = 0, b = 1, c = 0.001)
      )
      
      # Get coefficients
      pars <- model$m$getPars()
      
      # Define function to flatten data 
      uw <- function(x, b, a, c) {
        c*(x - a)^2 + b
      }
      
      # Flatten this hemisphere 
      data[mask_hemisphere, y_coord] <- uw(
        x = data[mask_hemisphere, x_coord], 
        b = data[mask_hemisphere, y_coord], 
        pars["a"], pars["c"]
        )
      
      return(data[mask_hemisphere,])
      
    }
    
    # Flatten each hemisphere
    coordinates[mask_right,] <- flatten(coordinates, mask_right, mask_right_L4)
    coordinates[mask_left,] <- flatten(coordinates, mask_left, mask_left_L4)
    
    # Test by plotting (flattened)
    plot_flattened <- plot_results(
      df, coordinates, mouse_num,
      paste("Transformed cortical layers (curve flattened), mouse", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_flattened = plot_flattened))
    
    # Translate so all points are positive 
    coordinates[,x_coord] <- coordinates[,x_coord] - min(coordinates[,x_coord])
    coordinates[,y_coord] <- coordinates[,y_coord] - min(coordinates[,y_coord])
    
    return(list(coord = coordinates, plot_list = plot_list))
    
  }

# Define helper function for binning coordinates and smoothing edges
coordinate_binning <- function(
    mouse_num,
    total_bins,             # Number of bins to use when binning data
    layer_names,            # e.g., c("L1", "L23", "L4", "L5", "L6a", "L6b")
    df,
    L1_removed = RemoveL1,  # Exclude transcripts in layer 1?
    nat_left = FALSE,       # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    nat_right = FALSE,      # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    verbose = TRUE
  ) {
    
    # Grab rows for the mouse
    mask <- df$mouse == mouse_num 
    mask_left <- mask & df[,"hemisphere"] == "left"
    mask_right <- mask & df[,"hemisphere"] == "right"
    
    # Apply coordinate transform to those rows
    if (verbose) snk.report...("Performing coordinate transform")
    coord_trans <- coordinate_transform(mouse_num, df, nat_left = nat_left, nat_right = nat_right)
    plot_list <- coord_trans$plot_list
    coord_trans <- coord_trans$coord
    
    # Find bin boundaries
    max_x <- max(coord_trans[,1])
    max_y <- max(coord_trans[,2])
    bin_ticks_x <- seq(0, max_x + 1, length.out = total_bins+1)
    bin_ticks_y <- seq(0, max_y + 1, length.out = total_bins+1)
    
    # Bin the transformed coordinates
    if (verbose) snk.report...("Binning transformed coordinates")
    x_bins <- findInterval(coord_trans[, 1], bin_ticks_x)
    y_bins <- findInterval(coord_trans[, 2], bin_ticks_y)
    
    # Add bin numbers to df
    df[mask,"x_bins_raw"] <- x_bins
    df[mask,"y_bins_raw"] <- y_bins
    df[mask,"x_bins"] <- x_bins
    df[mask,"y_bins"] <- y_bins
    
    # Apply nonlinear smoothing to bins 
    if (verbose) snk.report...("Smoothing bin edges with nonlinear transformation")
    stretch_to_fill <- function(
      df_, 
      hemisphere_mask,
      upper,
      layer_names
    ) {
      
      # Grab layer to lead stretch
      if (upper && L1_removed) layer_num <- 2
      else if (upper && !L1_removed) layer_num <- 1
      else layer_num <- length(layer_names)
      
      # Grab index
      mask_ <- hemisphere_mask & as.character(df_$layer) == layer_names[layer_num] 
      
      # Find max/min y point in layer
      # ... For smoothing laminar edges:
      Ly_by_x <- rep(NA,total_bins)
      Ly_by_x_abs <- rep(NA,total_bins)
      # ... For smoothing columnar edges:
      Lx_by_y <- rep(NA,total_bins)
      if (upper) {
        # Max point in column in upper layer
        Ly_by_x_vec <- tapply(X = df_$y_bins_raw[mask_], INDEX = df_$x_bins_raw[mask_], FUN = max, na.rm = TRUE)
        Ly_by_x[as.integer(names(Ly_by_x_vec))] <- Ly_by_x_vec
        # Max point in column
        Ly_by_x_abs_vec <- tapply(X = df_$y_bins_raw[hemisphere_mask], INDEX = df_$x_bins_raw[hemisphere_mask], FUN = max, na.rm = TRUE)
        Ly_by_x_abs[as.integer(names(Ly_by_x_abs_vec))] <- Ly_by_x_abs_vec
        # Max point in layer
        Lx_by_y_vec <- tapply(X = df_$x_bins_raw[hemisphere_mask], INDEX = df_$y_bins_raw[hemisphere_mask], FUN = max, na.rm = TRUE)
        Lx_by_y[as.integer(names(Lx_by_y_vec))] <- Lx_by_y_vec
      } else {
        # Min point in column in lower layer
        Ly_by_x_vec <- tapply(X = df_$y_bins_raw[mask_], INDEX = df_$x_bins_raw[mask_], FUN = min, na.rm = TRUE)
        Ly_by_x[as.integer(names(Ly_by_x_vec))] <- Ly_by_x_vec
        # Min point in column
        Ly_by_x_abs_vec <- tapply(X = df_$y_bins_raw[hemisphere_mask], INDEX = df_$x_bins_raw[hemisphere_mask], FUN = min, na.rm = TRUE)
        Ly_by_x_abs[as.integer(names(Ly_by_x_abs_vec))] <- Ly_by_x_abs_vec
        # Min point in layer
        Lx_by_y_vec <- tapply(X = df_$x_bins_raw[hemisphere_mask], INDEX = df_$y_bins_raw[hemisphere_mask], FUN = min, na.rm = TRUE)
        Lx_by_y[as.integer(names(Lx_by_y_vec))] <- Lx_by_y_vec
      }
      
      # Look for any empty columns 
      idx_empty <- which(is.na(Ly_by_x))
      idx_full <- which(!is.na(Ly_by_x))
      idx_empty_abs <- which(is.na(Ly_by_x_abs))
      idx_full_abs <- which(!is.na(Ly_by_x_abs))
      for (i in idx_empty) {
        if (i == 1) Ly_by_x[i] <- Ly_by_x[idx_full[which.min(idx_full > i)]]
        else Ly_by_x[i] <- Ly_by_x[idx_full[which.min(idx_full < i)]]
      }
      for (i in idx_empty_abs) {
        if (i == 1) Ly_by_x_abs[i] <- Ly_by_x_abs[idx_full_abs[which.min(idx_full_abs > i)]]
        else Ly_by_x_abs[i] <- Ly_by_x_abs[idx_full_abs[which.min(idx_full_abs < i)]]
      }
      
      # Look for any empty layers 
      idx_empty <- which(is.na(Lx_by_y))
      idx_full <- which(!is.na(Lx_by_y))
      for (i in idx_empty) {
        if (i == 1) Lx_by_y[i] <- Lx_by_y[idx_full[which.min(idx_full > i)]]
        else Lx_by_y[i] <- Lx_by_y[idx_full[which.min(idx_full < i)]]
      }
      
      # Check bounds, columns
      if (upper) Ly_by_x[Ly_by_x < Ly_by_x_abs] <- Ly_by_x_abs[Ly_by_x < Ly_by_x_abs]
      else Ly_by_x[Ly_by_x > Ly_by_x_abs] <- Ly_by_x_abs[Ly_by_x > Ly_by_x_abs]
      
      # Preserve edge jitter, columns
      if (upper) Ly_by_x <- Ly_by_x + c(abs(diff(Ly_by_x)),0)/2
      else Ly_by_x <- Ly_by_x - c(abs(diff(Ly_by_x)),0)/2
      
      # Check edge jitter, layers 
      if (upper) Lx_by_y <- Lx_by_y + c(abs(diff(Lx_by_y)),0)/2
      else Lx_by_y <- Lx_by_y - c(abs(diff(Lx_by_y)),0)/2
      
      # Perform smoothing
      if (upper) {
        mask_range <- hemisphere_mask & df_$y_bins_raw > total_bins/2
        mask_range_layer <- hemisphere_mask & df_$x_bins_raw > total_bins/2
        to_ <- c(total_bins/2,total_bins)
      } else {
        mask_range <- hemisphere_mask & df_$y_bins_raw < total_bins/2
        mask_range_layer <- hemisphere_mask & df_$x_bins_raw < total_bins/2
        to_ <- c(1,total_bins/2+1)
      }
      x_groups <- split(seq_len(nrow(df_)), df_$x_bins_raw)
      y_groups <- split(seq_len(nrow(df_)), df_$y_bins_raw)
      for (b in 1:total_bins) {
        # ... columns
        idx_x <- x_groups[[as.character(b)]]
        idx_x <- idx_x[mask_range[idx_x]]  # apply mask
        if (upper) from_ <- c(total_bins/2, Ly_by_x[b]) 
        else from_ <- c(Ly_by_x[b], total_bins/2) 
        weight <- abs(total_bins/2 - df_$y_bins_raw[idx_x]) / 
          abs(total_bins/2 - Ly_by_x[b])
        df_$y_bins[idx_x] <- as.integer(
          df_$y_bins_raw[idx_x] * (1 - weight) +
            weight * scales::rescale(
              df_$y_bins_raw[idx_x], 
              to = to_, 
              from = from_
              )
          )
        # ... layers
        idx_y <- y_groups[[as.character(b)]]
        idx_y <- idx_y[mask_range_layer[idx_y]]  # apply mask
        if (upper) from_ <- c(total_bins/2, Lx_by_y[b])
        else from_ <- c(Lx_by_y[b], total_bins/2)
        weight <- abs(total_bins/2 - df_$x_bins_raw[idx_y]) / 
          abs(total_bins/2 - Lx_by_y[b])
        df_$x_bins[idx_y] <- as.integer(
          df_$x_bins_raw[idx_y] * (1 - weight) +
            weight * scales::rescale(
              df_$x_bins_raw[idx_y], 
              to = to_, 
              from = from_
              )
          )
      }
      
      return(df_)
      
    }
    
    # stretch upper half, left
    df <- stretch_to_fill(df, mask_left, upper = TRUE, layer_names)
    # stretch lower half, left
    df <- stretch_to_fill(df, mask_left, upper = FALSE, layer_names)
    # stretch upper half, right
    df <- stretch_to_fill(df, mask_right, upper = TRUE, layer_names)
    # stretch lower half, right
    df <- stretch_to_fill(df, mask_right, upper = FALSE, layer_names)
    
    # Save transformed coordinates
    df[mask,c("x_trans","y_trans")] <- coord_trans
    
    # Make plot
    plot_nonlinear_smoothing <- plot_results(
      df, df[,c("x_bins","y_bins")], mouse_num,
      paste0("Transformed cortical layers, mouse ", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_nonlinear_smoothing = plot_nonlinear_smoothing))
    
    # Estimate layer boundaries
    if (verbose) snk.report...("Estimating layer boundaries")
    layer_boundary_bins <- rep(0, length(layer_names))
    for (lb_num in seq_along(layer_names)) {
      # "layer boundary" will mean the floor of the named layer; the floor of L6b is zero, by definition.
      if (any(as.character(df$layer) == layer_names[lb_num])) {
        minL <- min(df[mask_left & as.character(df$layer) == layer_names[lb_num],"y_bins"], na.rm = TRUE)
        minR <- min(df[mask_right & as.character(df$layer) == layer_names[lb_num],"y_bins"], na.rm = TRUE)
        layer_boundary_bins[lb_num] <- round(mean(minL,minR),0)
      } else {
        layer_boundary_bins[lb_num] <- NA
      }
    }
    # Set boundary for L6b as zero
    layer_boundary_bins[length(layer_names)] <- 0
    
    return(list(df = df, plot_list = plot_list, layer_boundary_bins = layer_boundary_bins))
    
  }

# Helper function to transform coordinates for each mouse and extract layer boundary estimates
cortical_coordinate_transform <- function(
    count_data,               # Data to transform
    total_bins,               # Number of bins into which to divide the spatial axes, e.g., 100
    keep_plots = FALSE,
    L1_removed = RemoveL1,    # Include data labeled as layer 1?
    nat_left = FALSE,         # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    nat_right = FALSE,        # This setting used based on manual visual checks to ensure consistent orientation, causes angle and coordinate sign flips in layer leveling
    verbose = TRUE
  ) {
    
    # Count_data input should have genes across in columns, not stacked in rows
    
    if (verbose) {
      snk.report("Transforming raw data into laminar and columnar coordinates")
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 0)
    }
    
    # Extract list and slice plots 
    slice_plots <- count_data$slice_plots
    count_data <- count_data$count_data
    
    # Set layer names
    layer_names <- c("L1", "L23", "L4", "L5", "L6a", "L6b")
    layer_boundary_bins <- array (0, dim = c(length(unique(count_data$mouse)),length(layer_names)))
    rownames(layer_boundary_bins) <- unique(count_data$mouse)
    plots <- list()
    
    # Loop through mice
    for (mouse_num in unique(count_data$mouse)) {
      if (verbose) snk.report(paste("Mouse number", mouse_num))
      
      # Grab slice plot
      assign(paste0("plot_list_m", mouse_num), list(slice_plot = slice_plots[[paste0("slice_plot", mouse_num)]]))
      
      # Transform and bin coordinates
      count_data <- coordinate_binning(mouse_num, total_bins, layer_names, count_data, L1_removed, nat_left, nat_right)
      
      # Extract layer boundary estimates and plots 
      layer_boundary_bins[mouse_num,] <- count_data$layer_boundary_bins
      assign(paste0("plot_list_m", mouse_num), c(get(paste0("plot_list_m", mouse_num)), count_data$plot_list))
      count_data <- count_data$df
      
      # Make histogram of cell distribution across cortical axis
      this_data <- count_data[count_data$mouse == mouse_num,]
      if (nrow(this_data) > 50000) {
        this_data <- this_data[sample(1:nrow(this_data), 50000),]
      }
      hist_distro_plot <- ggplot(this_data, aes(x=y_bins)) +
        geom_histogram(bins = total_bins, fill = "steelblue", color = "black", na.rm = TRUE) + 
        labs(x = "bin num (cortical axis)", y = "observations per bin", title = "Histograms of Laminar (y-axis) Observation Distribution") +
        theme_minimal()
      assign(paste0("plot_list_m", mouse_num), c(get(paste0("plot_list_m", mouse_num)), list(hist_distro_plot = hist_distro_plot)))
      
      # Save plots
      if (keep_plots) plots <- c(plots, list(get(paste0("plot_list_m", mouse_num))))
      else plots <- NULL
      
      # Combine and print plots of interest
      main_title <- textGrob(paste("Coordinate Transformation, mouse", mouse_num), gp = gpar(fontsize = 20, fontface = "bold"))
      make_plots <- get(paste0("plot_list_m", mouse_num))[c("slice_plot", "plot_nonlinear_smoothing", "hist_distro_plot")]
      make_plots <- do.call(arrangeGrob, c(make_plots, ncol = length(make_plots)))
      make_plots <- arrangeGrob(main_title, make_plots, ncol = 1, heights = c(0.05, 0.95))
      grid.arrange(make_plots)
      
    }
    
    # Check if this is Allen data
    if (!("trscrpt_gene_symb" %in% colnames(count_data))) {
      # ... Allen data formatted differently
      
      # Rearrange count_data columns 
      n <- which(colnames(count_data) == "y_coord") # last column before genes
      m <- which(colnames(count_data) == "x_trans") # first column after genes
      t <- ncol(count_data) # last column
      count_data <- count_data[,c(1:n,m:t,(n+1):(m-1))]
      
      # Make cell_num column on front 
      cell_num <- rownames(count_data)
      count_data <- cbind(cell_num, count_data)
    }
    
    return(
      list(
        df = count_data, 
        layer.boundary.bins = layer_boundary_bins,
        plots = plots
      )
    )
    
  }

# Functions for converting to count_data df for WSPmm modal ############################################################

# Function to convert to WSPmm format
create.count.data.WSPmm <- function(
    df.merfish,                                    # data frame of MERFISH data, produced by cortical_coordinate_transform function
    bin.dim = c("x_bins", "y_bins"),               # dimension by which to bin; must be one of these two; will collapse along the other
    gene.list,                                     # list of genes to include in count data, forms the "child" level of the model
    fixed.effect.names, 
    parent = NULL                                  # parent level for fixed effects; if NULL, will use "cortex"
  ) { 
    
    # Run checks 
    if (length(bin.dim) != 1 || !(bin.dim %in% c("x_bins","y_bins"))) stop("bin.dim must be one of 'x_bins' or 'y_bins'")
    if (!all(c("mouse", "cell_num", fixed.effect.names, parent) %in% colnames(df.merfish))) stop("df.merfish missing mouse, cell_num, parent, or fixed effect column")
    
    # Make count data column names and find its dimensions 
    num_genes <- length(gene.list)
    num_cells <- nrow(df.merfish)
    if (is.null(parent)) parent <- "cortex"
    numrow <- num_cells * num_genes
    
    # Make parent column 
    if (parent == "cortex") parent_col <- rep("cortex", numrow)
    else parent_col <- rep(df.merfish[,parent], num_genes)
    
    # Pre-allocate as much of the count data as possible
    count_data <- data.frame(
      count = rep(as.numeric(NA), numrow),
      bin = rep(df.merfish[,bin.dim], num_genes),
      parent = parent_col,
      gene = rep(gene.list, each = num_cells), 
      mouse = rep(df.merfish[,"mouse"], num_genes),
      cell_num = rep(df.merfish[,"cell_num"], num_genes)
    )
    colnames(count_data)[3] <- parent
    
    # Fill out the fixed-effect columns
    count_data_fe <- data.frame(
      matrix(rep(as.character(NA), numrow * length(fixed.effect.names)), nrow = numrow, ncol = length(fixed.effect.names))
    )
    colnames(count_data_fe) <- fixed.effect.names
    for (fe in fixed.effect.names) {
      count_data_fe[,fe] <- rep(df.merfish[,fe], num_genes)
    }
    
    # Combine into one df 
    count_data <- cbind(count_data, count_data_fe)
    
    # Fill in the count data
    for (i in 1:num_genes) {
      g <- gene.list[i]
      idx_initial <- (i-1)*num_cells+1
      idx_final <- i*num_cells
      count_data[idx_initial:idx_final, "count"] <- df.merfish[, g]
    }
    
    # return data frames and token_counts vector 
    return(count_data)
    
  }


# Function to convert to WSPmm format
create.count.data.WSPmm.allen <- function(
    df.merfish,                                    # data frame of MERFISH data, produced by the above functions
    bin.dim = c("x_bins", "y_bins"),               # dimension by which to bin; must be one of these two; will collapse along the other
    fixed.effect.names, 
    parent = NULL                                  # parent level for fixed effects; if NULL, will use "cortex"
  ) { 
    
    # Make count data column names and find its dimensions 
    num_genes <- length(gene.list)
    if (is.null(parent)) parent <- "cortex"
    numrow <- nrow(df.merfish)
    
    # Make parent column 
    if (parent == "cortex") parent_col <- rep("cortex", numrow)
    else parent_col <- rep(df.merfish[,parent], num_genes)
    
    # Grab data
    count_data <- data.frame(
      count = df.merfish$trscrpt_ct,
      bin = df.merfish[,bin.dim],
      parent = parent_col,
      gene = df.merfish$trscrpt_gene_symb,
      mouse = df.merfish$mouse,
      cell_num = df.merfish$cluster_id
    )
    colnames(count_data)[3] <- parent
    
    # Add age column 
    df.merfish$age <- as.factor("p5x")
    
    # Fill out the fixed-effect columns
    count_data_fe <- data.frame(
      matrix(rep(as.character(NA), numrow * length(fixed.effect.names)), nrow = numrow, ncol = length(fixed.effect.names))
    )
    colnames(count_data_fe) <- fixed.effect.names
    for (fe in fixed.effect.names) {
      count_data_fe[,fe] <- df.merfish[,fe]
    }
    
    # Combine into one df 
    count_data <- cbind(count_data, count_data_fe)
   
    # return data frames and token_counts vector 
    return(count_data)
    
  }

