
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
if (RemoveL1) {
  snk.report...("Setting flag to remove Layer 1 from the analysis")
} else {
  snk.report...("Setting flag to keep Layer 1 in the analysis")
}

# Load raw data ########################################################################################################

# Helper function to load and parse data
parse_hdf5 <- function(
    file_path,
    mouse_num = 0,
    raw = TRUE,
    ACx_only = TRUE,
    remove_L1 = RemoveL1
  ) {
    
    # Load data
    file <- H5File$new(file_path, mode = "r")
    
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
    if (raw) transcript_counts_raw <- t(file[["/raw/X"]][,])    # should have integer elements, not normalized
    else transcript_counts_raw <- t(file[["/X"]][,])            # should have integer elements, normalized for MapMyCells
    transcript_counts_raw <- transcript_counts_raw[,nonblanks]  # drop blanks
    colnames(transcript_counts_raw) <- gene_names               
    n_cells <- nrow(transcript_counts_raw)
    
    # Grab rows for cortical layers
    ROI_names <- list(
      LL1 = "/obs/ROI__Left Primary auditory area, layer 1",
      LL23 = "/obs/ROI__Left Primary auditory area, layer 2/3",
      LL4 = "/obs/ROI__Left Primary auditory area, layer 4",
      LL5 = "/obs/ROI__Left Primary auditory area, layer 5",
      LL6a = "/obs/ROI__Left Primary auditory area, layer 6a",
      LL6b = "/obs/ROI__Left Primary auditory area, layer 6b",
      RL1 = "/obs/ROI__Right Primary auditory area, layer 1",
      RL23 = "/obs/ROI__Right Primary auditory area, layer 2/3",
      RL4 = "/obs/ROI__Right Primary auditory area, layer 4",
      RL5 = "/obs/ROI__Right Primary auditory area, layer 5",
      RL6a = "/obs/ROI__Right Primary auditory area, layer 6a",
      RL6b = "/obs/ROI__Right Primary auditory area, layer 6b"
      )
    L1_idx <- file[[ROI_names[["LL1"]]]][] | file[[ROI_names[["RL1"]]]][]
    L23_idx <- file[[ROI_names[["LL23"]]]][] | file[[ROI_names[["RL23"]]]][]
    L4_idx <- file[[ROI_names[["LL4"]]]][] | file[[ROI_names[["RL4"]]]][]
    L5_idx <- file[[ROI_names[["LL5"]]]][] | file[[ROI_names[["RL5"]]]][]
    L6a_idx <- file[[ROI_names[["LL6a"]]]][] | file[[ROI_names[["RL6a"]]]][]
    L6b_idx <- file[[ROI_names[["LL6b"]]]][] | file[[ROI_names[["RL6b"]]]][]
    layer <- rep("notACx",n_cells)
    layer[L1_idx] <- "L1"
    layer[L23_idx] <- "L23"
    layer[L4_idx] <- "L4"
    layer[L5_idx] <- "L5"
    layer[L6a_idx] <- "L6a"
    layer[L6b_idx] <- "L6b"
    layer <- as.factor(layer)
    
    # Grab rows for left and right hemisphere
    left_idx <- file[[ROI_names[["LL1"]]]][] | 
      file[[ROI_names[["LL23"]]]][] | 
      file[[ROI_names[["LL4"]]]][] | 
      file[[ROI_names[["LL5"]]]][] | 
      file[[ROI_names[["LL6a"]]]][] | 
      file[[ROI_names[["LL6b"]]]][]
    right_idx <- file[[ROI_names[["RL1"]]]][] | 
      file[[ROI_names[["RL23"]]]][] | 
      file[[ROI_names[["RL4"]]]][] | 
      file[[ROI_names[["RL5"]]]][] | 
      file[[ROI_names[["RL6a"]]]][] | 
      file[[ROI_names[["RL6b"]]]][]
    hemisphere <- rep("notACx", n_cells)
    hemisphere[left_idx] <- "left"
    hemisphere[right_idx] <- "right"
    hemisphere <- as.factor(hemisphere)
    
    # Grab metadata and form columns
    age <- as.factor(rep(as.integer(sub("P", "", file[["/uns/info/metadata/age"]][])), n_cells))
    sex <- as.factor(rep(file[["/uns/info/metadata/sex"]][], n_cells))
    strain <- as.factor(rep(file[["/uns/info/metadata/strain"]][], n_cells))
    experience <- as.factor(rep(file[["/uns/info/metadata/treatment"]][], n_cells))
    
    # Grab spatial coordinates and form columns 
    # ... units are in microns (um)
    x_max <- file[["/obs/max_x"]][] 
    x_min <- file[["/obs/min_x"]][]
    x_cen <- file[["/obs/center_x"]][]
    y_max <- file[["/obs/max_y"]][] 
    y_min <- file[["/obs/min_y"]][]
    y_cen <- file[["/obs/center_y"]][]
    
    # Assign a number to the mouse
    mouse <- as.factor(rep(mouse_num,n_cells))
    
    # Make data frame
    which_cells <- TRUE 
    if (ACx_only) which_cells <- hemisphere != "notACx"
    transcript_counts <- data.frame(
      mouse, 
      celltype_MMC, cellsubclass_MMC, cellsupertype_MMC, 
      hemisphere, layer, 
      age, sex, 
      strain, experience, 
      x_max, x_min, x_cen, y_max, y_min, y_cen,
      transcript_counts_raw 
    )[which_cells,]
    
    # Replace cell type number with cell type name
    transcript_counts$celltype_MMC <- celltype_names[transcript_counts$celltype_MMC]
    transcript_counts$cellsubclass_MMC <- cellsubclass_names[transcript_counts$cellsubclass_MMC]
    transcript_counts$cellsupertype_MMC <- cellsupertype_names[transcript_counts$cellsupertype_MMC]
    
    # Collapse GABA cell types together 
    GABA_idx <- which(transcript_counts$celltype_MMC == "CTX-CGE GABA" | transcript_counts$celltype_MMC == "CTX-MGE GABA")
    transcript_counts$celltype_MMC[GABA_idx] <- "CTX-CGE/MGE GABA"
    
    # Convert to factor 
    transcript_counts$celltype_MMC <- as.factor(transcript_counts$celltype_MMC)
    
    # Remove Layer 1
    if (remove_L1) {
      L1_idx <- transcript_counts$layer == "L1"
      transcript_counts <- transcript_counts[!L1_idx,]
    }
    
    # Release file
    file$close_all()
    
    return(transcript_counts)
    
  }

# Load and parse data
make_count_data <- function(
    data_path, 
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
    for (f in seq_along(files)) {
      if (f < length(files)) cat(f, ", ", sep="")
      else cat(f, "\n")
      assign(paste0("mouse", f), parse_hdf5(file_path = files[f], mouse_num = f, raw = TRUE))
      ind_var_fields <- which(colnames(get(paste0("mouse", f))) %in% c(
        "mouse", 
        "celltype_MMC", "cellsubclass_MMC", "cellsupertype_MMC", 
        "hemisphere", "layer", 
        "age", "sex", 
        "strain", "experience",
        "x_max", "x_min", "x_cen", "y_max", "y_min", "y_cen"))
      noise_columns <- grep("_noise$", colnames(get(paste0("mouse", f))), value = FALSE)
      mean_rates <- c(mean_rates, mean(as.matrix(get(paste0("mouse", f))[,-c(ind_var_fields,noise_columns)])))
      cells_per_mouse <- c(cells_per_mouse, nrow(get(paste0("mouse", f))))
    }
    
    # Check the mean of each run to see if there's some reason to suspect a systematic difference
    if (verbose) {
      snk.print_vec("Mean transcript counts per cell for each mouse", mean_rates)
      snk.print_vec("Mean of means", mean(mean_rates))
      snk.print_vec("Standard deviation of means", sd(mean_rates))
    }
    
    # Combine data from all mice
    count_data <- mouse1
    for (f in 2:length(files)) {
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
    
    return(count_data)
    
  }

# Functions to transform into laminar and columnar coordinates #########################################################

# Define helper function for transforming coordinates
coordinate_transform <- function(
    mouse_num, 
    df,
    verbose = TRUE
  ) {
    
    if (verbose) snk.report...("Grabbing coordinates and defining layers and hemispheres")
    
    # Grab indexes
    idx <- df$mouse == mouse_num 
    idx_left <- df[idx,"hemisphere"] == "left"
    idx_right <- df[idx,"hemisphere"] == "right"
    
    # Define columns sets
    all_coord <- c("x_cen","y_cen")
    x_coord <- c("x_cen")
    y_coord <- c("y_cen")
    
    # Grab coordinates 
    coordinates <- df[idx,all_coord]
    
    # Define layer rows
    idx_left_L1 <- idx_left & df[idx,"layer"] == "L1"
    idx_right_L1 <- idx_right & df[idx,"layer"] == "L1"
    idx_left_L23 <- idx_left & df[idx,"layer"] == "L23"
    idx_right_L23 <- idx_right & df[idx,"layer"] == "L23"
    idx_left_L4 <- idx_left & df[idx,"layer"] == "L4"
    idx_right_L4 <- idx_right & df[idx,"layer"] == "L4"
    idx_left_L5 <- idx_left & df[idx,"layer"] == "L5"
    idx_right_L5 <- idx_right & df[idx,"layer"] == "L5"
    idx_left_L6a <- idx_left & df[idx,"layer"] == "L6a"
    idx_right_L6a <- idx_right & df[idx,"layer"] == "L6a"
    idx_left_L6b <- idx_left & df[idx,"layer"] == "L6b"
    idx_right_L6b <- idx_right & df[idx,"layer"] == "L6b"
    
    # Put layer rows in convenient lists
    # ... left hemisphere
    layer_rows_left <- list(
      L1 = idx_left_L1,
      L23 = idx_left_L23,
      L4 = idx_left_L4,
      L5 = idx_left_L5,
      L6a = idx_left_L6a,
      L6b = idx_left_L6b
    )
    # ... right hemisphere
    layer_rows_right <- list(
      L1 = idx_right_L1,
      L23 = idx_right_L23,
      L4 = idx_right_L4,
      L5 = idx_right_L5,
      L6a = idx_right_L6a,
      L6b = idx_right_L6b
    )
    
    # Plot original data
    range_limit <- range(c(df$x_cen, df$y_cen), na.rm = TRUE)
    plot_untransformed <- ggplot(df[df$mouse == mouse_num,], aes(x = x_cen, y = y_cen, color = layer)) +
      geom_point(na.rm = TRUE) +
      labs(
        x = "X Coordinate", 
        y = "Y Coordinate", 
        color = "Layer", 
        title = paste("Untransformed ACx layers, mouse", mouse_num)
        ) +
      theme_minimal() +
      scale_x_continuous(limits = range_limit) + 
      scale_y_continuous(limits = range_limit)
    plot_list <- list(plot_untransformed = plot_untransformed)
    
    # Helper function for future plotting
    plot_results <- function(
      df, coordinates, 
      idx, idx_right, idx_left, x_coord,
      plot_title,
      centered = TRUE
    ) {
      
      layer_right <- df$layer[idx][idx_right]
      df_right <- coordinates[idx_right,]
      if (centered) df_right[,x_coord] <- df_right[,x_coord] + abs(min(df_right[,x_coord])) + 50
      df_right <- cbind(df_right, layer_right)
      colnames(df_right)[ncol(df_right)] <- "layer"
      layer_left <- df$layer[idx][idx_left]
      df_left <- coordinates[idx_left,]
      if (centered) df_left[,x_coord] <- df_left[,x_coord] - abs(max(df_left[,x_coord])) - 50
      df_left <- cbind(df_left, layer_left)
      colnames(df_left)[ncol(df_left)] <- "layer"
      layers_both <- rbind(df_right, df_left)
      
      plot <- ggplot(layers_both, aes(x = x_cen, y = y_cen, color = layer)) +
        geom_point(na.rm = TRUE) +
        labs(x = "X Coordinate", y = "Y Coordinate", color = "Layer", title = plot_title) +
        theme_minimal()
      
      if (centered) {
        plot <- plot +
          ylim(-1500,1500) +
          xlim(-1500,1500)
      } else {
        plot <- plot + 
          scale_x_continuous(limits = range_limit - max(range_limit)/2) + 
          scale_y_continuous(limits = range_limit)
      }
      
      return(plot)
      
    }
    
    # Step 1: Ensure we're "facing" the slice from a consistent perspective (right on the right)
    if (verbose) snk.report...("Step 1, correcting perspective to ensure right hemisphere is on the right")
    correct_perspective <- function(coord) {
      
      # Center the points around the y-axis
      x_translation <- mean(coord[,x_coord])
      # Subtract from all x-coordinates 
      coord[,x_coord] <- coord[,x_coord] - x_translation
      
      # If mean right coordinate is less than mean left coordinate, reflect across midline 
      mean_right_x <- mean(coord[idx_right, x_coord])
      mean_left_x <- mean(coord[idx_left, x_coord])
      if (mean_right_x < mean_left_x) reflection_value <- -1
      else reflection_value <- 1
      # Multiply x-coordinates to reflect across midline
      coord[,x_coord] <- coord[,x_coord] * reflection_value
      
      return(coord)
      
    }
    coordinates <- correct_perspective(coordinates)
    
    # Test by plotting (corrected perspective)
    plot_perspective_correction <- plot_results(
      df, coordinates, 
      idx, idx_right, idx_left, x_coord,
      paste("Untransformed ACx layers (corrected perspective), mouse", mouse_num),
      centered = FALSE
    )
    plot_list <- c(plot_list, list(plot_perspective_correction = plot_perspective_correction))
    
    # Step 2: Center each patch (left and right) around the mean point of L5
    if (verbose) snk.report...("Step 2, centering each patch around the mean point of L5")
    recenter_around_L5 <- function(coord, idx_hemisphere, idx_layer) {
      
      mean_x <- mean(coord[idx_layer, x_coord])
      mean_y <- mean(coord[idx_layer, y_coord])
      
      coord[idx_hemisphere,x_coord] <- coord[idx_hemisphere,x_coord] - mean_x
      coord[idx_hemisphere,y_coord] <- coord[idx_hemisphere,y_coord] - mean_y
      
      return(coord)
      
    }
    coordinates <- recenter_around_L5(coordinates, idx_right, idx_right_L5)
    coordinates <- recenter_around_L5(coordinates, idx_left, idx_left_L5)
    
    # Test by plotting (recentered)
    plot_recenter <- plot_results(
      df, coordinates, 
      idx, idx_right, idx_left, x_coord,
      paste("Untransformed ACx layers (recentered), mouse", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_recenter = plot_recenter))
    
    # Helper functions for rotational linear transformations and leveling
    
    rot_matrix <- function(
      theta
    ) {
      return(matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, byrow = TRUE))
    }
    
    linear_transform <- function(
      coord, 
      matrix, 
      y_adj = 1
    ) {
      coord_cen <- cbind(coord[, c(x_coord, y_coord)])
      transformed <- as.matrix(coord_cen) %*% t(matrix)
      coord[, c(x_coord, y_coord)] <- transformed
      coord[,y_coord] <- coord[,y_coord] * y_adj
      return(coord)
    }
    
    level_layer <- function(
      data, 
      idx_hemisphere, 
      idx_layer, 
      verbose = FALSE, 
      reverse_angle = FALSE
    ) {
      # Fit linear model to the idx_layer coordinates to get the angle of the slice
      model <- lm(y_cen ~ x_cen, data = data[idx_layer,])
      # Find angle
      angle <- atan(model$coefficients[x_coord])
      if (reverse_angle) {
        angle <- pi + angle # should never be needed, but leaving in for now? 
      }
      if (verbose) cat("\nangle: ", angle*57.3)
      # level
      new_coord <- linear_transform(data[idx_hemisphere,], rot_matrix(angle))
      return(new_coord)
    }
    
    # Step 3: Rotate each patch so that L4 aligns with the x-axis and L1 (or L23, if L1 removed) is on top
    if (verbose) snk.report...("Step 3, rotating each patch so that L4 aligns with the x-axis and L1 (or L23, if L1 removed) is on top")
    level_around_lvl <- function(
      mouse_num, coord, 
      idx_hemisphere_right, idx_hemisphere_left,
      idx_lvl_layer_right, idx_lvl_layer_left,
      idx_ref_layer_right, idx_ref_layer_left
    ) {
     
      # level given layer
      coord[idx_hemisphere_right,] <- level_layer(coord, idx_hemisphere_right, idx_lvl_layer_right, verbose = FALSE)
      coord[idx_hemisphere_left,] <- level_layer(coord, idx_hemisphere_left, idx_lvl_layer_left, verbose = FALSE)
      
      # Check orientation, right 
      y_mean_layer_right <- mean(coord[idx_lvl_layer_right,y_coord])
      y_mean_ref_right <- mean(coord[idx_ref_layer_right,y_coord])
      if (y_mean_ref_right < y_mean_layer_right) coord[idx_hemisphere_right,y_coord] <- -coord[idx_hemisphere_right,y_coord]
      
      # Check orientation, left
      y_mean_layer_left <- mean(coord[idx_lvl_layer_left,y_coord])
      y_mean_ref_left <- mean(coord[idx_ref_layer_left,y_coord])
      if (y_mean_ref_left < y_mean_layer_left) coord[idx_hemisphere_left,y_coord] <- -coord[idx_hemisphere_left,y_coord]
      
      return(coord)
      
    }
    coordinates <- level_around_lvl(mouse_num, coordinates, idx_right, idx_left, idx_right_L4, idx_left_L4, idx_right_L23, idx_left_L23)
    
    # Test by plotting (L4 leveled)
    plot_level_L4 <- plot_results(
      df, coordinates, 
      idx, idx_right, idx_left, x_coord,
      paste("Transformed ACx layers (L4 leveled), mouse", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_level_L4 = plot_level_L4))
    
    # Step 4: Linearly transform each patch to vertically "straighten" them while preserving area
    if (verbose) snk.report...("Step 4, linearly transforming each patch to vertically straighten them while preserving area")
    
    # Helper functions for linear minimization of skew
    
    layer_x_min_diff <- function(
      df, 
      idx
    ) {
      x_mins <- c()
      for (layer in names(idx)) {
        if (any(idx[[layer]])) x_mins <- c(x_mins, min(df[idx[[layer]],x_coord]))
      }
      return(abs(min(x_mins)-max(x_mins)))
    }
    
    layer_x_max_diff <- function(
      df, 
      idx
    ) {
      x_maxs <- c()
      for (layer in names(idx)) {
        if (any(idx[[layer]])) x_maxs <- c(x_maxs, max(df[idx[[layer]],x_coord]))
      }
      return(abs(min(x_maxs)-max(x_maxs)))
    }
    
    # Define loss function
    transform_loss <- function(
      parameters, 
      coord, 
      idx_hemisphere, 
      idx_layers, 
      initial_height, 
      initial_width
    ) {
      
      # Convert parameter vector into matrix 
      y_adj <- parameters[5]
      parameters <- matrix(parameters[1:4], nrow = 2, byrow = TRUE) 
      
      # Apply the transformation to the coordinates
      coord[idx_hemisphere,] <- linear_transform(coord[idx_hemisphere,], parameters, y_adj)
      
      # Find current height and width 
      current_height <- max(coord[idx_hemisphere,y_coord]) - min(coord[idx_hemisphere,y_coord])
      current_width <- max(coord[idx_hemisphere,x_coord]) - min(coord[idx_hemisphere,x_coord])
      
      # Find distance from initial height and width
      distance <- sqrt((current_height - initial_height)^2 + (current_width - initial_width)^2)
      
      # Find and return loss
      return(
        distance + 
          layer_x_min_diff(coord, idx_layers) + 
          layer_x_max_diff(coord, idx_layers)
      )
      
    }
    
    # Find initial height and width 
    initial_height_left <- max(coordinates[idx_left,y_coord]) - min(coordinates[idx_left,y_coord])
    initial_width_left <- max(coordinates[idx_left,x_coord]) - min(coordinates[idx_left,x_coord])
    initial_height_right <- max(coordinates[idx_right,y_coord]) - min(coordinates[idx_right,y_coord])
    initial_width_right <- max(coordinates[idx_right,x_coord]) - min(coordinates[idx_right,x_coord])
    
    # Find linear transformation which minimizes loss for each hemisphere, by x-distance ("skew")
    parameters_left <- optim(
      par = c(1,0,0,1,1), 
      fn = transform_loss, 
      coord = coordinates, 
      idx_hemisphere = idx_left, 
      idx_layers = layer_rows_left,
      initial_height = initial_height_left,
      initial_width = initial_width_left
    )$par
    parameters_right <- optim(
      par = c(1,0,0,1,1), 
      fn = transform_loss, 
      coord = coordinates, 
      idx_hemisphere = idx_right, 
      idx_layers = layer_rows_right,
      initial_height = initial_height_right,
      initial_width = initial_width_right
    )$par
    linear_left <- matrix(parameters_left[1:4], nrow = 2, byrow = TRUE)
    linear_right <- matrix(parameters_right[1:4], nrow = 2, byrow = TRUE)
    y_adj_left <- parameters_left[5]
    y_adj_right <- parameters_right[5]
    
    # Apply the linear transformation to the coordinates
    coordinates[idx_left,] <- linear_transform(coordinates[idx_left,], linear_left, y_adj_left)
    coordinates[idx_right,] <- linear_transform(coordinates[idx_right,], linear_right, y_adj_right)
    
    # Test by plotting 
    plot_linear_skew <- plot_results(
      df, coordinates, 
      idx, idx_right, idx_left, x_coord,
      paste("Transformed ACx layers (linear skew), mouse", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_linear_skew = plot_linear_skew))
    
    # Step 5: Level all layers
    if (verbose) snk.report...("Step 5, leveling all layers")
    
    # Find each point's y-distance to the mean y-point of each layer
    left_transforms <- array(0, dim = c(nrow(coordinates),2))
    right_transforms <- array(0, dim = c(nrow(coordinates),2))
    for (layer in names(layer_rows_left)) { # same names for each hemisphere
      if (any(layer_rows_left[[layer]])) {
        left_mean <- mean(coordinates[layer_rows_left[[layer]],y_coord])
        right_mean <- mean(coordinates[layer_rows_right[[layer]],y_coord])
        left_distances <- abs(coordinates[layer_rows_left[[layer]],y_coord] - left_mean)
        right_distances <- abs(coordinates[layer_rows_right[[layer]],y_coord] - right_mean)
        left_distances <- (max(left_distances) - left_distances) / max(left_distances)
        right_distances <- (max(right_distances) - right_distances) / max(right_distances)
        left_transforms[layer_rows_left[[layer]],] <- as.matrix(coordinates[layer_rows_left[[layer]],]) * (1-left_distances) + 
          as.matrix(level_layer(coordinates, idx_left, layer_rows_left[[layer]]))[layer_rows_left[[layer]][idx_left],] * left_distances
        right_transforms[layer_rows_right[[layer]],] <- as.matrix(coordinates[layer_rows_right[[layer]],]) * (1-right_distances) + 
          as.matrix(level_layer(coordinates, idx_right, layer_rows_right[[layer]]))[layer_rows_right[[layer]][idx_right],] * right_distances
      }
    }
    coordinates[idx_left,] <- left_transforms[idx_left,]
    coordinates[idx_right,] <- right_transforms[idx_right,]
    
    # Test by plotting 
    plot_level_all <- plot_results(
      df, coordinates, 
      idx, idx_right, idx_left, x_coord,
      paste0("Transformed ACx layers (level all layers), mouse ", mouse_num)
    )
    plot_list <- c(plot_list, list(plot_level_all = plot_level_all))
    
    # Translate so all points are positive 
    coordinates[,x_coord] <- coordinates[,x_coord] - min(coordinates[,x_coord])
    coordinates[,y_coord] <- coordinates[,y_coord] - min(coordinates[,y_coord])
    
    return(list(coord = coordinates, plot_list = plot_list))
    
  }

# Define helper function for binning coordinates and smoothing edges
coordinate_binning <- function(
    mouse_num,
    total_bins, 
    layer_names,
    df,
    L1_removed = RemoveL1,
    verbose = TRUE
  ) {
    
    # Grab rows for the mouse
    idx <- df$mouse == mouse_num 
    idx_left <- idx & df[,"hemisphere"] == "left"
    idx_right <- idx & df[,"hemisphere"] == "right"
    
    # Apply coordinate transform to those rows
    if (verbose) snk.report...("Performing coordinate transform")
    coord_trans <- coordinate_transform(mouse_num, df)
    plot_list <- coord_trans$plot_list
    coord_trans <- coord_trans$coord
    
    # Find bin boundaries
    max_x <- max(coord_trans[,1])
    max_y <- max(coord_trans[,2])
    bin_ticks_x <- seq(0,max_x,length.out = total_bins+1)
    bin_ticks_y <- seq(0,max_y,length.out = total_bins+1)
    
    # Bin the transformed coordinates
    if (verbose) snk.report...("Binning transformed coordinates")
    x_bins <- c() 
    y_bins <- c() 
    for (r in 1:nrow(coord_trans)) {
      x_bins <- c(x_bins, sum(coord_trans[r,1] > bin_ticks_x))
      y_bins <- c(y_bins, sum(coord_trans[r,2] > bin_ticks_y))
    }
    x_bins[x_bins == 0] <- 1
    y_bins[y_bins == 0] <- 1
    
    # Add bin numbers to df
    df[idx,"x_bins_raw"] <- x_bins
    df[idx,"y_bins_raw"] <- y_bins
    df[idx,"x_bins"] <- x_bins
    df[idx,"y_bins"] <- y_bins
    
    # Apply nonlinear smoothing to bins 
    if (verbose) snk.report...("Smoothing bin edges with nonlinear transformation")
    stretch_to_fill <- function(
      df_, 
      hemisphere_idx,
      upper,
      layer_names
    ) {
      
      # Grab layer to lead stretch
      if (upper && L1_removed) layer_num <- 2
      else if (upper && !L1_removed) layer_num <- 1
      else layer_num <- length(layer_names)
      
      # Grab index
      idx_ <- hemisphere_idx & as.character(df_$layer) == layer_names[layer_num] 
      
      # Find max/min y point in layer
      # ... For smoothing laminar edges:
      Ly_by_x <- rep(NA,total_bins)
      Ly_by_x_abs <- rep(NA,total_bins)
      # ... For smoothing columnar edges:
      Lx_by_y <- rep(NA,total_bins)
      if (upper) {
        for (b in 1:total_bins) {
          # Max point in column in upper layer
          if (sum(idx_ & df_$x_bins_raw == b) > 0) Ly_by_x[b] <- max(df_$y_bins_raw[idx_ & df_$x_bins_raw == b], na.rm = TRUE)
          # Max point in column
          if (sum(hemisphere_idx & df_$x_bins_raw == b) > 0) Ly_by_x_abs[b] <- max(df_$y_bins_raw[hemisphere_idx & df_$x_bins_raw == b], na.rm = TRUE)
          # Max point in layer
          if (sum(hemisphere_idx & df_$y_bins_raw == b) > 0) Lx_by_y[b] <- max(df_$x_bins_raw[hemisphere_idx & df_$y_bins_raw == b], na.rm = TRUE)
        }
      } else {
        for (b in 1:total_bins) {
          # Min point in column in lower layer
          if (sum(idx_ & df_$x_bins_raw == b) > 0) Ly_by_x[b] <- min(df_$y_bins_raw[idx_ & df_$x_bins_raw == b], na.rm = TRUE)
          # Min point in column
          if (sum(hemisphere_idx & df_$x_bins_raw == b) > 0) Ly_by_x_abs[b] <- min(df_$y_bins_raw[hemisphere_idx & df_$x_bins_raw == b], na.rm = TRUE)
          # Min point in layer
          if (sum(hemisphere_idx & df_$y_bins_raw == b) > 0) Lx_by_y[b] <- min(df_$x_bins_raw[hemisphere_idx & df_$y_bins_raw == b], na.rm = TRUE)
        }
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
        idx_range <- hemisphere_idx & df_$y_bins_raw > total_bins/2
        idx_range_layer <- hemisphere_idx & df_$x_bins_raw > total_bins/2
        to_ <- c(total_bins/2,total_bins)
      } else {
        idx_range <- hemisphere_idx & df_$y_bins_raw < total_bins/2
        idx_range_layer <- hemisphere_idx & df_$x_bins_raw < total_bins/2
        to_ <- c(1,total_bins/2+1)
      }
      for (b in 1:total_bins) {
        # ... columns
        if (upper) from_ <- c(total_bins/2,Ly_by_x[b]) 
        else from_ <- c(Ly_by_x[b],total_bins/2) 
        df_$y_bins[idx_range & df_$x_bins_raw == b] <- as.integer(
          scales::rescale(
            df_$y_bins_raw[idx_range & df_$x_bins_raw == b], 
            to = to_,
            from = from_
          )
        )
        # ... layers
        if (upper) from_ <- c(total_bins/2,Lx_by_y[b]) 
        else from_ <- c(Lx_by_y[b],total_bins/2)
        df_$x_bins[idx_range_layer & df_$y_bins_raw == b] <- as.integer(
          scales::rescale(
            df_$x_bins_raw[idx_range_layer & df_$y_bins_raw == b], 
            to = to_,
            from = from_
          )
        )
      }
     
      return(df_)
      
    }
   
    # stretch upper half, left
    df <- stretch_to_fill(df, idx_left, upper = TRUE, layer_names)
    # stretch lower half, left
    df <- stretch_to_fill(df, idx_left, upper = FALSE, layer_names)
    # stretch upper half, right
    df <- stretch_to_fill(df, idx_right, upper = TRUE, layer_names)
    # stretch lower half, right
    df <- stretch_to_fill(df, idx_right, upper = FALSE, layer_names)
    
    # Save transformed coordinates
    df[idx,c("x_trans","y_trans")] <- coord_trans
    
    # Make data frames for plotting
    layer_right <- df$layer[idx_right]
    df_right <- df[idx_right,c("x_bins","y_bins")]
    df_right[,"x_bins"] <- df_right[,"x_bins"] + abs(min(df_right[,"x_bins"])) + 5
    df_right <- cbind(df_right, layer_right)
    colnames(df_right)[ncol(df_right)] <- "layer"
    layer_left <- df$layer[idx_left]
    df_left <- df[idx_left,c("x_bins","y_bins")]
    df_left[,"x_bins"] <- df_left[,"x_bins"] - abs(max(df_left[,"x_bins"])) - 5
    df_left <- cbind(df_left, layer_left)
    colnames(df_left)[ncol(df_left)] <- "layer"
    layers_both <- rbind(df_right, df_left)
    
    # Make plot
    plot_nonlinear_smoothing <- ggplot(layers_both, aes(x = x_bins, y = y_bins, color = layer)) +
      geom_point(na.rm = TRUE) +
      ylim(-75,175) +
      xlim(-125,125) +
      labs(x = "X Bin", y = "Y Bin", color = "Layer", title = paste0("Transformed ACx layers, mouse ", mouse_num)) +
      theme_minimal()
    plot_list <- c(plot_list, list(plot_nonlinear_smoothing = plot_nonlinear_smoothing))
    
    # Estimate layer boundaries
    if (verbose) snk.report...("Estimating layer boundaries")
    layer_boundary_bins <- rep(0, length(layer_names))
    for (lb_num in seq_along(layer_names)) {
      # "layer boundary" will mean the floor of the named layer; the floor of L6b is zero, by definition.
      if (any(as.character(df$layer) == layer_names[lb_num])) {
        minL <- min(df[idx_left & as.character(df$layer) == layer_names[lb_num],"y_bins"], na.rm = TRUE)
        minR <- min(df[idx_right & as.character(df$layer) == layer_names[lb_num],"y_bins"], na.rm = TRUE)
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
transform_and_extract <- function(
    count_data,
    total_bins,
    keep_plots = FALSE,
    verbose = TRUE
  ) {
    
    if (verbose) {
      snk.report("Transforming raw data into laminar and columnar coordinates")
      snk.horizontal_rule(reps = snk.simple_break_reps, end_breaks = 0)
    }
    
    layer_names <- c("L1", "L23", "L4", "L5", "L6a", "L6b")
    layer_boundary_bins <- array (0, dim = c(length(unique(count_data$mouse)),length(layer_names)))
    rownames(layer_boundary_bins) <- unique(count_data$mouse)
    plots <- list()
    
    for (mouse_num in unique(count_data$mouse)) {
      if (verbose) snk.report(paste("Mouse number", mouse_num))
      
      # Transform and bin coordinates
      count_data <- coordinate_binning(mouse_num, total_bins, layer_names, count_data)
      
      # Extract layer boundary estimates and plots 
      layer_boundary_bins[mouse_num,] <- count_data$layer_boundary_bins
      assign(paste0("plot_list_m",mouse_num), count_data$plot_list)
      count_data <- count_data$df
      
      # Make histogram of cell distribution across laminar axis
      hist_resids_plot <- ggplot(count_data[count_data$mouse == mouse_num,], aes(x=y_bins)) +
        geom_histogram(bins = total_bins, fill = "steelblue", color = "black", na.rm = TRUE) + 
        labs(x = "bin num (laminar axis)", y = "cells per bin", title = "Histograms of Laminar (y-axis) Cell Distribution") +
        theme_minimal()
      assign(paste0("plot_list_m",mouse_num), c(get(paste0("plot_list_m",mouse_num)), list(hist_resids_plot = hist_resids_plot)))
      
      # Save plots
      if (keep_plots) plots <- c(plots, list(get(paste0("plot_list_m",mouse_num))))
      else plots <- NULL
      
      # Combine and print plots of interest
      main_title <- textGrob(paste("Coordinate Transformation, mouse", mouse_num), gp = gpar(fontsize = 20, fontface = "bold"))
      make_plots <- get(paste0("plot_list_m",mouse_num))[c("plot_recenter", "plot_nonlinear_smoothing", "hist_resids_plot")]
      make_plots <- do.call(arrangeGrob, c(make_plots, ncol = length(make_plots)))
      make_plots <- arrangeGrob(main_title, make_plots, ncol = 1, heights = c(0.05, 0.95))
      grid.arrange(make_plots)
      
    }
    
    # Rearrange count_data columns 
    n <- which(colnames(count_data) == "y_cen") # last column before genes
    m <- which(colnames(count_data) == "x_trans") # first column after genes
    t <- ncol(count_data) # last column
    count_data <- count_data[,c(1:n,m:t,(n+1):(m-1))]
    
    # Make cell_num column on front 
    cell_num <- rownames(count_data)
    count_data <- cbind(cell_num, count_data)
    
    return(
      list(
        df = count_data, 
        layer_boundary_bins = layer_boundary_bins,
        plots = plots
      )
    )
    
  }

# Functions for converting to count_data df for WSPmm modal ############################################################

# Function to convert to WSPmm format
create.count.data.WSPmm <- function(
    df.merfish,                                    # data frame of MERFISH data, produced by transform_and_extract function
    bin.dim = c("x_bins", "y_bins"),               # dimension by which to bin; must be one of these two; will collapse along the other
    gene.list,                                     # list of genes to include in count data, forms the "child" level of the model
    fixed.effect.names, 
    parent = NULL,                                 # parent level for fixed effects; if NULL, will use "cortex"
    verbose = FALSE
  ) { 
    
    # Run checks 
    if (verbose) snk.report...("Running checks")
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
    
    # Combine into one df and add column names
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


