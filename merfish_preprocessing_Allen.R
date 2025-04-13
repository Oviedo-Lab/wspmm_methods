
# Load library to source Python code snippets
if ( !is.element("package:reticulate", search()) ) { # only try to load if not already loaded
  if(!require(reticulate)) {                         # if not installed, install
    install.packages("reticulate")                   # install if needed
    library(reticulate)                              # load
  }
}

# What version of python are we using? (For grabbing ROI masks)
# ... use_python("/usr/local/bin/python3.12")
use_virtualenv("~/envs/allen-env", required = TRUE)

# Orientation of reference-atlas space: 
# ... * Anterior -> Posterior * Superior -> Inferior * Left -> Right

# Set paths to files
data_path <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/Allen_data.csv"
data_path_h5 <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/C57BL6J-638850-raw.h5ad"
data_path_cellmeta <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/cell_metadata.csv"

make_slice_data <- function(
    slice = 21, 
    gene.list = c("Bcl11b", "Rorb") 
  ) {
    
    # Load sparse-matrix data
    file <- hdf5r::H5File$new(data_path_h5, mode = "r")
    #list.groups(file)
    #list.datasets(file)
    
    # Grab data from H5 file
    gene_codes <- file[["/X"]][["indices"]][]
    cell_id_idx <- file[["/X"]][["indptr"]][]
    cell_id_list <- file[["/obs/cell_label"]][]
    slice_id_list <- file[["/obs/brain_section_label/categories"]][]
    slice_num_list <- file[["/obs/brain_section_label/codes"]][]
    trscrpt_id_list <- file[["/var/transcript_identifier"]][]
    trscrpt_gene_symb_list <- file[["/var/gene_symbol"]][]
    trscrpt_gene_id_list <- file[["/var/gene_identifier"]][]
    run_lengths <- diff(cell_id_idx)
    gene_idx <- gene_codes + 1
    
    # Mask gene_idx for this slice
    gene_idx_mask = rep(slice_num_list, times = run_lengths) == slice
    gene_idx <- gene_idx[gene_idx_mask]
    
    # Mask cell_id for this slice
    slice_mask <- slice_num_list == slice
    slice_num_list <- slice_num_list[slice_mask]
    cell_id_list <- cell_id_list[slice_mask]
    run_lengths <- run_lengths[slice_mask]
    
    # Reorganize masked data into data frame
    num_initial_rows <- sum(gene_idx_mask)
    data <- data.frame(
      trscrpt_ct = file[["/X"]][["data"]][gene_idx_mask],
      trscrpt_gene_symb = trscrpt_gene_symb_list[gene_idx],
      slice_num = rep(slice_num_list, times = run_lengths),
      region = rep(NA, num_initial_rows),
      cell_id = rep(cell_id_list, times = run_lengths),
      slice_id = rep(slice_id_list[slice_num_list + 1], times = run_lengths),
      trscrpt_id = trscrpt_id_list[gene_idx],
      trscrpt_gene_id = trscrpt_gene_id_list[gene_idx]
    )
    
    # Remove blanks
    blank_mask <- grepl("Blank", data$trscrpt_gene_symb)
    data <- data[!blank_mask, ]
    
    # Get cell locations and remove unlabeled cells
    cell_meta <- data.table::fread(data_path_cellmeta)
    data <- data[data$cell_id %in% cell_meta$cell_label, ]
    
    # Keep only genes of interest
    data <- data[data$trscrpt_gene_symb %in% gene.list, ]
    
    # Grab cluster ids 
    data$cluster_id <- cell_meta$cluster_alias[match(data$cell_id, cell_meta$cell_label)]
    
    # Put spatial coordinates into data frame
    data <- merge(
      data,
      cell_meta[, c("cell_label", "x", "y", "z")],
      by.x = "cell_id",
      by.y = "cell_label",
      all.x = TRUE
    )
    
    # Reorganize columns
    data <- data[,c("trscrpt_ct", "trscrpt_gene_symb", "slice_num", "region", "x", "y", "z", "cluster_id", "cell_id", "slice_id", "trscrpt_id", "trscrpt_gene_id")]
    
    # Return data 
    return(data)
    
  }

slice_data <- list()
slice_range <- 30:40
for (s in seq_along(slice_range)) {
  cat("\nSlice data, slice", s)
  slice_data[[s]] <- make_slice_data(slice = slice_range[s])
}

# These dimensions go x, z, y? 
make_masks <- FALSE 
Layer_names <- c("L1", "L23", "L4", "L5", "L6a", "L6b")
if (make_masks) {
  source_python("Allen_CCF.py")
  masks_list <- generate_roi_masks()
  masks_list <- masks
  dim_names <- c("AntPost", "SupInf", "LefRig", "Layer")
  masks <- data.frame()
  for (m in seq_along(masks_list)) {
    masks_list[[m]] <- as.data.frame(masks_list[[m]])
    masks_list[[m]]$Layer <- Layer_names[m]
    colnames(masks_list[[m]]) <- dim_names
    masks <- rbind(masks, masks_list[[m]])
  }
  rm(masks_list)
  write.csv(masks, file = "masks.csv", row.names = FALSE)
} else {
  masks <- read.csv("masks.csv")
}

library(rgl)
mult <- 100
xm <- 1320/mult
ym <- 800/mult
zm <- 1140/mult

to_polar <- function(p) {
  if (is.null(dim(p))) {
    p <- array(p, dim = c(1,length(p)))
  }
  colnames(p) <- c("x", "y")
  r <- sqrt(p$x^2 + p$y^2)
  theta <- atan(abs(p$y/p$x))
  theta[p$x < 0 & p$y > 0] <- pi - theta[p$x < 0 & p$y > 0]
  theta[p$x < 0 & p$y < 0] <- pi + theta[p$x < 0 & p$y < 0]
  theta[p$x > 0 & p$y < 0] <- 2*pi - theta[p$x > 0 & p$y < 0]
  return(data.frame(r, theta))
}
to_cart <- function(p) {
  if (is.null(dim(p))) {
    p <- array(p, dim = c(1,length(p)))
  }
  colnames(p) <- c("r", "theta")
  x <- p$r * cos(p$theta)
  y <- p$r * sin(p$theta)
  return(data.frame(x, y))
}
translate <- function(p, p0) {
  colnames(p) <- c("x", "y")
  names(p0) <- c("x", "y")
  p$x <- p$x - p0["x"]
  p$y <- p$y - p0["y"]
  return(p)
}
twist <- function(p, p0, AL, AR) {
  p <- translate(p, p0)
  pp <- to_polar(p)
  pL <- data.frame(r = pp$r, theta = pp$theta + AL)
  pR <- data.frame(r = pp$r, theta = pp$theta + AR)
  pL <- to_cart(pL)
  pR <- to_cart(pR)
  max_right <- max(p$x)
  max_left <- min(p$x)
  mask_left <- p$x < 0
  mask_right <- !mask_left
  left_ratio <- p$x[mask_left] / max_left
  right_ratio <- p$x[mask_right] / max_right
  p$x[mask_left] <- p$x[mask_left] * (1 - left_ratio) + pL$x[mask_left] * left_ratio
  p$y[mask_left] <- p$y[mask_left] * (1 - left_ratio) + pL$y[mask_left] * left_ratio
  p$x[mask_right] <- p$x[mask_right] * (1 - right_ratio) + pR$x[mask_right] * right_ratio
  p$y[mask_right] <- p$y[mask_right] * (1 - right_ratio) + pR$y[mask_right] * right_ratio
  p <- translate(p, -p0)
  return(p)
}
pinch <- function(p, p0, gammaL, gammaR) {
  p <- translate(p, p0)
  pp <- to_polar(p)
  ppr2 <- pp$r^2
  pL <- data.frame(r = sqrt(ppr2 * gammaL), theta = pp$theta)
  pR <- data.frame(r = sqrt(ppr2 * gammaR), theta = pp$theta)
  pL <- to_cart(pL)
  pR <- to_cart(pR)
  max_right <- max(p$x)
  max_left <- min(p$x)
  mask_left <- p$x < 0
  mask_right <- !mask_left
  left_ratio <- p$x[mask_left] / max_left
  right_ratio <- p$x[mask_right] / max_right
  p$x[mask_left] <- p$x[mask_left] * (1 - left_ratio) + pL$x[mask_left] * left_ratio
  p$y[mask_left] <- p$y[mask_left] * (1 - left_ratio) + pL$y[mask_left] * left_ratio
  p$x[mask_right] <- p$x[mask_right] * (1 - right_ratio) + pR$x[mask_right] * right_ratio
  p$y[mask_right] <- p$y[mask_right] * (1 - right_ratio) + pR$y[mask_right] * right_ratio
  p <- translate(p, -p0)
  return(p)
}



plot_slice <- function(
    slicedata, 
    pinchL,
    pinchR,
    twistL,
    twistR,
    threeD = FALSE,
    slice_num = NULL
    ) {
  
  downsample_ratio <- 1/4
  
  # Find the z coordinates of these slices 
  z_coord <- rep(0, length(slicedata))
  
  
  if (threeD) {
    these_rows <- sample(c(1:nrow(masks)), size = as.integer(nrow(masks)*downsample_ratio/2), replace = FALSE)
    plot3d(
      x = masks$AntPost[these_rows],
      y = c(masks$SupInf + 180)[these_rows],
      z = masks$LefRig[these_rows],
      col = "steelblue",
      ylim = c(0, 1000), xlim = c(0, 1000), zlim = c(0, 1000),
    )
  }
  
  for (s in seq_along(slicedata)) {
    
    if (!is.null(slice_num)) {
      if (slice_range[s] != slice_num) {
        next
      }
    }
    
    cat("\nPlotting slice", slice_range[s])
    
    # Down-sample data
    slicedata_s <- slicedata[[s]][
      seq(
        1, nrow(slicedata[[s]]), 
        length.out = as.integer(nrow(slicedata[[s]])*downsample_ratio)
        ), ]
    
    # Convert slice coordinates into CCF
    slicedata_s$x <- (zm - slicedata_s$x) * mult
    slicedata_s$y <- (slicedata_s$y) * mult
    slicedata_s$z <- (xm - slicedata_s$z) * mult
    
    # Grab the z coordinate of the slice
    z_coord[s] <- as.integer(slicedata_s$z[1])
    
    # Plot slice 
    if (threeD) {
      plot3d(
        x = slicedata_s$z, 
        y = slicedata_s$y,
        z = slicedata_s$x,
        col = "gray",
        add = TRUE,
        aspect = TRUE
      )
    } else {
      
      # Find center 
      center <- c(
        min(slicedata_s$x) + (max(slicedata_s$x) - min(slicedata_s$x))/2, 
        min(slicedata_s$y) + (max(slicedata_s$y) - min(slicedata_s$y))/2
      )
      
      slicedata_s[,c("x", "y")] <- pinch(slicedata_s[,c("x", "y")], center, pinchL, pinchR)
      slicedata_s[,c("x", "y")] <- twist(slicedata_s[,c("x", "y")], center, twistL, twistR)
      
      # Plot background slice
      plot_title <- paste0("Slice ", slice_range[s], ", z = ", z_coord[s])
      plot(slicedata_s$x, slicedata_s$y, col = "gray", pch = 19, cex = 0.5, asp = 1, main = plot_title, ylim = c(0, 1000), xlim = c(0, 1200))
      
      # Grab z and set downsample size
      mask_z <- masks$AntPost == z_coord[s]
      ds_size <- 1000 
      superinf_adj <- 180
      
      # Downsample mask
      ds_array <- array(0, dim = c(ds_size*length(Layer_names), 3))
      for (i in seq_along(Layer_names)) {
        # Grab coordinates in this layer and z slice
        slidemask <- masks$Layer == Layer_names[i] & mask_z
        if (sum(slidemask) == 0) next
        
        # Subset mask to slice
        leftright <- masks[slidemask, "LefRig"]
        superinf <- masks[slidemask, "SupInf"] + superinf_adj
        
        # Down-sample 
        n_slidemask <- sum(slidemask)
        these_rows <- sample(c(1:n_slidemask), size = ds_size, replace = n_slidemask < ds_size)
        row_range <- ((i-1)*ds_size + 1):(i*ds_size)
        
        # Fill array
        ds_array[row_range, 1] <- leftright[these_rows]
        ds_array[row_range, 2] <- superinf[these_rows]
        ds_array[row_range, 3] <- i + 1
      }
      points(ds_array[,1], ds_array[,2], col = ds_array[,3], pch = 19, cex = 0.5)
      
      # Plot center point 
      names(center) <- c("x", "y")
      points(center["x"], center["y"], col = "black", pch = 19, cex = 3)
      
      # Find normal of each slice mask
      # ... In 2D plot, x is LefRig, y is SupInf
      # ... make masks
      left_side <- masks[mask_z, "LefRig"] < center[1]
      layer1 <- masks[mask_z, "Layer"] == "L1"
      layer6 <- masks[mask_z, "Layer"] == "L6a"
      # ... find top and bottom points
      find_topbot <- function(this_side, this_layer) {
        SupInf_ <- c(masks$SupInf[mask_z])[this_side & this_layer] + superinf_adj
        LefRig_ <- c(masks$LefRig[mask_z])[this_side & this_layer]
        Top_idx <- which.max(SupInf_)
        Top_ <- c(
          LefRig_[Top_idx], 
          SupInf_[Top_idx]
        )
        Bot_idx <- which.min(SupInf_)
        Bot_ <- c(
          LefRig_[Bot_idx], 
          SupInf_[Bot_idx]
        )
        out <- cbind(Top_, Bot_)
        colnames(out) <- c("Tp", "Bt")
        rownames(out) <- c("x", "y")
        out <- as.data.frame(out)
        return(out)
      }
      L1TB <- find_topbot(left_side, layer1)
      L6TB <- find_topbot(left_side, layer6)
      R1TB <- find_topbot(!left_side, layer1)
      R6TB <- find_topbot(!left_side, layer6)
      
      find_midpoint <- function(p1, p2) {
        names(p1) <- c("x", "y") 
        names(p2) <- c("x", "y")
        points(p1["x"], p1["y"], col = "black", pch = 19, cex = 1.5)
        points(p2["x"], p2["y"], col = "black", pch = 19, cex = 1.5)
        x_ <- p1["x"] + (p2["x"] - p1["x"])/2 
        y_ <- p1["y"] + (p2["y"] - p1["y"])/2
        p <- c(x_, y_)
        names(p) <- c("x", "y")
        return(p)
      }
      L1mid <- find_midpoint(L1TB$Tp, L1TB$Bt)
      L6mid <- find_midpoint(L6TB$Tp, L6TB$Bt)
      R1mid <- find_midpoint(R1TB$Tp, R1TB$Bt)
      R6mid <- find_midpoint(R6TB$Tp, R6TB$Bt)
      points(L1mid["x"], L1mid["y"], col = "black", pch = 19, cex = 2)
      points(L6mid["x"], L6mid["y"], col = "black", pch = 19, cex = 2)
      points(R1mid["x"], R1mid["y"], col = "black", pch = 19, cex = 2)
      points(R6mid["x"], R6mid["y"], col = "black", pch = 19, cex = 2)
      
      line_slopeInter <- function(origin, point) {
        intercept <- (point["x"]*origin["y"] - point["y"]*origin["x"])/(point["x"] - origin["x"])
        slope <- (point["y"] - intercept)/point["x"]
        param <- c(slope, intercept)
        names(param) <- c("slope", "intercept")
        return(param)
      }
      
      Lnorm_slopeInter <- line_slopeInter(L1mid, L6mid)
      Rnorm_slopeInter <- line_slopeInter(R1mid, R6mid)
      
      Lradius_slopeInter <- line_slopeInter(L1mid, center)
      Rradius_slopeInter <- line_slopeInter(R1mid, center)
      
      make_line <- function(slopeInter, x_range) {
        slope <- slopeInter["slope"]
        intercept <- slopeInter["intercept"]
        y_range <- slope*x_range + intercept
        line_ <- as.data.frame(cbind(x_range, y_range))
        colnames(line_) <- c("x", "y")
        return(line_)
      }
      
      Lnorm <- make_line(Lnorm_slopeInter, 1:1000)
      Rnorm <- make_line(Rnorm_slopeInter, 1:1000)
      
      Lradius <- make_line(Lradius_slopeInter, 1:1000)
      Rradius <- make_line(Rradius_slopeInter, 1:1000)
      
      
      lines(Lnorm$x, Lnorm$y, col = "black", lwd = 3)
      lines(Rnorm$x, Rnorm$y, col = "black", lwd = 3)
      
      lines(Lradius$x, Lradius$y, col = "black", lwd = 3)
      lines(Rradius$x, Rradius$y, col = "black", lwd = 3)
      
      make_circle <- function(origin, r) {
        names(origin) <- c("x", "y")
        theta <- seq(0, 2*pi, length.out=100)
        x <- origin["x"] + r * cos(theta)
        y <- origin["y"] + r * sin(theta)
        circ <- data.frame(x, y)
      }
      vec_dist <- function(p1, p2) {
        names(p1) <- c("x", "y")
        names(p2) <- c("x", "y")
        dist <- sqrt((p1["x"] - p2["x"])^2 + (p1["y"] - p2["y"])^2)
        return(dist)
      }
      # circ_Ltop <- make_circle(center, vec_dist(center, L1TB$Tp))
      # circ_Rtop <- make_circle(center, vec_dist(center, R1TB$Tp))
      # lines(circ_Ltop$x, circ_Ltop$y)
      # lines(circ_Rtop$x, circ_Rtop$y)
      
      
      
    }
    
  }
  
}
plot_slice(slice_data, pinchL = 1, pinchR = 1, twistL = 0, twistR = 0, slice_num = 35)
plot_slice(slice_data, pinchL = 0.95, pinchR = 1, twistL = 0, twistR = 0, slice_num = 35)
plot_slice(slice_data, pinchL = 0.95, pinchR = 1, twistL = -0.25, twistR = 0, slice_num = 35)

# our columnar axis is approx 1mm, these slices are 200um apart, should should be able to get 4 slices 
# Need to estimate the SupInf coordinate of our horizontal slices!!!??

