
# Setup ################################################################################################################

# Install and load DESeq2
if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")}
if (!require("apeglm", quietly = TRUE)) {BiocManager::install("apeglm")}
library(DESeq2) # vignette("DESeq2")

# Preprocess MERFISH data ##############################################################################################

# Load and parse raw visgen data from files
count_data <- make_count_data(data_path_cells, remove_L1 = TRUE, verbose = TRUE)

# Transform coordinates for each mouse into laminar and columnar axes and extract layer boundary estimates
count_data <- cortical_coordinate_transform(
  count_data = count_data, 
  total_bins = 100,        # Number of bins to use when binning data
  keep_plots = TRUE,       # Keep coordinate transformation plots? 
  L1_removed = TRUE,
  nat_right = TRUE,
  verbose = TRUE
)

# Unpack 
count_data <- count_data$df

# Run lm of non-replicated P7 mouse ####################################################################################

# Helper function to pull out the largest effects
weighted_top <- function(
    fe_result,
    top = 10,
    weight = TRUE
  ) {
    
    # Grab effect size and p-value
    effect_size <- fe_result$log2FoldChange
    effect_p <- fe_result$pvalue
    
    # Weight effect size squared by inverse of p-value 
    if (weight) weighted_effect_size <- effect_size^2 / effect_p
    else weighted_effect_size <- abs(effect_size)
    
    # Rank and grab top rows
    ranks <- order(weighted_effect_size, decreasing = TRUE)
    topten <- rownames(fe_result)[ranks[1:top]]
    
    return(topten)
    
  }

# Helper function to put data in proper form for lm
make_P7lm_data <- function(
    cdata
  ) {
    
    # Keep only P7 data
    cdata <- cdata[cdata$age == "7", ]
    
    # Grab variable columns
    mice <- cdata$mouse
    hemi <- cdata$hemisphere
    layer <- cdata$layer
    n_cells <- length(mice)
    
    # Grab gene columns
    last_var_col <- which(colnames(cdata) == "y_bins")
    gene_cols <- (last_var_col+1):ncol(cdata)
    n_genes <- length(gene_cols)
    gene <- colnames(cdata)[gene_cols]
    
    # Make full variable columns 
    mice <- rep(mice, n_genes)
    hemi <- rep(hemi, n_genes)
    layer <- rep(layer, n_genes)
    
    # Make full gene column and count column 
    gene <- rep(gene, each = n_cells)
    count <- rep(0, n_cells * n_genes)
    log2_count <- rep(0, n_cells * n_genes)
    
    # Make shell 
    data <- data.frame(
      mice = mice,
      hemi = hemi,
      layer = layer,
      gene = gene,
      count = count,
      log2_count = log2_count
    )
    
    # Fill count column 
    for (i in gene_cols) {
      data$count[data$gene == colnames(cdata)[i]] <- cdata[, i] 
    }
    
    # Take log plus 1
    data$log2_count <- log(data$count + 1, base = 2)
    
    return(data)
    
  }
p7data <- make_P7lm_data(count_data)

# Helper function to find log2FoldChanges of left > right hemi from lm fit of P7 data
run_P7 <- function(
    p7data
  ) {
   
    # Make data frame shell
    gene_list <- unique(p7data$gene) 
    n_genes <- length(gene_list)
    gene_col <- rep(gene_list, each = 6)
    coeffs_col <- rep(NA, n_genes * 6)
    pvals_col <- rep(NA, n_genes * 6)
    layer_col <- rep(NA, n_genes * 6)
    lm_results <- data.frame(
      genes = gene_col,
      log2FoldChange = coeffs_col,
      pvalue = pvals_col,
      layer = layer_col
    )
    
    # Fill with results
    for (g in gene_list) {
      # Prune data to just gene g
      p7data_g <- p7data[p7data$gene == g, ]
      # Test across entire hemisphere 
      p7model_g <- lm(log2_count ~ hemi, data = p7data_g)
      mask <- grepl("hemiright", names(p7model_g$coefficients))
      coeffs_all <- summary(p7model_g)$coefficients[, "Estimate"][mask]
      pvals_all <- summary(p7model_g)$coefficients[, "Pr(>|t|)"][mask]
      # Test across layers
      p7model_g <- lm(log2_count ~ hemi * layer, data = p7data_g)
      mask <- grepl("hemiright", names(p7model_g$coefficients))
      coeffs_l <- summary(p7model_g)$coefficients[, "Estimate"][mask]
      pvals_l <- summary(p7model_g)$coefficients[, "Pr(>|t|)"][mask]
      # Combine 
      coeffs <- c(coeffs_all, coeffs_l)
      pvals <- c(pvals_all, pvals_l)
      # Get layer names, assuming L23 is reference level of layer
      mask <- grepl("layer", names(coeffs_l))
      layer <- c("all", "L23", gsub("hemiright:", "", names(coeffs_l)[mask]))
      # Load results frame
      mask <- lm_results$genes == g
      lm_results$log2FoldChange[mask] <- coeffs
      lm_results$pvalue[mask] <- pvals
      lm_results$layer[mask] <- layer
    }
    
    # Find top differentially expressed genes
    top_fe_hemi_P7 <- data.frame()
    for (l in unique(lm_results$layer)) {
      lm_results_l <- lm_results[lm_results$layer == l, ]
      top_DE_genes <- weighted_top(lm_results_l, top = 20, weight = FALSE)
      lm_results_l <- lm_results_l[top_DE_genes, ]
      top_fe_hemi_P7 <- rbind(top_fe_hemi_P7, lm_results_l)
    }
    
    # Return results 
    return(top_fe_hemi_P7)
  }
top_fe_hemi_P7 <- run_P7(p7data)

# Save as csv
write.csv(
  top_fe_hemi_P7,
  file = "genes_of_interest_HemiDE_by_age_KEGG_P7.csv",
  row.names = FALSE
)

# Run DESeq2 ###########################################################################################################

# Helper function to put data in proper form for DESeq2
make_SESeq2_data <- function(
    cdata,
    ref_age = "12",
    layer = "all"
  ) {
    
    # Remove P7 and set reference level
    cdata <- cdata[cdata$age != "7", ]
    cdata$age <- relevel(cdata$age, ref = ref_age)
    
    # Grab variable levels
    last_var_col <- which(colnames(cdata) == "y_bins")
    mice <- unique(cdata$mouse)
    hemi <- unique(cdata$hemisphere)
    sample_names <- c(
      paste0(as.character(mice), "_", as.character(hemi)[1]),
      paste0(as.character(mice), "_", as.character(hemi)[2])
    )
    
    # Prune by layer
    if (layer != "all") {
      cdata <- cdata[cdata$layer == layer,]
    }
    
    # Get unique ages
    age <- c()
    for (m in mice) {
      age <- c(age, unique(as.character(cdata[cdata$mouse == m, "age"])))
    }
    
    # Make coldata
    coldata <- data.frame(
      hemi = rep(hemi, each = length(mice)), 
      age = rep(age, length(hemi)),
      mouse = rep(mice, length(hemi))
    )
    coldata$age <- as.factor(coldata$age)
    coldata$age <- relevel(coldata$age, ref = ref_age)
    rownames(coldata) <- sample_names
    
    # Make count matrix 
    gene_cols <- (last_var_col+1):ncol(cdata)
    cts <- array(NA, dim = c(length(gene_cols), nrow(coldata)))
    rownames(cts) <- colnames(cdata)[gene_cols]
    colnames(cts) <- rownames(coldata)
    for (i in gene_cols) {
      for (j in rownames(coldata)) {
        i_ <- i - last_var_col
        ct_mask <- cdata$mouse == coldata[j,"mouse"] & cdata$hemisphere == coldata[j,"hemi"]
        cts[i_,j] <- sum(cdata[ct_mask, i])
      }
    }
    
    return(list(cts = cts, coldata = coldata))
    
  }

# Helper function to run DESeq2
run_DESeq2 <- function(
    cdata,
    ref_age = "12",
    layer = "all",
    controls = TRUE
  ) {
    
    # Make data for DESeq2
    ddata <- make_SESeq2_data(cdata, ref_age, layer)
    cts <- ddata$cts
    coldata <- ddata$coldata
    
    # Construct a DESeqDataSet: 
    dds <- DESeqDataSetFromMatrix(
      countData = cts,
      colData = coldata,
      design = ~ age*hemi
    )
    
    # Run differential expression analysis, age and hemisphere
    gene_names <- rownames(cts)
    if (controls) dds <- estimateSizeFactors(dds, controlGenes = which(gene_names %in% c("Tubb2a", "Tuba1b", "Actl9", "Actbl2")))
    else dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds, fitType='local')
    
    # Apply shrinkage (ridge regression)
    resLFC_hemi <- lfcShrink(dds, coef = "hemi_right_vs_left", type = "apeglm")
    resLFC_hemi_results <- as.data.frame(resLFC_hemi@listData)
    
    return(
      list(
        cts = cts, 
        coldata = coldata,  
        resLFC_hemi = resLFC_hemi_results
      )
    )
    
  }

# Run DESeq2 on all layers
run_DESeq2_alllayers <- function(
    count_data
  ) {
    DESeq2_results <- list()
    layers <- c("all", unique(as.character(count_data$layer)))
    ages <- unique(as.character(count_data$age))
    controls <- c("ctrl", "no_ctrl")
    for (c in controls) {
      if (c == "ctrl") next
      DESeq2_results_c <- list()
      for (a in ages) {
        if (a == "7") next
        DESeq2_results_a <- list()
        for (l in layers) {
          DESeq2_results_a[[l]] <- run_DESeq2(
            cdata = count_data, 
            ref_age = a,
            layer = l,
            controls = "ctrl" == c
          )
        }
        DESeq2_results_c[[a]] <- DESeq2_results_a
      }
      DESeq2_results[[c]] <- DESeq2_results_c
    }
    return(DESeq2_results)
  }
DESeq2_results <- run_DESeq2_alllayers(count_data)

# Extract and analyze results ##########################################################################################

# Compile top differentially expressed genes from each layer
find_top_fe <- function(
    DESeq2_results,
    count_data
  ) {
    top_fe_hemi <- data.frame()
    layers <- c("all", unique(as.character(count_data$layer)))
    ages <- unique(as.character(count_data$age))
    controls <- c("ctrl", "no_ctrl")
    for (c in controls) {
      if (c == "ctrl") next
      top_fe_hemi_c <- data.frame()
      DESeq2_results_c <- DESeq2_results[[c]]
      for (a in ages) {
        if (a == "7") next
        top_fe_hemi_a <- data.frame()
        DESeq2_results_a <- DESeq2_results_c[[a]]
        for (l in layers) {
          fixed_effect_hemi <- DESeq2_results_a[[l]]$resLFC_hemi
          genes <- weighted_top(fixed_effect_hemi, top = 20, weight = FALSE)
          top_fe_hemi_l <- fixed_effect_hemi[genes,]
          layer <- rep(l, nrow(top_fe_hemi_l))
          ref_age <- rep(a, nrow(top_fe_hemi_l))
          controls <- rep(c, nrow(top_fe_hemi_l))
          top_fe_hemi_l <- cbind(
            genes,
            top_fe_hemi_l,
            layer,
            ref_age,
            controls
          )
          top_fe_hemi_a <- rbind(top_fe_hemi_a, top_fe_hemi_l)
        }
        # Label layer with largest effect for each gene
        top_fe_hemi_a$max <- rep(FALSE, nrow(top_fe_hemi_a))
        for (g in unique(top_fe_hemi_a$genes)) {
          mask <- top_fe_hemi_a$genes == g
          max_fc <- max(abs(top_fe_hemi_a$log2FoldChange[mask]))
          mask <- mask & abs(top_fe_hemi_a$log2FoldChange) == max_fc
          top_fe_hemi_a$max[mask] <- TRUE
        }
        top_fe_hemi_c <- rbind(top_fe_hemi_c, top_fe_hemi_a)
      }
      top_fe_hemi <- rbind(top_fe_hemi, top_fe_hemi_c)
    }
    rownames(top_fe_hemi) <- NULL
    return(top_fe_hemi)
  }
top_fe_hemi <- find_top_fe(DESeq2_results, count_data)

# Save as csv
write.csv(
  top_fe_hemi,
  file = "genes_of_interest_HemiDE_by_age_KEGG.csv",
  row.names = FALSE
)

# Clean up 
rm(count_data, p7data, DESeq2_results)

# Compare top genes across age and control #############################################################################

compare_results <- function() {
    
    p7genes <- unique(top_fe_hemi_P7$genes) # no controls, linear model
    gene_lists <- list()
    gene_lists[["p12genes_c"]] <- unique(top_fe_hemi[top_fe_hemi$ref_age == "12" & top_fe_hemi$controls == "ctrl", "genes"]) 
    gene_lists[["p12genes_nc"]] <- unique(top_fe_hemi[top_fe_hemi$ref_age == "12" & top_fe_hemi$controls == "no_ctrl", "genes"])
    gene_lists[["p18genes_c"]] <- unique(top_fe_hemi[top_fe_hemi$ref_age == "18" & top_fe_hemi$controls == "ctrl", "genes"])
    gene_lists[["p18genes_nc"]] <- unique(top_fe_hemi[top_fe_hemi$ref_age == "18" & top_fe_hemi$controls == "no_ctrl", "genes"])
    cat("\nDE across hemispheres, P7 without controls (lm):\n")
    print(sort(p7genes))
    cat("\nDE across hemispheres, P12 with controls:\n")
    print(sort(gene_lists[["p12genes_c"]]))
    cat("\nDE across hemispheres, P12 without controls:\n")
    print(sort(gene_lists[["p12genes_nc"]]))
    cat("\nDE across hemispheres, P18 with controls:\n")
    print(sort(gene_lists[["p18genes_c"]]))
    cat("\nDE across hemispheres, P18 without controls:\n")
    print(sort(gene_lists[["p18genes_nc"]]))
    
    compare_lists <- function(set1, set2) {
      # Common genes, P12 vs P18, with controls
      cnt <- 0 
      common <- c()
      for (g in gene_lists[[set1]]) {
        if (any(g == gene_lists[[set2]])) {
          cnt <- cnt + 1
          common <- c(common, g)
        }
      }
      cat("\n\nCommon genes,", set1, set2, "\n")
      print(sort(common))
      cat("\nNumber of genes in common:", cnt)
      cat("\nNumber of genes,", set1, length(gene_lists[[set1]]))
      cat("\nNumber of genes,", set2, length(gene_lists[[set2]]))
    }
    
    ages <- c("p12", "p18")
    ctrls <- c("genes_c", "genes_nc")
    paste0(rep(ages, length(ctrls)), rep(ctrls, each = length(ages)))
    
    # For each age, compare with and without controls
    for (a in ages) {
      compare_lists(paste0(a, ctrls[1]), paste0(a, ctrls[2]))
    }
    
    # For each ctrls, compare each age to the others 
    for (c in ctrls) {
      hit <- array(dim = c(0,2))
      for (a1 in ages) {
        for (a2 in ages) {
          if (a1 == a2) next
          this_pair <- c(a1, a2)
          repeat_pair <- FALSE
          if (nrow(hit) > 0) {
            for (r in 1:nrow(hit)) {
              if (hit[r,1] == this_pair[2] && hit[r,2] == this_pair[1]) {
                repeat_pair <- TRUE
              }
            }
          }
          if (repeat_pair) next
          compare_lists(paste0(a1, c), paste0(a2, c))
          hit <- rbind(hit, this_pair)
        }
      }
    }
    
  }
compare_results()




