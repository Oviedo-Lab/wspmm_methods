

data_path <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/Allen_data.csv"
data_path_h5 <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/C57BL6J-638850-raw.h5ad"
data_path_cellmeta <- "/Users/michaelbarkasi/Library/CloudStorage/OneDrive-WashingtonUniversityinSt.Louis/projects_Oviedo_lab/MERFISH/development_work/data/Allen_data/cell_metadata.csv"

# start_time <- Sys.time()
# Allen_data_raw <- data.table::fread(data_path)
# duration <- Sys.time() - start_time
# cat("Time to read in data: ", duration, units(duration), "\n")
# Allen_data_raw_slice25 <- Allen_data_raw[Allen_data_raw$slice_num == 25, ]
# head(Allen_data_raw)


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

slice <- 25
gene_idx_mask = rep(slice_num_list, times = run_lengths) == slice
gene_idx <- gene_idx[gene_idx_mask]

slice_mask <- slice_num_list == slice
slice_num_list <- slice_num_list[slice_mask]
cell_id_list <- cell_id_list[slice_mask]
run_lengths <- run_lengths[slice_mask]

# Reorganize into data frame
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

blank_mask <- grepl("Blank", data$trscrpt_gene_symb)
data <- data[!blank_mask, ]

cell_meta <- data.table::fread(data_path_cellmeta)
locatable_cells <- unique(cell_meta$cell_label)
data <- data[data$cell_id %in% locatable_cells, ]
gene.list <- c("Bcl11b", "Fezf2", "Rorb", "Satb2", "Nxph3", "Cux2", "Rorb") 
data <- data[data$trscrpt_gene_symb %in% gene.list, ]

data <- merge(
  data,
  cell_meta[, c("cell_label", "x", "y", "z")],
  by.x = "cell_id",
  by.y = "cell_label",
  all.x = TRUE
)

data <- data[,c("trscrpt_ct", "trscrpt_gene_symb", "slice_num", "region", "x", "y", "z", "cell_id", "slice_id", "trscrpt_id", "trscrpt_gene_id")]


plot(data$x, data$y, col = as.integer(factor(data$trscrpt_gene_symb)), pch = 19, cex = 0.5)


