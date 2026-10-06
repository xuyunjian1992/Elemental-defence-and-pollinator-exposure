# ============================================================
# Microbial co-occurrence network analysis
# Cleaned submission version based on the original analysis code
# Groups: Control (CK) and Cd treatment (Cd)
# ============================================================

# Required packages -------------------------------------------------------
library(microeco)
library(ggplot2)
library(rgexf)
library(igraph)

# Reproducibility ---------------------------------------------------------
set.seed(20240812)
theme_set(theme_bw())

# Input files -------------------------------------------------------------
# Run this script from the folder containing the "data" directory.
tax <- read.csv("data/Group1.csv", row.names = 1)

# Network analysis function ----------------------------------------------
# The analytical settings below are unchanged from the original code:
#   detection frequency > 20%
#   summed relative abundance >= 0.002
#   Spearman correlation
#   P < 0.05
#   |rho| >= 0.6
#   FDR adjustment
#   fast-greedy community detection

run_network_analysis <- function(otu_file, group_name) {

  # Read ASV absolute-abundance table
  ASV <- read.csv(otu_file, row.names = 1)

  # Filter ASVs detected in >20% of samples
  ASVs_detection <- ASV
  ASVs_detection[ASVs_detection > 0] <- 1
  Reserve_ASVs <- which(
    rowSums(ASVs_detection) / ncol(ASVs_detection) > 0.2
  )
  ASV3 <- ASV[Reserve_ASVs, , drop = FALSE]

  # Convert absolute abundance to relative abundance
  row_n <- as.numeric(length(row.names(ASV3)))
  col_n <- as.numeric(length(colnames(ASV3)))
  col_sum <- rep(colSums(ASV3), row_n)
  col_sum <- matrix(col_sum, nrow = col_n, ncol = row_n)
  ASV_relative <- ASV3 / t(col_sum)

  # Check column sums after relative-abundance conversion
  print(colSums(ASV_relative))

  # Retain ASVs with summed relative abundance >= 0.002 (0.2%)
  Total <- rowSums(ASV_relative)
  ASV_Shift <- which(Total >= 0.002)
  ASV1 <- ASV3[ASV_Shift, , drop = FALSE]
  ASV2 <- as.matrix(ASV1)

  # Prepare OTU table
  OTU <- data.frame(ASV2)

  # Create microeco dataset
  dataset <- microtable$new(
    otu_table = OTU,
    tax_table = tax
  )
  dataset$tidy_dataset()
  print(dataset)

  # Construct Spearman co-occurrence network
  network <- trans_network$new(
    dataset = dataset,
    cor_method = "spearman",
    use_WGCNA_pearson_spearman = TRUE
  )

  network$cal_network(
    COR_p_thres = 0.05,
    COR_cut = 0.6,
    COR_p_adjust = "fdr",
    add_taxa_name = c("phylum"),
    delete_unlinked_nodes = TRUE,
    COR_optimization = FALSE
  )

  # Module detection
  network$cal_module(method = "cluster_fast_greedy")

  # Generate node table without Zi/Pi node-role calculation
  network$get_node_table(node_roles = FALSE)

  # Module eigengenes
  network$cal_eigen()

  # Network topological attributes
  network$cal_network_attr()

  # Edge and adjacency-matrix information
  network$get_edge_table()
  network$get_adjacency_matrix()

  # Export network and result tables
  dir.create("output", showWarnings = FALSE)

  network$save_network(
    filepath = file.path("output", paste0(group_name, "_network.gexf"))
  )

  write.csv(
    network$res_node_table,
    file.path("output", paste0(group_name, "_node_table.csv")),
    row.names = FALSE
  )

  write.csv(
    network$res_edge_table,
    file.path("output", paste0(group_name, "_edge_table.csv")),
    row.names = FALSE
  )

  write.csv(
    as.matrix(network$res_adjacency_matrix),
    file.path("output", paste0(group_name, "_adjacency_matrix.csv"))
  )

  if (!is.null(network$res_network_attr)) {
    write.csv(
      as.data.frame(network$res_network_attr),
      file.path("output", paste0(group_name, "_network_attributes.csv")),
      row.names = TRUE
    )
  }

  # # Export node-role plot
  # pdf(
  #   file.path("output", paste0(group_name, "_taxa_roles.pdf")),
  #   width = 7,
  #   height = 6
  # )
  # print(network$plot_taxa_roles(use_type = 1))
  # dev.off()

  # Save complete R object for reproducibility
  saveRDS(
    network,
    file.path("output", paste0(group_name, "_network_object.rds"))
  )

  return(network)
}

# Run the two treatment groups -------------------------------------------
network_CK <- run_network_analysis(
  otu_file = "data/CKotu.csv",
  group_name = "CK"
)

network_Cd <- run_network_analysis(
  otu_file = "data/cdotu.csv",
  group_name = "Cd"
)

# End of script ----------------------------------------------------------
