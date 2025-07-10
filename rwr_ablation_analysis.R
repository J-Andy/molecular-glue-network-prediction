# R Script: Multilayer vs. PPI‐Only RWR Ablation (Protein‐only Ranking)
# Description:
#   This script demonstrates a minimal Random Walk with Restart (RWR) ablation analysis
#   to assess the impact of annotation layers (GO, Reactome) on the prioritization of
#   neo-substrate candidates (e.g., RBM5) in a multilayer network. It compares the full
#   multilayer network against a PPI-only network, then ranks proteins based solely on the proteome nodes.

# Required libraries
# install.packages(c("data.table","Matrix"))  # run once if missing
library(data.table)
library(Matrix)

# ---- USER SETTINGS ----
# Base project directory (adjust to your local path)
base_dir    <- "/molecular-glue-network-prediction-main"

# File paths for network layers
ppi_dir       <- file.path(base_dir, "multiplex/1")
react_file    <- file.path(base_dir, "multiplex/2/ReactomePathwaysRelation.tsv")
go_file       <- file.path(base_dir, "multiplex/3/GO.all.tsv")
pg_file       <- file.path(base_dir, "bipartite/1_3.tsv")  # protein→GO
pr_file       <- file.path(base_dir, "bipartite/1_2.tsv")  # protein→Reactome

# RWR parameters and seeds
target       <- "P52756"  # RBM5 UniProt ID
seeds        <- c("Q66K64","Q9BW61","Q16531","P62877","Q13620","Q13619")
alpha        <- 0.3         # restart probability
tol          <- 1e-10       # convergence threshold
max_iter     <- 1000        # maximum iterations

# ---- FUNCTION: Load edge list ----
load_edges <- function(path, is_directory=TRUE) {
  if (is_directory) {
    files <- list.files(path, full.names=TRUE)
    dt <- rbindlist(lapply(files, function(f) {
      tmp <- fread(f, header=FALSE)
      setnames(tmp, c("from","to")); tmp
    }))
  } else {
    dt <- fread(path, header=FALSE)
    setnames(dt, c("from","to"))
  }
  # Clean: remove NA/empty
  dt <- dt[complete.cases(dt) & nchar(from)>0 & nchar(to)>0]
  return(dt)
}

# Load layer data
ppi_dt   <- load_edges(ppi_dir, TRUE)
react_dt <- load_edges(react_file, FALSE)
go_dt    <- load_edges(go_file, FALSE)
pg_dt    <- load_edges(pg_file, FALSE)
pr_dt    <- load_edges(pr_file, FALSE)

# ---- FUNCTION: Random Walk with Restart ----
run_rwr <- function(edge_list, seeds, alpha, tol, max_iter) {
  # Build local node set
  nodes <- unique(c(edge_list$from, edge_list$to, seeds))
  idx   <- setNames(seq_along(nodes), nodes)
  N     <- length(nodes)
  # Construct adjacency matrix
  i <- idx[edge_list$from]; j <- idx[edge_list$to]
  M <- sparseMatrix(i=i, j=j, x=1, dims=c(N,N))
  M <- M + t(M)  # make undirected
  # Row-normalize to transition matrix
  Dinv <- Diagonal(x=1/rowSums(M))
  W    <- Dinv %*% M
  # Initialize probability
  p0 <- numeric(N); p0[idx[seeds]] <- 1/length(seeds)
  p  <- p0
  # RWR iterations
  for (k in seq_len(max_iter)) {
    p_new <- (1-alpha)*(W %*% p) + alpha*p0
    if (sum(abs(p_new - p)) < tol) break
    p <- p_new
  }
  # Return ranking
  res <- data.table(node=nodes, score=as.numeric(p_new))
  setorder(res, -score)
  return(res)
}

# ---- RUN: Full Multilayer vs. PPI-Only ----
full_edges <- rbindlist(list(ppi_dt, react_dt, go_dt, pg_dt, pr_dt))
r_full     <- run_rwr(full_edges, seeds, alpha, tol, max_iter)
r_ppi      <- run_rwr(ppi_dt,    seeds, alpha, tol, max_iter)

# ---- Filter to  proteins and compute ranks ----n#
# Define protein universe (from proteomics: adjust if needed)
proteins <- unique(c(ppi_dt$from, ppi_dt$to, pg_dt$from, pr_dt$from))

compute_rank <- function(ranking, proteins, target, topk=30, seeds) {
  prot_only <- ranking[node %in% proteins]
  recall    <- sum(sapply(seeds, function(s) which(prot_only$node==s) <= topk))
  trank     <- which(prot_only$node == target)
  return(list(recall=recall, rank=trank))
}

res_full <- compute_rank(r_full, proteins, target, 30, seeds)
res_ppi  <- compute_rank(r_ppi,  proteins, target, 30, seeds)

# ---- OUTPUT RESULTS ----
cat("=== Full Multilayer ===\n",
    sprintf("Recall@30 (seeds): %d/6\n", res_full$recall),
    sprintf("RBM5 rank:          %d\n\n", res_full$rank))
cat("=== PPI-Only Ablation ===\n",
    sprintf("Recall@30 (seeds): %d/6\n", res_ppi$recall),
    sprintf("RBM5 rank:          %d\n", res_ppi$rank))
