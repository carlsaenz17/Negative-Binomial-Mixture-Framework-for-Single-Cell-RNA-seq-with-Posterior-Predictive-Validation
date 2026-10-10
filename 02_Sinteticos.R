set.seed(1710)
##Creacion de la simulacion.
#Configuracion inicial
n_cells <- 4000
n_genes <- 15000

#Proporciones de los clusters
cluster_sizes <- c(
  900, 800, 800, 800,
  300, 300,
  50, 30,
  10, 5, 3, 2
)

stopifnot(sum(cluster_sizes) == n_cells)

group_prob <- cluster_sizes / sum(cluster_sizes)

#Generacion de los datos
simulate_dataset <- function(dropout = TRUE, de_prob = 0.15) {

  params <- newSplatParams()

  params <- setParams(params,
    nGenes = n_genes,
    batchCells = n_cells,

#Library size (ground truth)
    lib.loc = 11,
    lib.scale = 0.4,

#Mean expression
    mean.shape = 0.2,#antes 0.3
    mean.rate = 0.3,#antes 0.1

#Dropout
    dropout.type = ifelse(dropout, "experiment", "none"),
    dropout.mid = 3.5, #Antes era 0 y daba NB sparsity: 0.5679561
    dropout.shape = -1,

#Grupos
    group.prob = group_prob,

#Genes
    de.prob = de_prob,
    de.downProb = 0.5,
    de.facLoc = 1.0,
    de.facScale = 0.5
  )

  sim <- splatSimulate(
    params,
    method = "groups",
    verbose = FALSE
  )
  counts <- counts(sim)
  labels <- sim$Group

#Hacer clusters raros
  rare_clusters <- which(cluster_sizes <= 10)

  for (rc in rare_clusters) {

    cells_idx <- which(labels == rc)

    marker_genes <- sample(1:n_genes, 200)

    counts[marker_genes, cells_idx] <-
      counts[marker_genes, cells_idx] +
      matrix(
        rpois(length(cells_idx) * length(marker_genes), lambda = 15),
        nrow = length(marker_genes)
      )
  }

  list(
    counts = counts,
    labels = labels,
    size_factors = sim$ExpLibSize,  # ground truth real
    sparsity = mean(counts == 0)
  )
}

##Simulacion del NB
cat("Simulación NB...\n")
sim_nb <- simulate_dataset(dropout = FALSE)
cat("NB sparsity:", sim_nb$sparsity, "\n")

##Simulacion del ZINB
cat("Simulación ZINB...\n")
sim_zinb <- simulate_dataset(dropout = TRUE)
cat("ZINB sparsity:", sim_zinb$sparsity, "\n")

##Simulacion del NULL
cat("Simulación NULL...\n")
params_null <- newSplatParams()
params_null <- setParams(params_null,
  nGenes = n_genes,
  batchCells = n_cells,

  lib.loc = 11,
  lib.scale = 0.4,

  mean.shape = 0.3,
  mean.rate = 0.1,

  dropout.type = "none",

  group.prob = 1,     # un solo cluster
  de.prob = 0
)

sim_null_obj <- splatSimulate(
  params_null,
  method = "groups",
  verbose = FALSE
)

sim_null <- list(
  counts = counts(sim_null_obj),
  labels = rep(1, n_cells),
  size_factors = sim_null_obj$ExpLibSize,
  sparsity = mean(counts(sim_null_obj) == 0)
)
cat("NULL sparsity:", sim_null$sparsity, "\n")

##Asignacion para poder usar la matriz de conteos.
cells <- SingleCellExperiment(
list(counts = sim_nb$counts)  #Para el modelo NB
#list(counts = sim_zinb$counts) #PAra el modelo ZINB
#list(counts = sim_null$counts) #PAra el modelo NULL
  )

colData(cells)$true_group <- sim_nb$labels
#colData(cells)$true_group <- sim_zinb$labels
#colData(cells)$true_group <- sim_null$labels

#Size-factors
clusters <- quickCluster(cells)
cells <- computeSumFactors(cells, clusters = clusters)

sf <- sizeFactors(cells)

sf[is.na(sf)] <- 1
sf[sf <= 0] <- 1

#Creacion de la matriz de conteos
matrix_counts <- t(counts(cells))  # células × genes
