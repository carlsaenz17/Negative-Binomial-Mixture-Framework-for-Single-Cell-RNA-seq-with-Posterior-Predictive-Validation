#Semilla de esta seccion
set.seed(1710)
# Estadisticas base
counts_mat <- counts(cells)
mean_per_gene <- rowMeans(counts_mat)
var_per_gene  <- rowVars(counts_mat)

## Residuos de Pearson
mu_gene <- rowSums(counts_mat) / sum(sf)
mu_ij <- mu_gene %o% sf
df_nb_fit <- data.frame(
  mu  = mean_per_gene,
  var = var_per_gene
)
df_nb_fit <- df_nb_fit[df_nb_fit$mu > 0.1, ]
fit_alpha <- lm(I(var - mu) ~ I(mu^2) - 1,
                data = df_nb_fit)
alpha_global <- max(coef(fit_alpha)[1], 0)

#Residuos de pearson
var_nb_ij <- mu_ij + alpha_global * mu_ij^2
pearson_mat <- (as.matrix(counts_mat) - mu_ij) /
               sqrt(var_nb_ij + 1e-8)
var_pearson <- rowVars(pearson_mat)
hvg_pearson <- names(sort(var_pearson,
                          decreasing = TRUE)[1:2000])
pca_pearson <- prcomp(t(pearson_mat[hvg_pearson, ]),
                      center = TRUE,
                      scale. = FALSE)

##Seleccion de la mejor dimenson latente
cat("Estimando dimensión latente por 5-fold CV...\n")
mat_input <- t(pearson_mat[hvg_pearson, ])
mat_input <- scale(mat_input, center = TRUE, scale = FALSE)
n_cells <- nrow(mat_input)
max_k <- 50
folds <- createFolds(1:n_cells, k = 5)
calc_rmse_test <- function(pca_obj, X_test, k) {
  scores_test <- predict(pca_obj, newdata = X_test)[, 1:k, drop = FALSE]
  reconstructed <- scores_test %*% t(pca_obj$rotation[, 1:k, drop = FALSE])
  error <- X_test - reconstructed
  sqrt(mean(error^2))
}
rmse_cv <- matrix(NA, nrow = 5, ncol = max_k)
for (i in 1:5) {
  test_idx  <- folds[[i]]
  train_idx <- setdiff(1:n_cells, test_idx)
  X_train <- mat_input[train_idx, ]
  X_test  <- mat_input[test_idx, ]
  pca_train <- prcomp(X_train, center = FALSE, scale. = FALSE)
  for (k in 1:max_k) {
    rmse_cv[i, k] <- calc_rmse_test(pca_train, X_test, k)
  }
}
rmse_mean <- colMeans(rmse_cv)
best_k <- which.min(rmse_mean)

#Criterio parsimonioso (1% del mínimo)
threshold <- min(rmse_mean) * 1.01
best_k_parsimonious <- which(rmse_mean <= threshold)[1]
cat("Dimensión óptima (mínimo RMSE):", best_k, "\n")
cat("Dimensión parsimoniosa (1% regla):", best_k_parsimonious, "\n")

#PCA final con la mejor dimension
pca_final <- prcomp(mat_input,
                    center = FALSE,
                    scale. = FALSE)
reduced_dims <- pca_final$x[, 1:best_k_parsimonious]

cat("Pipeline completo ejecutado correctamente.\n")

## Clustering Leiden
#Aplicamos Leiden sobre reduced_dims
vector_clusters <- clusterRows(
  reduced_dims,
  NNGraphParam(
    k = 20,
    cluster.fun = "leiden",
    cluster.args = list(resolution = 0.5)
  )
)

#Ver el resultado
table(vector_clusters)

## Asignaciones
cells_leiden <- aggregateAcrossCells(cells, ids = vector_clusters, statistics = "mean")

#Extraemos la matriz de promedios (Clusters en filas, Genes en columnas)
mu_leiden <- t(assays(cells_leiden)$counts)
rownames(mu_leiden) <- paste0("C", rownames(mu_leiden))
colnames(mu_leiden) <- rownames(cells)

#Labels por célula
labels <- factor(colData(cells)$true_group)

#Tabla cluster vs label
tab <- table(vector_clusters, labels)

#Etiqueta dominante por cluster
cluster_label <- apply(tab, 1, function(x) {
  names(which.max(x))
})

#Pureza del cluster
cluster_purity <- apply(tab, 1, function(x) {
  max(x) / sum(x)
})

#Tamaños de cluster
cluster_counts <- table(vector_clusters)
n_total <- sum(cluster_counts)

#Tabla final
df_clusters <- data.frame(
  Cluster = names(cluster_counts),
  Etiqueta = cluster_label,
  Celulas = as.numeric(cluster_counts),
  Porcentaje = round(as.numeric(cluster_counts) / n_total * 100, 2),
  Pureza = round(cluster_purity, 2)
)

#Ordenar por abundancia
df_clusters <- df_clusters[order(df_clusters$Celulas, decreasing = TRUE), ]

cat("\nIdentificación de Clusters Leiden (Basado en labels reales):\n")
print(df_clusters, row.names = FALSE)
