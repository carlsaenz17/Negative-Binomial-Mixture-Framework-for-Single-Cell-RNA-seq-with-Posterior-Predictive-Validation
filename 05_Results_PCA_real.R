#Extraemos la matriz de promedios (Clusters en filas, Genes en columnas)
mu_leiden <- t(assays(cells_leiden)$counts)
rownames(mu_leiden) <- paste0("C", rownames(mu_leiden))
colnames(mu_leiden) <- rownames(cells)

##Esta seccion se puede omitir cuando el dataset NO tiene labels
#Aplicamos transformación logarítmica y Z-score (Normalización por gen), esto permite que los genes de la firma que son específicos resalten más (para datos como PBMC4K).
mu_leiden_log <- log1p(mu_leiden)
z_scores <- scale(mu_leiden_log) # Centra y escala cada columna (gen)

#Matriz de score
score_matrix_PCA <- matrix(0,
                           nrow = nrow(mu_leiden),
                           ncol = length(cell_group_ensembe))
rownames(score_matrix_PCA) <- rownames(mu_leiden)
colnames(score_matrix_PCA) <- names(cell_group_ensembe)

#Z-score
for (i in seq_len(nrow(mu_leiden))) {
  for (tipo in names(cell_group_ensembe)) {
    genes <- intersect(cell_group_ensembe[[tipo]], colnames(mu_leiden))
    if (length(genes) > 0) {
      # El score es el promedio de Z-scores de los genes de la firma en ese cluster
      score_matrix_PCA[i, tipo] <- mean(z_scores[i, genes])
    }
  }
}

#Construir tabla
cluster_counts <- table(vector_clusters)
n_total <- sum(cluster_counts)
best_match_idx <- apply(score_matrix_PCA, 1, which.max)

df_clusters <- data.frame(
  Cluster = rownames(score_matrix_PCA),
  Etiqueta = colnames(score_matrix_PCA)[best_match_idx],
  Celulas = as.numeric(cluster_counts),
  Porcentaje = round(as.numeric(cluster_counts) / n_total * 100, 2),
  Score_Z = round(apply(score_matrix_PCA, 1, max), 2)
)

#Ordenar por abundancia
df_clusters <- df_clusters[order(df_clusters$Celulas, decreasing = TRUE), ]
cat("\nIdentificación de Clusters Leiden (Basado en Z-scores):\n")
print(df_clusters, row.names = FALSE)

##Esta es la seccion para cuando los labels si existen:
# #Labels por célula (descomentar desde aqui)
# labels <- factor(cells$label)
# #labels <- factor(cells$true_group)
#
# #Tabla cluster vs label
# tab <- table(vector_clusters, labels)
#
# #Etiqueta dominante por cluster
# cluster_label <- apply(tab, 1, function(x) {
#   names(which.max(x))
# })
#
# #Pureza del cluster
# cluster_purity <- apply(tab, 1, function(x) {
#   max(x) / sum(x)
# })
#
# #Tamaños de cluster
# cluster_counts <- table(vector_clusters)
# n_total <- sum(cluster_counts)
#
# #Tabla final
# df_clusters <- data.frame(
#   Cluster = names(cluster_counts),
#   Etiqueta = cluster_label,
#   Celulas = as.numeric(cluster_counts),
#   Porcentaje = round(as.numeric(cluster_counts) / n_total * 100, 2),
#   Pureza = round(cluster_purity, 2)
# )
#
# #Ordenar por abundancia
# df_clusters <- df_clusters[order(df_clusters$Celulas, decreasing = TRUE), ]
#
# cat("\nIdentificación de Clusters Leiden (Basado en labels reales):\n")
# print(df_clusters, row.names = FALSE)


# #Visualización con Heatmap
# pheatmap(score_matrix_PCA,
#          cluster_rows = TRUE,
#          cluster_cols = FALSE,
#          display_numbers = TRUE,
#          number_format = "%.2f",
#          main = "Especificidad: Z-score Promedio por Firma",
#          color = colorRampPalette(c("navy", "white", "firebrick"))(100))
