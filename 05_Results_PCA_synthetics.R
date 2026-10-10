cluster_counts <- table(vector_clusters)

pi_leiden <- as.numeric(cluster_counts) /
             sum(cluster_counts)

entropy_pi <- -sum(pi_leiden * log(pi_leiden + 1e-12))

K_eff <- exp(entropy_pi)

cat("K efectivo (Leiden):", K_eff, "\n")

ARI_pca <- adjustedRandIndex(
  labels,
  vector_clusters
)

NMI_pca <- NMI(
  labels,
  vector_clusters
)

confusion_pca <- table(
  labels,
  vector_clusters
)

##Continuidad y Trustworthiness
# Espacio original (alta dimensión)
X_high <- t(pearson_mat[hvg_pearson, ])

# Espacio reducido (PCA)
X_low <- reduced_dims

#Funciones
trustworthiness_vera <- function(X_original, X_reducido, k = 20) {

  N <- nrow(X_original)

  D_orig <- as.matrix(dist(X_original))
  D_red  <- as.matrix(dist(X_reducido))

  diag(D_orig) <- Inf
  diag(D_red)  <- Inf

  penalty <- 0

  for(i in seq_len(N)) {

    vecinos_red <- order(D_red[i, ])[1:k]

    orden_orig <- order(D_orig[i, ])
    rango_orig <- integer(N)
    rango_orig[orden_orig] <- seq_len(N)

    intrusos <- vecinos_red[rango_orig[vecinos_red] > k]

    if(length(intrusos) > 0)
      penalty <- penalty + sum(rango_orig[intrusos] - k)
  }

  1 - (2 * penalty) /
    (N * k * (2 * N - 3 * k - 1))
}

continuity_vera <- function(X_original, X_reducido, k = 20) {

  N <- nrow(X_original)

  D_orig <- as.matrix(dist(X_original))
  D_red  <- as.matrix(dist(X_reducido))

  diag(D_orig) <- Inf
  diag(D_red)  <- Inf

  penalty <- 0

  for(i in seq_len(N)) {

    vecinos_orig <- order(D_orig[i, ])[1:k]

    orden_red <- order(D_red[i, ])
    rango_red <- integer(N)
    rango_red[orden_red] <- seq_len(N)

    omitidos <- vecinos_orig[rango_red[vecinos_orig] > k]

    if(length(omitidos) > 0)
      penalty <- penalty + sum(rango_red[omitidos] - k)
  }

  1 - (2 * penalty) /
    (N * k * (2 * N - 3 * k - 1))
}

#trust <- trustworthiness_vera(X_high, X_low, k = 20)
cont  <- continuity_vera(X_high, X_low, k = 20)

#cat("Trustworthiness:", round(trust, 4), "\n")
cat("Continuity:", round(cont, 4), "\n")



print(list(
  K_efectivo       = K_eff,
  ARI              = ARI_pca,
  NMI              = NMI_pca,
  Trustworthiness  =trust#,
  Continuity       =cont
))
