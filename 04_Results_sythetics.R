#La diferencia con el bloque de BaronPancreas es minima, pero prefiero dejarlo asi para poder mantener todo separado, aparte que si tenia controlados mis scripts de forma local.

##Estructura del modelo
r_full  <- mixture$r
mu_full <- mixture$mu
phi      <- mixture$phi
n_cells <- nrow(r_full)
K        <- ncol(r_full)
G        <- ncol(matrix_counts)
hard_clusters <- max.col(r_full)
pi_est <- colSums(r_full)/n_cells
entropy_pi <- -sum(pi_est*log(pi_est+1e-12))
K_eff <- exp(entropy_pi)

##Incertidumbre
entropy_r <- -rowSums(r_full*log(r_full+1e-12))
max_prob <- apply(r_full,1,max)
r_sorted <- t(apply(r_full,1,sort,decreasing=TRUE))
gap <- r_sorted[,1]-r_sorted[,2]

##Anotacion biologica
#Esto es para la simulacion no nula
# labels <- factor(cells$label)
# Y_label <- model.matrix(~ labels - 1)

#Esto es para las simulaciones nulas
labels <- factor(cells$true_group)
Y_label <- diag(nlevels(labels))[as.integer(labels), , drop = FALSE]
colnames(Y_label) <- levels(labels)

#Composición de cada cluster
cluster_type_prob <- t(r_full)%*%Y_label
cluster_type_prob <- cluster_type_prob / rowSums(cluster_type_prob)
cluster_type_prob[is.na(cluster_type_prob)] <- 0

#Probabilidad de tipo celular por célula
cell_type_prob <- r_full %*% cluster_type_prob

#Entropía biológica
entropy_type <- -rowSums(cell_type_prob*log(cell_type_prob+1e-12))
entropy_type_norm <- entropy_type / log(ncol(cell_type_prob))

#Abundancias
cell_type_abundance <- colSums(cell_type_prob)
tabla_tipos <- data.frame(
  Tipo = colnames(cell_type_prob),
  Masa = cell_type_abundance,
  Porcentaje =
    round(cell_type_abundance/n_cells*100,2)
)

##Log-verosimilitud
loglik_cell_nb <- function(y,mu,phi,pi,sf){
  n <- nrow(y)
  K <- length(pi)
  logp <- matrix(0,n,K)
  for(k in 1:K){
    mu_k <- sweep(matrix(mu[k,],n,ncol(y),byrow=TRUE),1,sf,"*")
    logp[,k] <- rowSums(dnbinom(y,size=phi,mu=mu_k,log=TRUE)) +
      log(pi[k]+1e-12)
  }
  matrixStats::rowLogSumExps(logp)
}

ll_cell <- loglik_cell_nb(matrix_counts,mu_full,phi,pi_est,sf)

ll_total <- sum(ll_cell)
ll_per_cell <- ll_total/n_cells
ll_per_gene <- ll_total/(n_cells*G)

##Rareza
thr_ll <- quantile(ll_cell,0.01)
rare_ll <- which(ll_cell < thr_ll)
thr_entropy <- quantile(entropy_r,0.99)
rare_entropy <- which(entropy_r > thr_entropy)
rare_cells <- intersect(rare_ll,rare_entropy)
rarity_score <- scale(-ll_cell) + scale(entropy_r)
top_rare <- order(rarity_score,decreasing=TRUE)[1:50]

##Ajuste del modelo
ARI_nb <- adjustedRandIndex(
  labels,
  hard_clusters
)

NMI_nb <- NMI(
  labels,
  hard_clusters
)

##Geometria
#Trustworthiness
trustworthiness <- function(Y_original, Y_reducido, k = 20) {
  N <- nrow(Y_original)

  D_orig <- as.matrix(dist(Y_original))
  D_red  <- as.matrix(dist(Y_reducido))

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
  1 - (2 * penalty) / (N * k * (2 * N - 3 * k - 1))
}
#Continuty
continuity <- function(Y_original, Y_reducido, k = 20) {

  N <- nrow(Y_original)

  D_orig <- as.matrix(dist(Y_original))
  D_red  <- as.matrix(dist(Y_reducido))

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

sf_norm <- sf / mean(sf)
sf_norm[sf_norm <= 0] <- 1
# Alta Dimensión Pura: Solo conteos divididos por el factor de tamaño
Y_high <- sweep(matrix_counts, 1, sf_norm, "/")
#Baja Dimensión: Raíz cuadrada de 'r' extraída del modelo
Y_low <- sqrt(r_full)

trust <- trustworthiness(Y_high, Y_low, k = 20)
cont  <- continuity(Y_high, Y_low, k = 20)

print(list(
GEOMETRIA = list(
Trustworthiness = trust,
Continuity      = cont
),

LABELS = list(
ARI = ARI_nb,
NMI = NMI_nb
),

COMPLEJIDAD = list(
K_efectivo     = K_eff,
Entropia_pi    = entropy_pi,
Proporciones_pi = pi_est
),

INCERTIDUMBRE = list(
Entropia_media = mean(entropy_r),
Entropia_mediana = median(entropy_r),
Confianza_media = mean(max_prob),
Gap_medio = mean(gap)
),


AJUSTE = list(
LogLik_total = ll_total,
LogLik_por_celula = ll_per_cell,
LogLik_por_gen = ll_per_gene
),

BIOLOGIA = list(
Numero_tipos = nrow(tabla_tipos),
Entropia_tipo_media = mean(entropy_type_norm),
Abundancias = tabla_tipos
),

RAREZA = list(
N_rare_ll = length(rare_ll),
N_rare_entropy = length(rare_entropy),
N_outliers = length(rare_cells)
)
))
