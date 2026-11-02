##Modelo de mezclas
NB_mixture <- function(y,
                       size_factors,
                       K_max,
                       alpha = 1.0,
                       max_iter = 150,
                       tol = 1e-6,
                       min_phi = 1e-3,
                       max_phi = 1e3,
                       block_size = 2000,
                       prune_threshold = 1e-4,
                       lambda_mu = 1e-2,
                       update_phi_freq = 5,
                       verbose = TRUE
) {

  n <- nrow(y)
  G <- ncol(y)
  K <- K_max

  if (verbose) cat(sprintf("NB-MAPEM (Strict Math): n=%d, G=%d, K_max=%d\n", n, G, K))

  #Normalización de profundidad
  size_factors <- size_factors / mean(size_factors)
  size_factors[size_factors <= 0] <- 1

  #Inicialización Bayes Empírico de Dispersión
  gene_means <- colMeans(y)
  gene_vars <- matrixStats::colVars(y)

  phi_g <- numeric(G)
  mask <- (gene_vars > gene_means + 1e-5) & (gene_means > 1e-5)
  phi_g[mask] <- (gene_means[mask]^2) / (gene_vars[mask] - gene_means[mask])
  phi_g[!mask] <- max_phi

  trend <- median(phi_g[is.finite(phi_g)])
  phi_g <- 0.7 * phi_g + 0.3 * trend
  phi_g <- pmax(pmin(phi_g, max_phi), min_phi)

  #Inicialización Poisson Warm-up
  init_r_nb <- function(y_mat, sf, K_init) {
    gene_v <- matrixStats::colVars(y_mat)
    # Usar los genes más variables para encontrar los clústers iniciales
    top_genes <- order(gene_v, decreasing = TRUE)[1:min(500, ncol(y_mat))]
    y_sub <- y_mat[, top_genes, drop=FALSE]

    seed_idx <- sample(1:nrow(y_mat), K_init)
    mu_init <- sweep(y_sub[seed_idx, , drop=FALSE], 1, sf[seed_idx], "/")

    logp <- matrix(0, nrow(y_mat), K_init)
    for(k in 1:K_init) {
      mu_ik <- outer(sf, mu_init[k, ])
      logp[, k] <- rowSums(dpois(y_sub, lambda = mu_ik, log = TRUE))
    }
    logp <- logp - matrixStats::rowMaxs(logp)
    r_init <- exp(logp)
    r_init / rowSums(r_init)
  }

  r <- init_r_nb(y, size_factors, K)

#Inicialización MAP de mu
  Nk_sf <- as.vector(crossprod(r, size_factors))
  mu_kg <- sweep(crossprod(r, y) + lambda_mu, 1, Nk_sf + lambda_mu, "/")
  mu_kg[is.na(mu_kg)] <- 1e-9
  mu_kg <- pmax(mu_kg, 1e-9)
  # Inicialización del Stick-breaking
  a_param <- rep(1, K)
  b_param <- rep(alpha, K)

  obj_trace <- numeric(max_iter)

  row_log_sum_exp <- function(mat) {
    m <- matrixStats::rowMaxs(mat)
    out <- m + log(rowSums(exp(mat - m)))
    out[is.infinite(m)] <- -1e308
    return(out)
  }

# CICLO EM
  for(iter in 1:max_iter) {

    #PASO E
    # E[log pi_k] basado en proceso de Dirichlet
    Elogv <- digamma(a_param) - digamma(a_param + b_param)
    Elog1mv <- digamma(b_param) - digamma(a_param + b_param)

    Elogpi <- numeric(K)
    cum <- 0
    for(k in 1:K) {
      Elogpi[k] <- Elogv[k] + cum
      if(k < K) cum <- cum + Elog1mv[k]
    }

    #Verosimilitud por Bloques
    logp_mat <- matrix(0, n, K)
    n_blocks <- ceiling(G / block_size)

    for(block_idx in 1:n_blocks) {
      g_idx <- ((block_idx-1)*block_size + 1):min(block_idx*block_size, G)
      y_sub <- y[, g_idx, drop=FALSE]
      mu_sub <- mu_kg[, g_idx, drop=FALSE]
      phi_vec <- rep(phi_g[g_idx], each=n)

      for(k in 1:K) {
        mu_ik <- outer(size_factors, mu_sub[k, ])
        ll_vec <- dnbinom(as.vector(y_sub), size = phi_vec, mu = as.vector(mu_ik), log = TRUE)
        logp_mat[, k] <- logp_mat[, k] + rowSums(matrix(ll_vec, nrow=n))
      }
    }

    #Responsabilidades (q(Z))
    log_r <- sweep(logp_mat, 2, Elogpi, "+")
    lse <- row_log_sum_exp(log_r)
    r <- exp(log_r - lse)

    #FUNCIÓN OBJETIVO VARIACIONAL (ELBO + MAP)
    ll_expected <- sum(lse)

    #Kernel matemático estricto del Prior Gamma(lambda_mu + 1, lambda_mu)
    prior_mu <- sum(lambda_mu * log(mu_kg) - lambda_mu * mu_kg)

    #Prior Beta del Stick-Breaking
    prior_pi <- sum((alpha - 1) * Elog1mv)

    obj <- ll_expected + prior_mu + prior_pi
    obj_trace[iter] <- obj

    #PODA ADAPTATIVA (mezcla sobreajustada)
    pi_est <- colSums(r) / n
    active <- (pi_est > prune_threshold) & !is.na(pi_est)

    if (sum(active) < K && sum(active) > 1) {
      r <- r[, active, drop=FALSE]
      r <- r / rowSums(r)
      mu_kg <- mu_kg[active, , drop=FALSE]
      a_param <- a_param[active]
      b_param <- b_param[active]
      K <- sum(active)
      if(verbose) cat(sprintf("  [Pruning] Componentes eliminados. K efectivo = %d\n", K))
    }

    if (verbose) {
      cat(sprintf("Iter %3d | Obj %.6e | K_eff %d\n", iter, obj, K))
    }

    # Convergencia
    if (iter > 2) {
      delta <- abs(obj_trace[iter] - obj_trace[iter-1]) / (abs(obj_trace[iter-1]) + 1e-10)
      if (delta < tol) {
        if(verbose) cat("Convergencia alcanzada.\n")
        break
      }
    }

    # PASO M (Parámetros)
    Nk <- colSums(r) + 1e-10
    a_param <- 1 + Nk
    b_param <- alpha + c(rev(cumsum(rev(Nk)))[-1], 0)

    # MAP exacto de mu_kg
    Nk_sf <- as.vector(crossprod(r, size_factors))
    mu_kg <- sweep(crossprod(r, y) + lambda_mu, 1, Nk_sf + lambda_mu, "/")
    mu_kg <- pmax(mu_kg, 1e-9)

    # Actualizamos periódicamente calculando momentos por bloques para ahorrar RAM
    if (iter %% update_phi_freq == 0) {
      phi_new <- numeric(G)

      for(block_idx in 1:n_blocks) {
        g_idx <- ((block_idx-1)*block_size + 1):min(block_idx*block_size, G)
        y_sub <- y[, g_idx, drop=FALSE]

        # mu_mix = sum_k r_ik * (s_i * mu_kg)
        mu_mix <- (r %*% mu_kg[, g_idx, drop=FALSE]) * size_factors

        mean_g <- colMeans(mu_mix)
        var_g <- colMeans((y_sub - mu_mix)^2)

        valid <- var_g > mean_g

        phi_new[g_idx[valid]] <- (mean_g[valid]^2) / (var_g[valid] - mean_g[valid])
        phi_new[g_idx[!valid]] <- max_phi
      }

      trend_iter <- median(phi_new[is.finite(phi_new)])
      phi_g <- 0.7 * phi_new + 0.3 * trend_iter
      phi_g <- pmax(phi_g, min_phi)
    }
  }

  #Normalización final de las esperanzas
  v_mean <- a_param / (a_param + b_param)
  Epi <- numeric(K)
  remaining <- 1
  for(k in 1:K){
    if(k < K) {
      Epi[k] <- v_mean[k] * remaining
      remaining <- remaining * (1 - v_mean[k])
    } else {
      Epi[k] <- remaining
    }
  }

  return(list(
    responsibilities = r,
    mu = mu_kg,
    phi = phi_g,
    pi = Epi,
    K_final = K,
    obj_trace = obj_trace[1:iter]
  ))
}


##Predictiva posterior.
#(tratar de usar el paquete matrixstats)
ppc <- function(Y_real,
                model, size_factors,
                B = 10) {

  # --- Funciones Internas ---
  col_vars <- function(x) {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::colVars(x)
    } else {
      apply(x, 2, var)
    }
  }

  outliers <- function(Y) {
    thr <- quantile(Y, 0.99, na.rm = TRUE)
    mean(Y > thr, na.rm = TRUE)
  }

  over_dispersion <- function(Y) {
    mean(col_vars(Y) / (colMeans(Y) + 1e-8), na.rm = TRUE)
  }

  variability <- function(Y) {
    var(rowSums(Y))
  }

  zero_frac <- function(Y) {
    mean(Y == 0)
  }

  simulate_from_model <- function(model, size_factors) {
    n_cells <- length(size_factors)
    K <- ncol(model$r)
    G <- ncol(model$mu)

    # Muestreo del estado latente (mezcla)
    k_vec <- apply(model$r, 1, function(p) sample.int(K, 1, prob = p))

    # Matriz de medias esperadas
    mu_mat <- model$mu[k_vec, ] * size_factors

    # Alineación del parámetro de dispersión
    phi_vec <- if(length(model$phi) > 1) rep(model$phi, each = n_cells) else model$phi

    # OPTIMIZACIÓN DE RAM: Enteros y retorno directo
    matrix(
      as.integer(rnbinom(n_cells * G, size = phi_vec, mu = as.vector(mu_mat))),
      nrow = n_cells, ncol = G
    )
  }

  # --- Estadísticos Reales ---
  G <- ncol(Y_real)
  real_mean <- colMeans(Y_real)
  real_var  <- col_vars(Y_real)
  real_zero <- colMeans(Y_real == 0)

  T_real <- c(
    disp    = over_dispersion(Y_real),
    tail    = outliers(Y_real),
    zero    = zero_frac(Y_real),
    lib_var = variability(Y_real)
  )

  # --- Simulación y Acumulación ---
  mean_acc <- numeric(G)
  var_acc  <- numeric(G)
  zero_acc <- numeric(G)

  T_rep <- matrix(0, nrow = B, ncol = length(T_real))
  colnames(T_rep) <- names(T_real)

  for(b in seq_len(B)) {
    cat(sprintf("\rReplicando %d / %d...", b, B))

    Yrep <- simulate_from_model(model, size_factors)

    mean_acc <- mean_acc + colMeans(Yrep)
    var_acc  <- var_acc  + col_vars(Yrep)
    zero_acc <- zero_acc + colMeans(Yrep == 0)

    T_rep[b, ] <- c(
      disp    = over_dispersion(Yrep),
      tail    = outliers(Yrep),
      zero    = zero_frac(Yrep),
      lib_var = variability(Yrep)
    )

    # SALVAVIDAS DE RAM: Liberar memoria explícitamente
    rm(Yrep)
    gc()
  }
  cat("\n¡Simulación completada!\n")

  # --- Resumen Final ---
  rep_mean <- mean_acc / B
  rep_var  <- var_acc / B
  rep_zero <- zero_acc / B

  p_values <- colMeans(sweep(T_rep, 2, T_real, `>=`))
  gene_error <- abs(real_mean - rep_mean)

  simplex <- if(exists("analyze_simplex_fit", mode = "function")) {
    analyze_simplex_fit(Y_real, model, size_factors)
  } else {
    NULL
  }

  res <- list(
    gene_stats = list(
      real = data.frame(mean = real_mean, var = real_var, zero = real_zero),
      rep  = data.frame(mean = rep_mean, var = rep_var, zero = rep_zero)
    ),
    discrepancy = list(
      T_real = T_real,
      T_rep_mean = colMeans(T_rep),
      p_values = p_values
    ),
    gene_diagnostic = list(
      error = gene_error,
      worst_genes = order(gene_error, decreasing = TRUE)[1:min(50, G)]
    ),
    simplex = simplex
  )

  class(res) <- "ppc_result"
  return(res)
}

#Extraccion de metricas
print_thesis_tables <- function(res,
                                dataset_name = "Dataset") {
  cat("\n======================================================\n")
  cat("RESULTADOS NUMÉRICOS PARA:", dataset_name, "\n")
  cat("======================================================\n\n")

  cat("--- Cuadro 1: Relación Media-Varianza ---\n")
  cat(sprintf("Observado: %.8f\n", mean(res$gene_stats$real$var, na.rm = TRUE)))
  cat(sprintf("Modelo:    %.8f\n\n", mean(res$gene_stats$rep$var, na.rm = TRUE)))

  cat("--- Cuadro 2: Proporción de Ceros ---\n")
  cat(sprintf("Observado: %.8f\n", mean(res$gene_stats$real$zero, na.rm = TRUE)))
  cat(sprintf("Modelo:    %.8f\n\n", mean(res$gene_stats$rep$zero, na.rm = TRUE)))

  cat("--- Cuadro 3: Expresión Media ---\n")
  cat(sprintf("Observado: %.8f\n", mean(res$gene_stats$real$mean, na.rm = TRUE)))
  cat(sprintf("Modelo:    %.8f\n", mean(res$gene_stats$rep$mean, na.rm = TRUE)))
  cat("======================================================\n")
}

#Graficos
plot_mean_variance <- function(res,
                               title = "Relación media-varianza",
                               subtitle = "Posterior Predictive Check") {

  pseudo_count <- 1e-4

  df_real <- data.frame(
    Mean = res$gene_stats$real$mean + pseudo_count,
    Variance = res$gene_stats$real$var + pseudo_count,
    Grupo = "Observados"
  )

  df_rep <- data.frame(
    Mean = res$gene_stats$rep$mean + pseudo_count,
    Variance = res$gene_stats$rep$var + pseudo_count,
    Grupo = "Modelo"
  )

  df <- rbind(df_real, df_rep)
  df$Grupo <- factor(df$Grupo, levels = c("Modelo", "Observados"))

  p <- ggplot(df, aes(x = Mean, y = Variance, color = Grupo)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(title = title, subtitle = subtitle, x = "Expresión media", y = "Varianza", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

plot_mean_zero <- function(res,
                           title = "Fracción de ceros vs Expresión media",
                           subtitle = "Posterior Predictive Check") {

  pseudo_count <- 1e-4

  df_real <- data.frame(
    Mean = res$gene_stats$real$mean + pseudo_count,
    Zero = res$gene_stats$real$zero,
    Grupo = "Observados"
  )

  df_rep <- data.frame(
    Mean = res$gene_stats$rep$mean + pseudo_count,
    Zero = res$gene_stats$rep$zero,
    Grupo = "Modelo"
  )

  df <- rbind(df_real, df_rep)
  df$Grupo <- factor(df$Grupo, levels = c("Modelo", "Observados"))

  p <- ggplot(df, aes(x = Mean, y = Zero, color = Grupo)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, x = "Expresión media", y = "Proporción de ceros", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

plot_mean_expression_density <- function(res,
                                         title = "Distribución de la expresión media",
                                         subtitle = "Posterior Predictive Check") {

  pseudo_count <- 1e-4
  mean_real <- res$gene_stats$real$mean + pseudo_count
  mean_rep  <- res$gene_stats$rep$mean + pseudo_count

  df <- data.frame(
    Expresion = c(mean_real, mean_rep),
    Grupo = factor(
      rep(c("Observados", "Modelo"), times = c(length(mean_real), length(mean_rep))),
      levels = c("Modelo", "Observados")
    )
  )

  p <- ggplot(df, aes(x = Expresion, fill = Grupo, color = Grupo)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(title = title, subtitle = subtitle, x = "Expresión media", y = "Densidad", fill = NULL, color = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

run_scalability <- function(Y, fit_fn, fractions = c(0.2, 0.4, 0.6, 0.8, 1.0), B = 5) {
  N <- nrow(Y)
  G <- ncol(Y)
  results <- list()
  idx <- 1

  cat("Iniciando análisis de escalabilidad...\n")
  cat(sprintf("Total células (N): %d | Total genes (G): %d\n\n", N, G))

  for (alpha in fractions) {
    n_cells <- floor(alpha * N)

    for (b in seq_len(B)) {
      cat(sprintf("Fracción: %.2f (n = %d) | Réplica %d de %d... ", alpha, n_cells, b, B))

      # 1. Submuestreo aleatorio de células
      sampled_cells <- sample(seq_len(N), n_cells, replace = FALSE)
      Y_sub <- Y[sampled_cells, , drop = FALSE]

      # 2. Medición del tiempo y ajuste del modelo
      t0 <- proc.time()

      # Aquí se ejecuta el modelo (debe devolver objective, K_active, K_effective)
      fit_res <- fit_fn(Y_sub)

      t1 <- proc.time()

      # 3. Cálculos
      exec_time <- as.numeric((t1 - t0)["elapsed"])
      norm_obj  <- fit_res$objective / (n_cells * G)

      # 4. Guardar resultados iteración
      results[[idx]] <- data.frame(
        alpha = alpha,
        n = n_cells,
        rep = b,
        time = exec_time,
        norm_obj = norm_obj,
        K_active = fit_res$K_active,
        K_effective = fit_res$K_effective
      )

      idx <- idx + 1
      cat("Hecho.\n")
    }
  }

  df_res <- do.call(rbind, results)

#Resumen

#a) Regresión de ley de potencias para el tiempo: log(T) = log(C) + beta*log(n)
  fit_time <- lm(log(time) ~ log(n), data = df_res)
  beta_est <- coef(fit_time)["log(n)"]
  r2_est   <- summary(fit_time)$r.squared

#b) Regresión lineal para la estabilidad de la función objetivo: L_norm = intercept + m*n
  fit_obj   <- lm(norm_obj ~ n, data = df_res)
  slope_obj <- coef(fit_obj)["n"]
  pval_obj  <- summary(fit_obj)$coefficients["n", "Pr(>|t|)"]

  # c) Promedios en la máxima fracción muestral (n = N)
  max_n <- max(df_res$n)
  df_max <- subset(df_res, n == max_n)
  k_fin_avg <- mean(df_max$K_active)
  k_eff_avg <- mean(df_max$K_effective)

  summary_tab <- data.frame(
    beta = beta_est,
    R2 = r2_est,
    pendiente_L = slope_obj,
    p_valor = pval_obj,
    K_final = k_fin_avg,
    K_efectivo = k_eff_avg
  )
  rownames(summary_tab) <- NULL

  res_out <- list(data = df_res, summary = summary_tab)
  class(res_out) <- "scalability_result"
  return(res_out)
}

#Tablas y graficas
print.scalability_result <- function(x, dataset_name = "Dataset", ...) {
  cat("\n======================================================================\n")
  cat(sprintf("Cuadro resumen: Propiedades de escalamiento (%s)\n", dataset_name))
  cat("======================================================================\n")
  cat(sprintf("Beta (exponente de tiempo) : %.3f\n", x$summary$beta))
  cat(sprintf("R^2 (regresión de tiempo)  : %.3f\n", x$summary$R2))
  cat(sprintf("Pendiente L_norm (por n)   : %e\n", x$summary$pendiente_L))
  cat(sprintf("p-valor (estabilidad L)    : %.4f\n", x$summary$p_valor))
  cat(sprintf("K_final (promedio en N)    : %.2f\n", x$summary$K_final))
  cat(sprintf("K_efectivo (promedio en N) : %.2f\n", x$summary$K_efectivo))
  cat("======================================================================\n")
}


plot.scalability_result <- function(x, ...) {
  df <- x$data
  tema_base <- theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11))

# 1. Arriba Izquierda: Tiempo vs n (Escala log-log)
  p1 <- ggplot(df, aes(x = n, y = time)) +
    geom_point(alpha = 0.5, color = "#2c3e50") +
    geom_smooth(method = "lm", color = "#e74c3c", se = TRUE, fill = "#e74c3c", alpha = 0.2) +
    scale_x_log10(labels = scales::math_format(10^.x), breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    scale_y_log10(labels = scales::math_format(10^.x), breaks = scales::trans_breaks("log10", function(x) 10^x)) +
    labs(title = "Escalamiento computacional", x = "Número de células (n)", y = "Tiempo de ejecución (s)") +
    tema_base

# 2. Arriba Derecha: L_norm vs n
  p2 <- ggplot(df, aes(x = n, y = norm_obj)) +
    geom_point(alpha = 0.5, color = "#2980b9") +
    geom_smooth(method = "loess", color = "#e74c3c", se = TRUE, alpha = 0.2) + # Loess muestra mejor la estabilización
    labs(title = bquote("Estabilidad de la función objetivo"~tilde(L)),
         x = "Número de células (n)", y = bquote(tilde(L) == L / (nG))) +
    tema_base

# 3. Abajo Izquierda: K_efectivo vs n
  p3 <- ggplot(df, aes(x = n, y = K_effective)) +
    geom_point(alpha = 0.5, color = "#27ae60") +
    geom_smooth(method = "loess", color = "#2c3e50", se = TRUE) +
    labs(title = "Número efectivo de componentes", x = "Número de células (n)", y = bquote(K[efectivo])) +
    tema_base

# 4. Abajo Derecha: K_final (activos) vs n
  p4 <- ggplot(df, aes(x = n, y = K_active)) +
    geom_point(alpha = 0.5, color = "#8e44ad") +
    geom_smooth(method = "loess", color = "#2c3e50", se = TRUE) +
    labs(title = "Número de componentes activos", x = "Número de células (n)", y = bquote(K[final])) +
    tema_base

# Combinar los gráficos usando patchwork
  plot_final <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title = "Propiedades empíricas de escalamiento del modelo",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )

  return(plot_final)
}
