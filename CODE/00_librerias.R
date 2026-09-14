suppressPackageStartupMessages({

  # --- 1. ESTRUCTURAS DE DATOS Y MANEJO DE MATRICES ---
  library(SingleCellExperiment) # Contenedor principal para datos de RNA-seq single-cell en Bioconductor
  library(Matrix)               # Operaciones con matrices dispersas (dgCMatrix)
  library(MatrixGenerics)       # Funciones genéricas para resúmenes de matrices (rowVars, colMeans)
  library(matrixStats)          # Cálculo eficiente de estadísticas por filas/columnas en matrices densas

  # --- 2. ACCESO Y OBTENCIÓN DE DATASETS DE PRUEBA ---
  library(scRNAseq)             # Colección de datasets públicos anotados de single-cell desde Bioconductor
  library(TENxPBMCData)         # Carga directa de datasets de PBMC procesados por 10x Genomics
  library(org.Hs.eg.db)         # Mapeo y anotación de genes humanos (Ensembl, Entrez, Symbol)
  library(splatter)             # Simulación sintética de datos RNA-seq de single-cell (útil para benchmarking)

  # --- 3. CONTROL DE CALIDAD Y PROCESAMIENTO INICIAL ---
  library(scuttle)              # Métricas de QC (conteo de mitocondriales, filtrado de células) y normalización
  library(scDblFinder)          # Detección y eliminación automática de dobletes (cell doublets)
  library(zinbwave)             # Modelado del exceso de ceros (zero-inflation) y reducción de dimensión
  library(scry)                 # Selección de genes altamente variables mediante devianzas de Poisson

  # --- 4. ANALISIS PIPELINE SINGLE-CELL (CORE) ---
  library(scran)                # Normalización (deconvolution pooling) y selección de genes variables (HVGs)
  library(scater)               # Control de calidad, visualización QC y diagnóstico de expresiones
  library(Seurat)               # Pipeline integral alternativo para procesamiento y análisis scRNA-seq
  library(sctransform)          # Normalización mediante regresión de Poisson corregida por varianza (Seurat)

  # --- 5. REDUCCIÓN DE DIMENSIONALIDAD ---
  library(uwot)                 # Implementación rápida en C++ del algoritmo UMAP
  library(umap)                 # Implementación alternativa del algoritmo UMAP
  library(Rtsne)                # Implementación del algoritmo t-SNE (Barnes-Hut)
  library(dimRed)               # Framework unificado para métodos de reducción de dimensionalidad

  # --- 6. CLUSTERING Y REDES ---
  library(bluster)              # Algoritmos de clustering para scRNA-seq (KNN graphs, Louvain, Walktrap)
  library(cluster)              # Algoritmos clásicos de clustering (PAM, Silhouette, Hierarchical)
  library(igraph)               # Construcción y análisis de grafos y redes complejas

  # --- 7. MÉTRICAS DE EVALUACIÓN Y BENCHMARKING ---
  library(mclust)               # Provee la función adjustedRandIndex() (ARI) para comparar clusters
  library(aricode)              # Provee la función NMI() (Normalized Mutual Information) y ARI rápida
  library(coRanking)            # Evaluación de calidad de reducciones de dimensión (preservación de distancia)
  library(entropy)              # Estimación de entropía e información mutua en distribuciones de datos
  library(caret)                # Métricas de validación cruzada y entrenamiento de modelos de Machine Learning
  library(MASS)                 # Funciones estadísticas avanzadas (LDA, regresiones robustas)

  # --- 8. MANIPULACIÓN DE DATOS (TIDYVERSE) ---
  library(dplyr)                # Manipulación eficiente de data frames (filter, mutate, summarize)
  library(tidyr)                # Remodelación y limpieza de tablas (pivot_longer, pivot_wider)
  library(reshape2)             # Transformación de datos entre formatos anchos y largos (melt, dcast)

  # --- 9. VISUALIZACIÓN Y FIGURAS ---
  library(ggplot2)              # Sistema principal para la generación de gráficos declarativos
  library(patchwork)            # Combinación modular de múltiples objetos ggplot2 en una sola figura
  library(pheatmap)             # Generación de mapas de calor estáticos personalizables
  library(RColorBrewer)         # Paletas de colores optimizadas para datos discretos y continuos
  library(viridis)              # Paletas de colores continuas accesibles (daltónicos/perceptuales)
  library(scales)               # Formateo y transformación de ejes/escalas en gráficos (log, porcentaje)
  library(ggforce)              # Extensiones avanzadas para ggplot2 (elipses de concentración, zooms)
  library(grid)                 # Funciones gráficas de bajo nivel para manipulación de layouts

  # --- 10. COMPUTACIÓN EN PARALELO ---
  library(BiocParallel)         # Soporte para evaluación en paralelo en Bioconductor (multi-core)
})
