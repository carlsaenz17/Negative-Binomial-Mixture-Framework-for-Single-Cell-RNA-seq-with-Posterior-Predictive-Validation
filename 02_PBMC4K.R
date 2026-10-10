##Genes mitocondriales para el control de calidad.
genes_mito_humanos <- c(
  "ENSG00000198888","ENSG00000198763","ENSG00000198804","ENSG00000198712",
  "ENSG00000228253","ENSG00000198899","ENSG00000198938","ENSG00000198840",
  "ENSG00000212907","ENSG00000198886","ENSG00000198786","ENSG00000198695",
  "ENSG00000198727"
)

##Lista de firmas celulares para agrupamiento despues del modelo.
cell_group_ensembl <- list(

  ## 1. CD4+ T (naive/helper)
  "CD4+ T Cells" = c(
    "ENSG00000010610", # CD4
    "ENSG00000168685", # IL7R
    "ENSG00000138795"  # LEF1
  ),

  ## 2. CD8+ T (citotóxicas puras)
  "CD8+ T Cells" = c(
    "ENSG00000153563", # CD8A
    "ENSG00000171858", # CD8B
    "ENSG00000100450"  # GZMK (memoria/early cytotoxic)
  ),

  ## 3. NK (más específico)
  "NK Cells" = c(
    "ENSG00000143315", # KLRF1
    "ENSG00000115523", # GNLY
    "ENSG00000180644"  # PRF1
  ),

  ## 4. NK CD56 bright (separadas funcionalmente)
  "NK Bright" = c(
    "ENSG00000149294", # NCAM1
    "ENSG00000162591", # KLRB1
    "ENSG00000100450"  # GZMK
  ),

  ## 5. B Cells
  "B Cells" = c(
    "ENSG00000012124", # CD19
    "ENSG00000156738", # MS4A1
    "ENSG00000105369"  # CD79A
  ),

  ## 6. Plasmablasts (más específicos)
  "Plasmablasts" = c(
    "ENSG00000196611", # PRDM1
    "ENSG00000109046", # MZB1
    "ENSG00000211899"  # IGHG1
  ),

  ## 7. Monocitos clásicos
  "Classical Monocytes" = c(
    "ENSG00000170458", # CD14
    "ENSG00000143546", # S100A8
    "ENSG00000163220"  # S100A9
  ),

  ## 8. Monocitos no clásicos
  "Non-classical Monocytes" = c(
    "ENSG00000196126", # FCGR3A
    "ENSG00000203747", # MS4A7
    "ENSG00000167601"  # LILRB1
  ),

  ## 9. Dendritic Cells (cDC más específicas)
  "Dendritic Cells" = c(
    "ENSG00000179639", # FCER1A
    "ENSG00000120217", # HLA-DRA
    "ENSG00000106565"  # ITGAX (CD11c)
  ),

  ## 10. pDCs (limpias)
  "pDCs" = c(
    "ENSG00000198171", # CLEC4C
    "ENSG00000171513", # LILRA4
    "ENSG00000140614"  # TCF4
  ),

  ## 11. Tregs (más robusto)
  "Tregs" = c(
    "ENSG00000049768", # FOXP3
    "ENSG00000134460", # IL2RA
    "ENSG00000122886"  # CTLA4
  ),

  ## 12. Exhausted CD8+
  "Exhausted CD8+ T" = c(
    "ENSG00000188389", # PDCD1
    "ENSG00000135077", # LAG3
    "ENSG00000181019"  # HAVCR2
  ),

  ## 13. MAIT (más específica)
  "MAIT Cells" = c(
    "ENSG00000162591", # KLRB1
    "ENSG00000115415", # SLC4A10
    "ENSG00000112115"  # IL18RAP
  ),

  ## 14. Platelets
  "Platelets" = c(
    "ENSG00000163737", # PPBP
    "ENSG00000164692"  # PF4
  )
)

cat("Cargando PBMC4k...\n")
cells <- TENxPBMCData("pbmc4k")

##Anotacion de genes
symbols <- mapIds(org.Hs.eg.db,keys = rownames(cells),keytype = "ENSEMBL",column = "SYMBOL")
rowData(cells)$symbol <- symbols

##Quitar genes que se consideran solo ruido estadistico
#Control de calidad sobre los genes mitocondriales, esto se toma de la lista "genes_mito_humanos"
is_mito <- rownames(cells) %in% genes_mito_humanos
#Control de calidad para genes ribosomales
is_ribo <- grepl("^RPL|^RPS", rowData(cells)$symbol)
#Control de calidad para la hemoglobina
is_hb <- grepl("^HB[ABDG]", rowData(cells)$symbol)


##Genes con poca expresion.
cells$sum <- colSums(counts(cells))
cells$detected <- colSums(counts(cells) > 0)

##Filtrado de genes.
#Porcentaje de genes mitocondriales
mito_counts <- colSums(counts(cells)[is_mito, , drop = FALSE])
cells$subsets_Mito_percent <- (mito_counts / cells$sum) * 100
#Porcentaje de genes ribosomales
ribo_counts <- colSums(counts(cells)[is_ribo, , drop = FALSE])
cells$subsets_Ribo_percent <- (ribo_counts / cells$sum) * 100
#Outliers
#lib_out  <- isOutlier(cells$sum, type = "lower", log = TRUE, nmads = 3)
lib_out  <- isOutlier(cells$sum, type = "both", log = FALSE, nmads = 3)
gene_out <- isOutlier(cells$detected, type = "lower", log = FALSE, nmads = 3)
mito_out <- isOutlier(cells$subsets_Mito_percent, type = "higher", nmads = 3)
ribo_out <- isOutlier(cells$subsets_Ribo_percent, type = "higher", nmads = 3)
#Combinar todos los filtros
discard <- lib_out | gene_out | mito_out | ribo_out
#Filtrado de celulas
cells <- cells[, !discard]
cat("Células retenidas tras QC completo:", ncol(cells), "\n")
#Eliminacion de genes tecnicos
remove_genes <- is_mito | is_ribo | is_hb
cells <- cells[!remove_genes, ]
cat("Genes retenidos:", nrow(cells), "\n")
##Deteccion de celulas dobles
cells <- scDblFinder(cells)
cells <- cells[, cells$scDblFinder.class == "singlet"]
# ##Filtrado de genes poco expresados despues de eliminar el ruido tecnico.
# cells <- cells[rowSums(counts(cells) > 0) >= 10, ]

##Size factors (muy importantes)
cells <- computeLibraryFactors(cells)
sf <- sizeFactors(cells)
# sf <- colSums(counts(cells))
# sf <- sf / mean(sf)

##Creacion de la matriz de conteos.
matrix_counts <- as.matrix(t(counts(cells)))  # cells × genes

##Entrenamiento del modelo
mixture <- NB_mixture_MAPEM(
  y = matrix_counts,
  size_factors = sf,
  K = 30,
  max_iter = 150,
  block_size = 2000
)
