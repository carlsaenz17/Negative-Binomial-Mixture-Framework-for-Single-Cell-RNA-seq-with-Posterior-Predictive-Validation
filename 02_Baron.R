cat("Cargando Baron Pancreas Data...\n")
cells <- BaronPancreasData()

##Identificación de genes por Patrones de Texto (Symbols).
#En Baron(Humano), los mitocondriales empiezan por "MT-"
is_mito <- grepl("^MT-", rownames(cells), ignore.case = TRUE)
#Ribosomales (RPL y RPS)
is_ribo <- grepl("^RP[LS]", rownames(cells), ignore.case = TRUE)
#Hemoglobina(menos común en páncreas, pero útil para detectar contaminación de sangre)
is_hb <- grepl("^HB[ABDG]", rownames(cells), ignore.case = TRUE)

##Cálculo de métricas de calidad con scater
cells <- addPerCellQC(cells, subsets = list(Mito = is_mito, Ribo = is_ribo, HB = is_hb))

## Detección de Outliers
lib_out   <- isOutlier(cells$sum, type = "lower", log = TRUE, nmads = 3)
gene_out  <- isOutlier(cells$detected, type = "lower", log = TRUE, nmads = 3)

#En páncreas, el porcentaje mitocondrial puede ser más alto que en sangre. Si filtramos de forma muy estricta, perdemos células sanas.
mito_out  <- isOutlier(cells$subsets_Mito_percent, type = "higher", nmads = 3)

##Filtrado de Genes de Bajo Valor Informativo
#Eliminamos genes que no se expresan en casi ninguna célula
keep_feature <- rowSums(counts(cells) > 0) > 5
cells <- cells[keep_feature, ]

##Filtrado de células
discard <- lib_out | gene_out | mito_out
cells <- cells[, !discard]

##Detección de celulas dobles (Crucial en tejidos sólidos como páncreas)
#scDblFinder funciona mejor si le damos una pista de los clusters iniciales
cells <- scDblFinder(cells)
cells <- cells[, cells$scDblFinder.class == "singlet"]

cat("Células tras eliminar dobletes:", ncol(cells), "\n")

##Limpieza de Genes Técnicos
#Eliminamos Mito y HB de la matriz final para que no dominen el análisis posterior. Nota: A veces se dejan los ribosomales a menos que causen mucho ruido.
remove_genes <- is_mito | is_hb
cells <- cells[!remove_genes[keep_feature], ]

##Size factors (muy importantes)
cells <- computeLibraryFactors(cells)
sf <- sizeFactors(cells)

##Creacion de la matriz de conteos.
matrix_counts <- as.matrix(t(counts(cells)))  # cells × genes

#Proporcion de grupos celulares
# Conteo por label
 table(cells$label)
#
# Porcentaje por label
100 * prop.table(table(cells$label))
