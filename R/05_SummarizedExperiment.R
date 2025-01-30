# Para correr la el paquete de SummarizedExperiment.
# Es útil para almacenar datos de experimentos biológicos.
library("SummarizedExperiment")

# Número de genes - filas
nrows <- 200

# Número de muestras - columnas
ncols <- 6

# Comando para mantener estáticos los números generados de manera aleatoria.
set.seed(20210223)
# Matriz con los números generados de manera aleatoria.
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)

# Se generan los datos de todos los genes.
rowRanges <- GRanges(
    rep(c("chr1", "chr2"), c(50, 150)),
    IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
    strand = sample(c("+", "-"), 200, TRUE),
    feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

# Se generan los datos de las muestras.
colData <- DataFrame(
    Treatment = rep(c("ChIP", "Input"), 3),
    row.names = LETTERS[1:6]
)

# Se lleva a cabo la creación del objeto SummarizedExperiment.
# Donde se junta todo lo que se llevó a cabo con anterioridad.
rse <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    rowRanges = rowRanges,
    colData = colData
)

# Se comprueba el objeto.
rse
dim(rse)
dimnames(rse)
head(assay(rse))
rowRanges(rse)
rowData(rse)
colData(rse)

# Ejercicio 1:
rse[1:2, ]
# RESPUESTA: Imprime las primeras dos filas de todas las columnas del objeto rse.

# Ejercicio 2:
rse[, c("A", "D", "F")]
# RESPUESTA: Imprime todas las filas de las columnas correspondientes a las letras A, D y F.

# La función identical verifica si ambas salidas son iguales.
identical(rse[, c("A", "D", "F")], rse[, c(1, 4, 6)])
# Con la respuesta TRUE se verifica que son iguales.

# En caso de que no fueran iguales el código se detiene.
stopifnot(rse[, c("A", "D", "F")] == rse[, c(1, 4, 6)])

# Se juntan ambas funciones.
stopifnot(identical(rse[, c("A", "D", "F")], rse[, c(1, 4, 6)])

# Descargar la librería iSEE.
library("iSEE")
iSEE::iSEE(rse)
# Se llevan a cabo las gráficas con iSEE de los datos generados de manera azarosa.

# Data set - descarga de datos
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer

# Myestra las imágenes y gráficas generadas con iSEE de los datos descargados.
iSEE::iSEE(sce_layer)

# Código para ejecutar el PCA:
## The following list of commands will generate the plots created in iSEE
## Copy them into a script or an R session containing your SingleCellExperiment.
## All commands below refer to your SingleCellExperiment object as `se`.

se <- sce_layer
colormap <- ExperimentColorMap()
se <- iSEE::cleanDataset(se)
colormap <- synchronizeAssays(colormap, se)
all_contents <- list()

################################################################################
# Defining brushes
################################################################################

all_active <- list()
all_active[['ReducedDimensionPlot1']] <- list()
all_active[['FeatureAssayPlot1']] <- list()
all_active[['ColumnDataPlot1']] <- list()
all_active[['RowDataPlot1']] <- list()
all_active[['SampleAssayPlot1']] <- list()

################################################################################
## Reduced dimension plot 1
################################################################################


red.dim <- reducedDim(se, "PCA");
plot.data <- data.frame(X=red.dim[, 1], Y=red.dim[, 2], row.names=colnames(se));

plot.data$ColorBy <- colData(se)[, "layer_guess_reordered_short"];

# Avoid visual biases from default ordering by shuffling the points
set.seed(76);
plot.data <- plot.data[sample(nrow(plot.data)),,drop=FALSE];

dot.plot <- ggplot() +
    geom_point(aes(x=X, y=Y, color=ColorBy), alpha=1, plot.data, size=4) +
    labs(x="Dimension 1", y="Dimension 2", color="layer_guess_reordered_short", title="PCA") +
    coord_cartesian(xlim=range(plot.data$X, na.rm=TRUE),
        ylim=range(plot.data$Y, na.rm=TRUE), expand=TRUE) +
    scale_color_manual(values=colDataColorMap(colormap, "layer_guess_reordered_short", discrete=TRUE)(7), na.value='grey50', drop=FALSE) +
    scale_fill_manual(values=colDataColorMap(colormap, "layer_guess_reordered_short", discrete=TRUE)(7), na.value='grey50', drop=FALSE) +
    guides(colour = guide_legend(override.aes = list(size=1)), fill = guide_legend(override.aes = list(size=1))) +
    theme_bw() +
    theme(legend.position='bottom', legend.box='vertical', legend.text=element_text(size=9), legend.title=element_text(size=11),
            axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12))

# Saving data for transmission
all_contents[['ReducedDimensionPlot1']] <- plot.data

# Código para ejecutar el heatmap:


################################################################################
## Complex heatmap 1
################################################################################


.chosen.rows <- c("ENSG00000168314", "ENSG00000183036", "ENSG00000197971");
.chosen.columns <- colnames(se);
plot.data <- assay(se, "logcounts")[.chosen.rows, .chosen.columns, drop=FALSE]
plot.data <- as.matrix(plot.data);

plot.data <- plot.data - rowMeans(plot.data)
plot.data <- plot.data / apply(plot.data, 1, sd)

.assay_colors <- c("purple", "black", "yellow")
.assay_colors <- circlize::colorRamp2(breaks = c(-2.17978728210008, 0, 2.66762445840324), colors = .assay_colors)

# Keep all samples to compute the full range of continuous annotations
.column_data <- colData(se)[, "layer_guess_reordered_short", drop=FALSE]
.column_data[["Selected points"]] <- iSEE::multiSelectionToFactor(list(), colnames(se))

.column_col <- list()

.color_values <- .column_data[["layer_guess_reordered_short"]]
.color_values <- setdiff(unique(.color_values), NA)
.col_colors <- colDataColorMap(colormap, "layer_guess_reordered_short", discrete=TRUE)(length(.color_values))
if (is.null(names(.col_colors))) names(.col_colors) <- levels(factor(.color_values))
.column_col[["layer_guess_reordered_short"]] <- .col_colors

.column_col[["Selected points"]] <- iSEE::columnSelectionColorMap(colormap, levels(.column_data[["Selected points"]]))

.column_data <- .column_data[colnames(plot.data), , drop=FALSE]
.column_data <- as.data.frame(.column_data, optional=TRUE)
.column_annot_order <- order(.column_data[["Selected points"]], .column_data[["layer_guess_reordered_short"]])
.column_data <- .column_data[.column_annot_order, , drop=FALSE]
plot.data <- plot.data[, .column_annot_order, drop=FALSE]
.column_annot <- ComplexHeatmap::columnAnnotation(df=.column_data, col=.column_col, annotation_legend_param=list(direction="horizontal", nrow=10))

hm <- ComplexHeatmap::Heatmap(matrix=plot.data, col=.assay_colors,
    top_annotation=.column_annot, cluster_rows=TRUE,
    clustering_distance_rows="spearman", clustering_method_rows="ward.D2",
    cluster_columns=FALSE, name="logcounts (centered, scaled)",
    show_row_names=TRUE, show_column_names=TRUE,
    row_names_gp=grid::gpar(fontsize=8),
    column_names_gp=grid::gpar(fontsize=5),
    heatmap_legend_param=list(direction="horizontal"))

ComplexHeatmap::draw(hm, heatmap_legend_side="bottom", annotation_legend_side="bottom")

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()

# Ambas imágenes se encuentran como imágenes en el Git.

# Explora en con un heatmap la expresión de los genes MOBP, MBP y PCP4. Si hacemos un clustering (agrupamos los genes), ¿cúales genes se parecen más?
# RESPUESTA: MOBP y MBP se parecen más.

# ¿En qué capas se expresan más los genes MOBP y MBP?
# RESPUESTA: En la capa WM.
