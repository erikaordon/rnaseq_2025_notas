## ----'start', message=FALSE-----------------------------------
## Load recount3 R package
library("recount3")


## ----'quick_example'------------------------------------------
## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()

# Para saber el tipo de objeto
class(human_projects)

# que tan grande es
dim(human_projects)

# Inicio y final del objeto
head(human_projects)
tail(human_projects)

head(human_projects[order(human_projects$n_samples, decreasing = TRUE), ])
head(human_projects[order(human_projects$n_samples, decreasing = FALSE), ])

summary(human_projects$n_samples)
table(human_projects$n_samples)

## Encuentra tu proyecto de interés. Aquí usaremos
## SRP009615 de ejemplo
proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
)

class(proj_info)
dim(proj_info)

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)
## Explora el objeto RSE
rse_gene_SRP009615

# 175 variables de información de las muestras.

# Saber la cantidad de genes y muuestras respectivamente
nrow(rse_gene_SRP009615)
ncol(rse_gene_SRP009615)

# Explorar metadatos
metadata(rse_gene_SRP009615)

rowRanges(rse_gene_SRP009615)
rowData(rse_gene_SRP009615)

sort(table(seqnames(rowRanges(rse_gene_SRP009615)))), decreasing = TRUE)

# Forma interactiva del subset.
proj_info_interactive <- interactiveDisplayBase::display(human_projects)

# Nombres de assays.
assayNames(rse_gene_SRP009615)

## ----"interactive_display", eval = FALSE----------------------
# ## Explora los proyectos disponibles de forma interactiva
# proj_info_interactive <- interactiveDisplayBase::display(human_projects)
# ## Selecciona un solo renglón en la tabla y da click en "send".
#
# ## Aquí verificamos que solo seleccionaste un solo renglón.
# stopifnot(nrow(proj_info_interactive) == 1)
# ## Crea el objeto RSE
# rse_gene_interactive <- create_rse(proj_info_interactive)


## ----"tranform_counts"----------------------------------------
## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"--------------------------------------
## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
# Expander la información.

rse_gene_SRP009615$sra.sample_attributes

rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)
colData(rse_gene_SRP009615)[
    ,
    grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

# Ejercicio:
iSEE::iSEE(rse_gene_SRP009615)
