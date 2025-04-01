#Cargar librerías
library(readxl)
library(tidyverse)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(iSEE)

#Cargar datos
raw_data <- read_excel("C:/Users/Judit/Downloads/GastricCancer_NMR.xlsx", sheet = "Data")
metadata <- read_excel("C:/Users/Judit/Downloads/GastricCancer_NMR.xlsx", sheet = "Peak")

#Limpiar datos, eliminamos los metabolitos que tengan más del 20% missing

# Primero identificamos los metabolitos con >10% de datos faltantes en metadata
metabolitos_a_conservar <- metadata$Name[metadata$Perc_missing <= 10]

# Filtramos ambos dataframes para mantener solo los metabolitos con ≤10% de datos faltantes
metadata<- metadata[metadata$Perc_missing  <= 10, ]
# Verificar qué columnas de raw_data corresponden a metabolitos
columnas_metabolitos <- grep("^M", colnames(raw_data), value = TRUE)
metabolitos_a_conservar_en_raw <- columnas_metabolitos[columnas_metabolitos %in% metabolitos_a_conservar]

# Crear un dataframe filtrado con las columnas de interés
raw_data_filtrado <- raw_data %>%
  dplyr::select(SampleID, SampleType, Class, all_of(metabolitos_a_conservar_en_raw))
raw_data <- raw_data_filtrado






#Seleccionar solo las columnas relevantes (SampleID y mediciones)
  samples <- raw_data %>%
  dplyr::select(SampleID, starts_with("M"))

#Transponer los datos para que SampleID sean columnas
df_transposed <- samples %>%
pivot_longer(cols = -SampleID, names_to = "Measurement", values_to = "Value") %>%
  pivot_wider(names_from = SampleID, values_from = Value)

#Convertir en matriz de conteos
counts_matrix <- df_transposed %>%
  column_to_rownames(var = "Measurement") %>%
  as.matrix()

#Crear `rowData` con metadatos de mediciones
rowData_df <- metadata %>%
  dplyr::select(-Idx) %>%
  column_to_rownames(var = "Name")  

#Crear `colData` con metadatos de muestras
colData_df <- raw_data %>%
  dplyr::select(SampleID, SampleType, Class) %>%
  column_to_rownames(var = "SampleID")

#Crear `SummarizedExperiment`
se <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  colData = colData_df,
  rowData = rowData_df
)

# Verificar que el objeto se creó correctamente
se
#SumarizedExperimente observaciones
assays(se)$counts
rowRanges(se)
colData(se)

###########################################
#Visualizar Summarized Experiment con iSEE# (se abre una ventana en el navegador)
###########################################

iSEE(se, initial = list(
  FeatureAssayPlot(
    PanelId = 1L,  # Convertimos el identificador en entero
    Assay = "counts",  
    XAxis = "Column data",  
    XAxisColumnData = "Class"
  ),
  RowDataTable(PanelId = 2L)  # También corregimos aquí
))

#Para visualizarlo en el informe se empleó el siguiente código
se <- iSEE::cleanDataset(se)
colormap <- synchronizeAssays(ExperimentColorMap(), se)

set.seed(100)
plot.data <- data.frame(Y = assay(se, "counts")["M1", ], X = factor(colData(se)[, "Class"])) |> 
  subset(!is.na(Y)) |> 
  transform(GroupBy = X, jitteredX = iSEE::jitterViolinPoints(X, Y, width = 0.4, method = 'quasirandom')) 

set.seed(124)
plot.data <- plot.data[sample(nrow(plot.data)), , drop = FALSE]

dot.plot <- ggplot(plot.data, aes(x = X, y = Y)) +
  geom_violin(aes(group = GroupBy), alpha = 0.2, scale = 'width', width = 0.8) +
  geom_point(aes(x = jitteredX), alpha = 1, color = '#000000', size = 1) +
  labs(x = "Class", y = "M1 (counts)", title = "M1 vs Class") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm = TRUE), expand = TRUE) +
  theme_bw() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), title = element_text(size = 12))

dot.plot

#################################################################
###Analisis de componentes pricipales con la librería MixOmics###
#################################################################


library(mixOmics)
X <- raw_data[, -c(1:4)]
Y <- raw_data$Class

# Realizar PCA usando mixOmics
pca_modelo <- pca(X, ncomp = 5, scale = TRUE)

# Gráfico de scores
plotIndiv(pca_modelo, group = Y,ind.names = FALSE, ellipse = TRUE, legend = TRUE,
          title = "PCA Scores: PC1 vs PC2")