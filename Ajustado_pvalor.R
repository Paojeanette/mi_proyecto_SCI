#ANALISIS DE EXPRESION DIFERENCIAL - RESPUESTA INMUNE EN SCI 

#Instalo librerias necesarias
install.packages("pacman")
install.packages("BiocManager")

# Cargar/ instalar los paquetes
pacman::p_load(
  dplyr,
  ggplot2,
  pheatmap,
  ggpubr,
  msigdbr,
  GEOquery,
  PCAtools,
  limma,
  reshape2,
  clusterProfiler,
  org.Hs.eg.db,   #Humano (Homo sapiens)
  enrichplot,
  readxl,               # Para leer Excel (los datos descargados)
  data.table,           # Lectura de tablas grandes
  matrixStats,
  AnnotationDbi         # Herramientas de anotación de genes
)



#Cargo los datos desde GEO con el ID del estudio
gse_id <- "GSE226238"       

print(gse_id)

#Recupero los datos de RNAseq mediante el ID de GSE 
gse_query <- getGEO(gse_id) 

print(gse_query)

#Cuando uso getGEO(), el resultado es una LISTA que puede contener uno o varios objetos.
#Cada elemento de esa lista suele ser un objeto de tipo "ExpressionSet".

#[[getGEO(gse_id) devuelve una lista de objetos ExpressionSet (uno por plataforma/GPL).
#A veces es una lista de longitud 1 (solo una plataforma) y otras veces trae varios]]-
length(gse_query)   #En este caso trae solo 1 elemento(1 plataforma)

# Un ExpressionSet es una estructura de Bioconductor que guarda:
#  La matriz de expresión (genes x muestras)
#  La metadata de las muestras
#  La información de los genes/sondas

gse_data <- gse_query[[1]] #Obtengo el primer objeto de ExpressionSet.

# Dentro de un ExpressionSet puedo explorar:
# pData(objeto) -> metadata de las muestras
pData(gse_data)[1:5,1:10]

# exprs(objeto) -> matriz de expresión (genes x muestras)
exprs(gse_data)[1:5,1:10] #Error in exprs(gse_data)[1:5, 1:10] : subscript out of bounds -- Esto significa que la matriz de expresion esta vacia. Los datos se guardaron en material complementario

dim(exprs(gse_data))   # Muestra el tamaño de la matriz de expresión


# fData(objeto) -> información de las sondas/genes -> Error in `[.data.frame`(fData(gse_data), 1:5, 1:10)- No hay info
fData(gse_data)[1:5,1:10]

# El objeto gse_data no trae anotaciones de features/genes.
# Cuando pongo print (gse_query) veo: assayData: 0 features, 46 samples ... featureData: none
# Eso significa que exprs(gse_data) no tiene filas (0 features) y no se cargó info de sondas/genes.
# Por eso fData(gse_data) devuelve un dataframe vacío (sin columnas).
# pData(gse_data) sí funciona porque hay metadatos de muestras (46).
# exprs(gse_data) está vacío porque la series matrix no incluye datos de expresión.
# fData(gse_data) está vacío porque no hay anotaciones de sondas.
# Los datos estan  almacenados en material suplementario ---> Siguiente paso, descargarlo

# metadatos disponibles por muestra
# Si voy a gse_data y veo los slots puedo encontrar phenoData -> data
# Sino puedo hacerlo mediante codigo de R, usando los nombres de las colnames de pData
colnames(pData(gse_data))           # nombres de variables #characteristics_ch1.2 ->Trae los tiempos


# Esto es para saber si los datos guardados en ExpressionSet son microarray o RNA-seq
annotation(gse_data) # "GPL18573"  -> plataforma de secuenciación (RNA-seq)
unique(pData(gse_data)$library_strategy)  

unique(pData(gse_data)$instrument_model) 



# Esto es una descripcion de experimentalData
experimentData(gse_data)


featureNames(gse_data)      # character(0) porque hay 0 features
sampleNames(gse_data)[1:5]  # nombres de muestras [1] "GSM7068702" "GSM7068703" "GSM7068704" "GSM7068705" "GSM7068706"

#NOTA::Las funciones de Biobase no se pueden aplicar a un string --> gse_id es solo "GSE226238", no un ExpressionSet.
#Por eso tengo que aplicarlas a gse_data !!



#############################################################################################
#             DESCARGO MATERIAL SUPLEMENTARIO - AHI ESTA,EN ESTE CASO, LA MATRIZ DE CONTEOS
##############################################################################################


#DESCARGAR ARCHIVOS SUPLEMENTARIOS
getGEOSuppFiles(gse_id)

#Ver qué archivos se descargaron
list.files(gse_id)   # "GSE226238_Morrison_et_al_processed_data.xlsx" -> Me descarga la matriz en excel

# Definir la ruta al archivo (está dentro de la carpeta con el ID del estudio)
ruta_excel <- file.path("GSE226238", "GSE226238_Morrison_et_al_processed_data.xlsx")

# Leer la primera hoja del Excel
datos <- read_excel(ruta_excel)

# Ver las primeras filas 
head(datos)


##############################################################################
#             LEER EL EXCEL Y PREPARAR LA MATRIZ (genes x muestras)
##############################################################################

##Pasos para manipular los datos de excel y preparar la matriz de conteos para el analisis
# 1) Paso a data frame para manipular mas fácil.
df <- as.data.frame(datos)

class(datos)  # da "tbl_df" "tbl" "data.frame"
class(df)     # "data.frame"

# 2) Pongo los nombres de muestra como rownames
rownames(df) <- df[["Library Name"]]   # filas = muestras
df[["Library Name"]] <- NULL           # quitar la columna de nombres

#ESTA PARTE NO LA HAGO, PORQUE ELIMINA EL NUMERO DE TRANSCRIPTO
# 3) Dejo solo el símbolo del gen (quitar lo que está antes de ">" que es el identificador de transcripto RefSeq mRNA)
#    Ej: "NM_000016>ACADM" -> "ACADM"
#colnames(df) <- sub(".*>", "", colnames(df)). 

#AGREGO ESTA PARTE PARA CONSERVAR TRANSCRIPTO Y POR OTRO LADO SEPARAR NOMBRE DE GEN PARA DICCIONARIO
#3) Mantengo las etiquetas originales (NM_xxx>GENE) como nombres de columnas
#    Esto conserva el ID del transcripto y el símbolo del gen juntos
#    Ej: "NM_000016>ACADM" se mantiene igual
#    (NO se borra nada, solo confirmo que están bien como vienen)

# 4) Creo un diccionario transcript_id -> gene_symbol
dict <- data.frame(
  transcript_id = colnames(df),                      # NM_xxx>GENE
  gene_symbol   = sub(".*>", "", colnames(df)),      # Solo el nombre del gen
  stringsAsFactors = FALSE
)
write.csv(dict, "diccionario_transcript_gene.csv", row.names = FALSE)

# 5) Me aseguro que todos los valores sean numéricos
df[] <- lapply(df, as.numeric)

# 6) Convertir a MATRIZ con transcriptos en filas y muestras en columnas
matriz <- t(as.matrix(df))

# 4) Me aseguro que todo sea numérico
df[] <- lapply(df, as.numeric)

# 5) Convertir a MATRIZ con ::: GENES en FILAS y MUESTRAS en COLUMNAS
#.   (Esto es xq algunos paquetes esperan este formato)
matriz <- t(as.matrix(df))

# Chequeo el contenido
dim(matriz)            # genes x muestras
matriz[1:5, 1:5]       # ver la matriz
head(matriz)


rownames(matriz)[1:5]       # primeros genes (símbolos)
colnames(matriz)[1:5]       # primeras muestras [1] "ID1v0"  "ID1v3"  "ID1V6"  "ID1v12" "ID2v0" 

#Como veo que hay mezcla en los nombres con minuscula y mayuscula, pongo todo en minuscula
colnames(matriz) <- tolower(colnames(matriz))

#Verifico que se modifico correctamente
colnames(matriz)[1:5] 



##############################################################################
#                  Extracción los datos clinicos que son de mi interés
##############################################################################

###Volviendo a pData <-  Para ver información clinica puedo:
#Extraer la información clinica (metadatos por muestras) -> Renombro para que sea mas sencillo, es igual que arriba pero resumo pasos)
clin_data <-  pData(gse_data)
colnames(clin_data) # Es igual que :: colnames(pData(gse_data)) 

head(clin_data$title)  #Etiquetas de muestras

head(clin_data$characteristics_ch1.1) #Grupos (SCI / Control)

head(clin_data$characteristics_ch1.2) #Tiempos post lesion


library(dplyr)


# Uso dplyr para hacer: select + rename + mutate 
sample_info <- clin_data %>%
  # 1) selecciono columnas de interés
  dplyr::select(title, `characteristics_ch1.1`, `characteristics_ch1.2`) %>%
  
  # 2) renombro (solo nombres de las columnas seleccionadas, no valores) Ej: RENOMBRO a 'title' como 'sample_label', idem con los otros
  dplyr::rename(
    sample_label = title,                        # "ID1v0, SCI, …". 
    group_raw    = `characteristics_ch1.1`,      # "group: SCI"/"group: CTL"
    treatment_raw= `characteristics_ch1.2`       # "treatment: Acute/3mpi/6mpi/12mpi"
  ) %>%
  
  
  # 3) limpio valores para sacar la parte de adelante y creo variables (Con MUTATE)
  # quito el prefijo "group: " y dejo solo "SCI"/"CTL"; lo paso a mayúsculas (por consistencia)
  # redefino group_raw como group
  # quito el prefijo "treatment: " y dejo "Acute"/"3mpi"/"6mpi"/"12mpi"
  # redefino tratment_raw como tratment.
  dplyr::mutate(
    sample_label = tolower(sub(",.*","", sample_label)),          # "ID1v0, ..." -> "id1v0"  (que tenga solo el id)
    group        = toupper(sub("^group:\\s*", "", group_raw)),    # "group: SCI" -> "SCI"    (que tenga solo SCI o CTL (no group: SCI))
    treatment    = sub("^treatment:\\s*", "", treatment_raw)      # "treatment: 3mpi" -> "3mpi"
  )





# Construir 'condition':
# creo un ÚNICO factor 'condition':
#   - para controles (group=="CTL") pongo "CTL"
#   - para pacientes SCI uso el valor de 'treatment' (Acute, 3mpi,6mpi,12mpi)
# ordeno los niveles del factor como quiero analizarlos y graficarlos
condition <- sample_info$treatment
condition[sample_info$group == "CTL"] <- "CTL"
condition <- factor(condition, levels = c("CTL","Acute","3mpi","6mpi","12mpi"))

# 4) Me quedo solo con las columnas finales en el orden deseado
sample_info <- sample_info %>%
  dplyr::select(sample_label, group, treatment) %>%
  dplyr::mutate(condition = condition)

# 5) alinear con la matriz
#ALINEO filas al ORDEN de las columnas de la matriz de expresión ---
# (filtra solo las muestras que están en la matriz y reordena para que
#  sample_label == colnames(matriz) en la misma posición)
sample_info <- sample_info %>%
  filter(sample_label %in% colnames(matriz)) %>%
  arrange(match(sample_label, colnames(matriz)))

# 6) chequeo : las etiquetas deben coincidir 1 a 1 con las columnas de 'matriz'
stopifnot(identical(sample_info$sample_label, colnames(matriz)))


##############################################################################
#                  CREO UN DICCIONARIO (Para usar despúes)
##############################################################################
library(AnnotationDbi)
library(org.Hs.eg.db)

matriz_gene <- matriz

## Usar los gene_symbol del diccionario
syms <- dict$gene_symbol

# solo símbolos únicos 
# y descartar NAs en los símbolos
syms_u <- syms[!duplicated(syms) & !is.na(syms)]


map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = syms_u,
  keytype = "SYMBOL",
  columns = c("SYMBOL","ENTREZID","ENSEMBL", "GENENAME")
)


# Eliminar filas duplicadas exactas (mismo contenido en todas las columnas)
map <- map[!duplicated(map), ]

# Garantizar una sola fila por SYMBOL (quedarse con la primera aparición)
map <- map[!is.na(map$SYMBOL) & !duplicated(map$SYMBOL), ]

# Reordenar el resultado para que siga el orden de primera aparición en syms_u
map <- map[match(syms_u, map$SYMBOL), ]

# Guardar 
write.csv(map, "gene_mapping_SYMBOL_ENTREZ_ENSEMBL.csv", row.names = FALSE)

# Veo contenido
head(map, 10)   # ver los primeros 10 genes
map[, c("SYMBOL","GENENAME","ENSEMBL","ENTREZID")] |> head()
colSums(is.na(map))   # cantidad de NAs por columna

# Ejemplo: revisar un símbolo concreto
map[map$SYMBOL == "PRKCD", ]


#Solo es un diccionario, no tiene conteos. Elimino duplicados de Symbol-


##############################################################################
#                        COMIENZO ANALISIS ESTADISTICO DESCRIPTIVO
##############################################################################

#Verifico si hay NA (en los conteos, expresion de genes, NO en los nombres de genes)
anyNA(matriz)         #Esta es la matriz de conteo (tiene duplicados)
sum(is.na(matriz))

#Uso apply para calcular mean,median, sd, min y max 
#Tengo que poner el dataframe + filas (1) columnas(2) + function(x) y escribo funciones a aplicar

#En las filas tengo genes
summary_stats_genes <- apply(matriz,1, function(x) {
  c(minimo = min(x), maximo= max(x), mediana = median(x), media = mean(x), DS = sd (x))})
head(summary_stats_genes) [1:5,1:5]

#En las columnas tengo muestras
summary_stats_muestras <- apply(matriz,2, function(x) {
  c(minimo = min(x), maximo= max(x), mediana = median(x), media = mean(x), DS = sd (x))})
head(summary_stats_muestras) [1:5,1:5]


#Grafico de dispersion de las intensidades globales (media) para ver la calidad de las muestras
sample_intensities <-  apply(matriz,2,function(x) media = mean(x))

df_gse_data <- data.frame(Muestras = colnames(matriz),             #Para extraer nombre de columnas(muestras)
                          Intensidad_media = sample_intensities)  #Para extraer valores calculados de la media

grafico_dispersion <- ggplot(df_gse_data, aes(x = Muestras, y = Intensidad_media)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title = "Intensidades Globales de las Muestras", x = "Muestras", y= "Intensidad Media")

print(grafico_dispersion)


# Ordenar para ver las más bajas y más altas
sort(sample_intensities)[1:5]      # 5 muestras con media más baja
tail(sort(sample_intensities), 5)  # 5 muestras con media más alta

#Calculo el percentil 90 de expresionde los genes, para usar como threshold y calcular cuantos genes estan por encima o debajo de ese valor
gene_means <- apply(matriz,1,function(x) Media = mean(x))
rowMeans(matriz)
threshold <- quantile(gene_means, probs= 0.9)
high_expr_genes <- gene_means > threshold
low_expr_genes <-  gene_means < threshold
sum(high_expr_genes)
sum(low_expr_genes)


# Distribución de intensidades promedio por muestra, agrupada por condición (Intensidad media por muestra)
sample_intensities <- colMeans(matriz)

df_intensidades <- data.frame(
  Muestra   = colnames(matriz),
  Intensidad= sample_intensities,
  Condicion = sample_info$condition
)

ggplot(df_intensidades, aes(x = Condicion, y = Intensidad, color = Condicion)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = "Intensidad promedio por condición")

############################################################
#                     MAS GRAFICOS
############################################################

#No tengo que hacer transformación logaritmoca, la matriz ya esta transformada
counts_data <- matriz


library (reshape2)

#Creo un boxplot con ggplot
counts_df <-  as.data.frame(counts_data)

counts_df$Gene <- rownames(counts_data)

counts_long <- reshape2 :: melt(counts_df, id.vars = "Gene", variable.name = "Sample",value.name = "Expression")


ggplot(counts_long, aes(x= Sample, y= Expression,)) +
  geom_boxplot()+
  theme_minimal()+
  labs(x = "Muestras", y= "Valores de Expresión", title = "Boxplot de Counts Data") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(counts_long, aes(x= Expression,))+
  geom_histogram(bins = 50, fill ="red", color = "green")+
  labs(title= "Distribución de la Expresión Génica", x = "Valores de Expresión", y= "Frecuencia")





############################################################
#                     PCA
############################################################

#Inicialmente me dio error:
#Error in PCAtools::pca(counts_data, metadata = sample_info) : 
#  'colnames(mat)' is not identical to 'rownames(metadata)'.

colnames(counts_data)

rownames(sample_info)

#Asi hago la correción para que coincidan
rownames(sample_info) <- sample_info$sample_label

identical(colnames(counts_data), rownames(sample_info))

#Me da nuevamente error, por los nombres duplicados
#Error in `.rowNamesDF<-`(x, value = value) : 
#duplicate 'row.names' are not allowed


# Alinear metadata por muestras (columnas)
rownames(sample_info) <- sample_info$sample_label
stopifnot(identical(colnames(counts_data), rownames(sample_info)))

# Calculo PCA
fit <- PCAtools::pca(counts_data, metadata = sample_info)

# Grafico PCA por tiempo post-lesión
library(ggplot2)

pca_plot <- PCAtools::biplot(
  fit,
  lab = sample_info$sample_label,   # etiquetas con los IDs de muestra
  colby = "condition",              # colores por condición (CTL, Acute, 3mpi, 6mpi, 12mpi)
  legendPosition = "right"
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")

print(pca_plot)


#Por Grupo (SCI / Control)
pca_plot_2 <- PCAtools::biplot(
  fit,
  lab = sample_info$sample_label,   # etiquetas con los IDs de muestra
  colby = "group",              # grupo (SCI / Control)
  legendPosition = "right"
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")

print(pca_plot_2)

######################################################################
#                     MAPA DE CALOR - De correlacion de muestras
######################################################################

comtx <- cor(counts_data)

print(comtx)

#Por grupo
annotation_df <- data.frame(group = sample_info$group)

rownames(annotation_df) <- row.names(sample_info)

pheatmap::pheatmap(comtx, annotation_col = annotation_df, main = "Mapa de correlaciones", angle_col = 45)

#Por condicion(agudo,cronico)
annotation_df_2 <- data.frame(group = sample_info$condition)

rownames(annotation_df_2) <- row.names(sample_info)

pheatmap::pheatmap(comtx, annotation_col = annotation_df_2, main = "Mapa de correlaciones", angle_col = 45)


######################################################################
#                     MAPA DE CALOR - De Expresion de Genes
######################################################################
#Selecciono los 50 genes con mayor variablidad
#Con función apply calculo la varianza por fila 
#######################################################
#                     MAPA DE CALOR - De Expresion de Genes
######################################################################
# Selecciono los 50 genes con mayor variabilidad
# Con función apply calculo la varianza por fila 

# --- Alinear y comprobar que sample_info coincide con las columnas de counts_data ---
stopifnot(all(colnames(counts_data) %in% sample_info$sample_label))

# Reordeno sample_info para que el orden sea exactamente el de las columnas de counts_data
sample_info <- sample_info[match(colnames(counts_data), sample_info$sample_label), , drop = FALSE]
stopifnot(identical(colnames(counts_data), sample_info$sample_label))

# Me aseguro de que annotation_col tenga el mismo orden que las columnas de la matriz
annotation_col <- data.frame(Condicion = sample_info$condition, row.names = sample_info$sample_label)

# 1) Calculo la variabilidad de cada transcripto (varianza por fila)
#    Uso na.rm = TRUE por seguridad ante posibles NA
gene_variability <- apply(counts_data, 1, function(x) var(x, na.rm = TRUE))

# 2) Selecciono los 50 transcriptos con mayor variabilidad (sin pasarme del total de filas)
k <- min(50L, nrow(counts_data))
top_genes <- names(sort(gene_variability, decreasing = TRUE))[seq_len(k)]
print(gene_variability[top_genes])

# 3) Extraigo esos transcriptos de la matriz de conteo
toplot <- counts_data[top_genes, , drop = FALSE]

# 4) Paso a data.frame para poder agregar anotaciones
toplot_df <- as.data.frame(toplot, stringsAsFactors = FALSE)

# 5) Agrego columna transcript_id (los rownames originales de la matriz)
toplot_df$transcript_id <- rownames(toplot_df)

# 6) Uno con el diccionario (transcripto ↔ gen) para agregar el símbolo de gen
#    Uso match para conservar el orden de 'toplot_df' y evitar reordenamientos del merge
#    'dict' debe tener columnas: transcript_id y gene_symbol
toplot_df$gene_symbol <- dict$gene_symbol[match(toplot_df$transcript_id, dict$transcript_id)]

# Ahora toplot_df tiene:
#   ... columnas de expresión por muestra ... | transcript_id | gene_symbol
# Muestro las primeras columnas (ajustá el rango si querés ver más/menos)
head(toplot_df[, c("transcript_id", "gene_symbol"), drop = FALSE])

# 7) Armo matriz para el heatmap:
#    - Selecciono exactamente las columnas de muestra en el mismo orden que counts_data
#    - Uso transcript_id como rownames para asegurar unicidad
sample_cols <- colnames(counts_data)
toplot_matrix <- as.matrix(toplot_df[, sample_cols, drop = FALSE])
rownames(toplot_matrix) <- toplot_df$transcript_id

# 8) Opcional: preparo anotaciones de fila para mostrar símbolos de gen en vez de transcriptos
row_annot <- data.frame(Gene = toplot_df$gene_symbol, row.names = toplot_df$transcript_id, check.names = FALSE)

# Comprobación clave: las columnas de la matriz deben coincidir con la anotación de columnas
stopifnot(identical(colnames(toplot_matrix), rownames(annotation_col)))

# 9) Hago el mapa de calor con pheatmap
pheatmap::pheatmap(
  toplot_matrix,
  annotation_row = row_annot,         
  annotation_col = annotation_col,    
  main = "Top 50 transcriptos más variables",
  scale = "row",                      # estandariza por fila para comparar patrones
  show_rownames = FALSE,              
  cluster_cols = TRUE, 
  cluster_rows = TRUE
)



####################################################################################
####################################################################################
#                  ANALISIS DE EXPRESION DIFERENCIAL   — LIMMA 
#####################################################################################

library(limma)

# 1) Creo la matriz de diseño sin intercepto
desing <- model.matrix(~ 0 + condition, data = sample_info)  
colnames(desing) <- levels(sample_info$condition)            

# Limpio nombres: si empiezan con números, make.names() les agrega "X"
colnames(desing) <- make.names(colnames(desing))
print(colnames(desing))
# → "CTL" "Acute" "X3mpi" "X6mpi" "X12mpi"

# Renombro para que quede claro
colnames(desing) <- c("CTL", "Acute", "X3mpi", "X6mpi", "X12mpi")

# Ajusto modelo lineal inicial
fit  <- lmFit(counts_data, desing)
head(fit$coefficients)


# 2) Defino CONTRASTES múltiples (comparaciones entre grupos)
contrast_matrix <- makeContrasts(
  # Comparaciones con control
  Acute_vs_CTL = Acute - CTL,
  M3_vs_CTL    = X3mpi - CTL,
  M6_vs_CTL    = X6mpi - CTL,
  M12_vs_CTL   = X12mpi - CTL,
  
  # Comparaciones entre tiempos  
  M3_vs_Acute  = X3mpi - Acute,
  M6_vs_Acute  = X6mpi - Acute,
  M12_vs_Acute = X12mpi - Acute,
  
  M6_vs_M3     = X6mpi - X3mpi,
  M12_vs_M3    = X12mpi - X3mpi,
  M12_vs_M6    = X12mpi - X6mpi,
  
  levels = desing
)

# Aplico contrastes y estadística bayesiana
fit2_allcontrasts <- eBayes(contrasts.fit(fit, contrast_matrix))


############################################################
# Conteo de genes UP y DOWN por contraste vs CTL
############################################################

# Defino el umbral como en el paper
lfc_thr <- log2(1.5)  # ≈ 0.58
padj_thr <- 0.05      # o 0.001 si querés más estricto

# Lista de contrastes a revisar
contrastes <- c("Acute_vs_CTL", "M3_vs_CTL", "M6_vs_CTL", "M12_vs_CTL")

# Función para contar genes
contar_genes <- function(contraste) {
  res <- topTable(fit2_allcontrasts, coef = contraste, number = Inf, sort.by = "P")
  up   <- sum(res$logFC >=  lfc_thr & res$adj.P.Val <= padj_thr, na.rm = TRUE)
  down <- sum(res$logFC <= -lfc_thr & res$adj.P.Val <= padj_thr, na.rm = TRUE)
  data.frame(Contraste = contraste, UP = up, DOWN = down)
}

# Aplico a cada contraste
tabla_resumen <- do.call(rbind, lapply(contrastes, contar_genes))

print(tabla_resumen)



###  ACA VOY MODIFICANDO LOS CONTRASTES ######
# Extraigo resultados completos de un contraste puntual (ejemplo: M12_vs_CTL)
res <- topTable(fit2_allcontrasts, coef = "M12_vs_CTL", number = Inf, sort.by = "P")

# Agrego transcript_id como columna
res$transcript_id <- rownames(res)

# Uno con el diccionario para tener el gene_symbol
res <- merge(res, dict, by = "transcript_id", all.x = TRUE)

# Filtrar duplicados: quedarme con el transcripto con menor p-valor por gen
library(dplyr)
res_filtrado <- res %>%
  group_by(gene_symbol) %>%
  slice_min(order_by = adj.P.Val, with_ties = FALSE) %>%
  ungroup()

# Ordenar por p-valor
res_filtrado <- res_filtrado[order(res_filtrado$adj.P.Val), ]

# Guardar resultados filtrados
write.csv(res_filtrado, "DEG_filtrados_por_pvalor.csv", row.names = FALSE)




# 3) CONTRASTE colapsado CTL vs SCI_all
#    Colapso todos los SCI en un solo grupo
sample_info$condition_merged <- factor(
  ifelse(sample_info$group == "CTL", "CTL", "SCI_all"),
  levels = c("CTL","SCI_all")
)

# Defino sujeto (ID del individuo) para bloquear medidas repetidas
sample_info$subject <- sub("v.*$", "", sample_info$sample_label)

# Matriz de diseño sin intercepto (~0 + factor)
design2 <- model.matrix(~ 0 + condition_merged, data = sample_info)
colnames(design2) <- levels(sample_info$condition_merged)  # c("CTL","SCI_all")

# Estimo correlación intra-sujeto (importante para repetidas medidas)
corfit2 <- duplicateCorrelation(counts_data, design = design2, block = sample_info$subject)

# Ajusto modelo lineal con bloqueo por sujeto
fit_blocked <- lmFit(
  counts_data,
  design2,
  block       = sample_info$subject,
  correlation = corfit2$consensus.correlation
)

# Defino contraste SCI_all vs CTL
contr <- makeContrasts(SCI_all_vs_CTL = SCI_all - CTL, levels = design2)

# Aplico contraste y estadística bayesiana
fit2_collapsed <- eBayes(contrasts.fit(fit_blocked, contr))

# Extraigo tabla de resultados del contraste colapsado
tt <- topTable(fit2_collapsed, coef = "SCI_all_vs_CTL", number = Inf)

# Marco dirección biológica (up / down / NS)
lfc_thr  <- 1
fdr_thr  <- 0.05
tt$direction <- ifelse(tt$adj.P.Val <= fdr_thr & tt$logFC >=  lfc_thr, "UP",
                       ifelse(tt$adj.P.Val <= fdr_thr & tt$logFC <= -lfc_thr, "DOWN", "NS"))

# Guardo resultados
write.csv(tt, "DEG_SCI_all_vs_CTL.csv", row.names = FALSE)


############################################################
#                     — VOLCANOS
############################################################
colnames(fit2_allcontrasts$coefficients)


############################################################
# Volcano con criterios del paper (Morrison et al. 2023)
############################################################

# Elegí el contraste
contraste <- "M12_vs_CTL"

# Criterios del paper
lfc_thr <- log2(1.5)   # ≈ 0.58
padj_thr <- 0.05       

# Extraigo resultados completos del contraste
res_paper <- topTable(fit2_allcontrasts, coef = contraste, number = Inf, sort.by = "P")

# Agrego transcript_id
res_paper$transcript_id <- rownames(res_paper)

# Uno con el diccionario de genes
res_paper <- merge(res_paper, dict, by = "transcript_id", all.x = TRUE)

# Calculo -log10(FDR)
res_paper$neglog10 <- -log10(pmax(res_paper$adj.P.Val, .Machine$double.xmin))

# Clasificación según criterios
res_paper$clase <- "NS"
res_paper$clase[res_paper$logFC >=  lfc_thr & res_paper$adj.P.Val <= padj_thr] <- "UP"
res_paper$clase[res_paper$logFC <= -lfc_thr & res_paper$adj.P.Val <= padj_thr] <- "DOWN"

# Genes a rotular (top 15 por FDR)
to_label <- res_paper[order(res_paper$adj.P.Val), ]
to_label <- head(to_label, 15)

# Volcano plot
p_volcano_paper <- ggplot(res_paper, aes(x = logFC, y = neglog10, color = clase)) +
  geom_point(alpha = 0.6, size = 1.3) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = to_label,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    title = paste("Volcano (criterios Morrison et al.) —", contraste),
    x = "log2 Fold Change",
    y = expression(-log[10]("(FDR)")),
    color = "Regulación"
  ) +
  theme_minimal()

print(p_volcano_paper)


############################################################
# Volcano exploratorio (menos estricto)
############################################################

# Elegí el contraste que quieras mirar
contraste <- "M12_vs_CTL"

# Umbrales más suaves
lfc_thr <- 0.5       # antes usabas 1
padj_thr <- 0.1      # antes usabas 0.05

# Extraigo resultados sin filtrar transcriptos
res_expl <- topTable(fit2_allcontrasts, coef = contraste, number = Inf, sort.by = "P")

# Agrego transcript_id como columna
res_expl$transcript_id <- rownames(res_expl)

# Uno con el diccionario para ver símbolos de gen
res_expl <- merge(res_expl, dict, by = "transcript_id", all.x = TRUE)

# Calculo -log10 del FDR
res_expl$neglog10 <- -log10(pmax(res_expl$adj.P.Val, .Machine$double.xmin))

# Clasificación menos estricta
res_expl$clase <- "NS"
res_expl$clase[res_expl$logFC >=  lfc_thr & res_expl$adj.P.Val <= padj_thr] <- "UP"
res_expl$clase[res_expl$logFC <= -lfc_thr & res_expl$adj.P.Val <= padj_thr] <- "DOWN"

# Genes a rotular (top 20 por FDR)
to_label <- res_expl[order(res_expl$adj.P.Val), ]
to_label <- head(to_label, 20)

# Volcano plot
p_volcano_expl <- ggplot(res_expl, aes(x = logFC, y = neglog10, color = clase)) +
  geom_point(alpha = 0.6, size = 1.3) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = to_label,
    aes(label = gene_symbol),
    size = 2.8,
    max.overlaps = 20
  ) +
  labs(
    title = paste("Volcano exploratorio —", contraste),
    x = "log2 Fold Change",
    y = expression(-log[10]("(adj. p-value)")),
    color = "Regulación"
  ) +
  theme_minimal()

print(p_volcano_expl)


############################################################
############################################################
# GSEA enfocado en Inflamación / CV / Hueso (robusto)
# - Usa fit2_allcontrasts y dict transcript→gene
# - Filtra MSigDB por keywords para 3 buckets temáticos
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(limma)
})

# --------- 1) Ranking por gen para un contraste dado (tu misma lógica) ---------
build_rank_vector <- function(coef_name,
                              fit_obj = fit2_allcontrasts,
                              dict_tbl = dict,
                              metric = c("t","logFC")) {
  metric <- match.arg(metric)
  res <- limma::topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P")
  res$transcript_id <- rownames(res)
  res$gene_symbol   <- dict_tbl$gene_symbol[match(res$transcript_id, dict_tbl$transcript_id)]
  
  # Colapso a nivel gen: me quedo con el transcripto de menor FDR por símbolo
  res_gene <- res %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = adj.P.Val, with_ties = FALSE) %>%
    ungroup()
  
  rank_vec <- if (metric == "t") res_gene$t else res_gene$logFC
  names(rank_vec) <- res_gene$gene_symbol
  rank_vec <- rank_vec[!is.na(rank_vec)]
  
  # Si hubiera duplicados (raro), me quedo con el de mayor |valor|
  if (any(duplicated(names(rank_vec)))) {
    ord <- order(abs(rank_vec), decreasing = TRUE)
    rank_vec <- rank_vec[ord]
    rank_vec <- rank_vec[!duplicated(names(rank_vec))]
  }
  sort(rank_vec, decreasing = TRUE)
}

# --------- 2) Cargar MSigDB (robusto) ---------
get_msig_tbl <- function(collections = c("H","C2","C5"), include_c7 = FALSE) {
  cols <- collections
  if (include_c7) cols <- c(cols, "C7")
  out <- dplyr::bind_rows(lapply(cols, function(cat) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = cat)
  }))
  message("MSigDB cargado: ", nrow(out), " filas; columnas: ", paste(colnames(out), collapse=", "))
  out
}

# --------- 3) Filtro por keywords (robusto a ausencia de gs_cat/gs_subcat) ---------
filter_msig_by_keywords <- function(msig_tbl, keywords, restrict_subcats = NULL) {
  stopifnot(is.data.frame(msig_tbl), "gs_name" %in% colnames(msig_tbl), "gene_symbol" %in% colnames(msig_tbl))
  tbl <- msig_tbl
  
  has_cat    <- "gs_cat"    %in% colnames(tbl)
  has_subcat <- "gs_subcat" %in% colnames(tbl)
  
  # Restringir por subcategorías solo si las columnas existen
  if (!is.null(restrict_subcats) && has_cat && has_subcat) {
    tbl <- tbl %>%
      dplyr::filter(gs_cat %in% names(restrict_subcats) &
                      gs_subcat %in% unlist(restrict_subcats))
  }
  
  name_lc <- tolower(tbl$gs_name)
  patt <- paste0(keywords, collapse = "|")
  tbl2 <- tbl[grepl(patt, name_lc, perl = TRUE), ]
  
  dplyr::distinct(tbl2, gs_name, gene_symbol, .keep_all = TRUE)
}

# --------- 4) Runner de GSEA por bucket (usa el filtro robusto) ---------
run_gsea_bucket <- function(rank_vec, msig_tbl, keywords, bucket_name,
                            minGS = 15, maxGS = 500, p_cut = 0.05, save_prefix = NULL,
                            restrict_subcats = NULL) {
  
  msig_f <- filter_msig_by_keywords(msig_tbl, keywords, restrict_subcats)
  
  if (nrow(msig_f) == 0) {
    message("Sin conjuntos para bucket: ", bucket_name)
    return(invisible(NULL))
  }
  
  t2g <- msig_f %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::distinct()
  
  gsea_res <- clusterProfiler::GSEA(
    geneList      = rank_vec,
    TERM2GENE     = t2g,
    minGSSize     = minGS,
    maxGSSize     = maxGS,
    pvalueCutoff  = p_cut,
    pAdjustMethod = "BH",
    eps           = 0,
    seed          = TRUE
  )
  
  if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    message("GSEA sin términos significativos para bucket: ", bucket_name)
    return(invisible(NULL))
  }
  
  tab <- gsea_res@result %>% dplyr::arrange(p.adjust)
  print(utils::head(tab, 10))
  
  if (!is.null(save_prefix)) {
    utils::write.csv(tab, paste0(save_prefix, "_", bucket_name, "_GSEA.csv"), row.names = FALSE)
  }
  
  ups   <- tab %>% dplyr::filter(NES > 0) %>% dplyr::arrange(dplyr::desc(NES)) %>% dplyr::pull(ID) %>% utils::head(2)
  downs <- tab %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES)                  %>% dplyr::pull(ID) %>% utils::head(2)
  
  for (id in ups) {
    p <- enrichplot::gseaplot(gsea_res, geneSetID = id, title = paste0(bucket_name, " — ", id))
    print(p)
    if (!is.null(save_prefix)) ggsave(paste0(save_prefix, "_", bucket_name, "_UP_", id, ".png"),
                                      p, width = 7, height = 5, dpi = 300)
  }
  for (id in downs) {
    p <- enrichplot::gseaplot(gsea_res, geneSetID = id, title = paste0(bucket_name, " — ", id))
    print(p)
    if (!is.null(save_prefix)) ggsave(paste0(save_prefix, "_", bucket_name, "_DOWN_", id, ".png"),
                                      p, width = 7, height = 5, dpi = 300)
  }
  
  toplot <- tab %>%
    dplyr::mutate(Direction = ifelse(NES >= 0, "UP", "DOWN"),
                  Term = factor(Description))
  p_bubble <- ggplot2::ggplot(toplot, ggplot2::aes(x = Direction, y = Term, size = abs(NES), color = p.adjust)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_continuous(name = "FDR", trans = "reverse") +
    ggplot2::labs(title = paste("GSEA —", bucket_name), x = "", y = "") +
    ggplot2::theme_light()
  print(p_bubble)
  if (!is.null(save_prefix)) ggsave(paste0(save_prefix, "_", bucket_name, "_bubble.png"),
                                    p_bubble, width = 8, height = 10, dpi = 300)
  
  invisible(list(result = gsea_res, table = tab,
                 ups = ups, downs = downs, bubble = p_bubble))
}

# --------- 5) Buckets temáticos (ajustá libremente estas listas) ---------
kw_inflam <- c(
  "inflamm", "nfk", "tnfa", "il1", "il2", "il6", "jak", "stat",
  "interferon", "ifna", "ifng", "toll", "tlr", "inflammasome",
  "cytokine", "chemokine", "complement", "acute_phase"
)

kw_cardio <- c(
  "cholesterol", "lipid", "lipoprotein", "fatty",
  "atherosclerosis", "endotheli", "angiogen",
  "platelet", "coagulation", "hemostasis",
  "smooth_muscle", "vascul", "oxidative_phosphorylation", "reactive_oxygen"
)

kw_bone <- c(
  "bone", "osteoclast", "osteoblast", "ossification",
  "resorption", "mineralization", "rankl", "rank", "opg"
)

# (Opcional) restringir subcategorías si EXISTEN en tu msig_tbl
# - C2: CP:REACTOME / CP:WIKIPATHWAYS
# - C5: GO:BP
restrict_subcats <- list(
  "C2" = c("CP:REACTOME", "CP:WIKIPATHWAYS"),
  "C5" = c("GO:BP")
)

# ===================== EJECUCIÓN =====================

# Elegí el contraste (como en colnames(fit2_allcontrasts$coefficients))
contraste_gsea <- "M12_vs_CTL"
# Métrico para rankear (recomendado "t")
metric <- "t"

# 1) Ranking
rank_vec <- build_rank_vector(contraste_gsea, metric = metric)

# 2) MSigDB combinado
msig_tbl <- get_msig_tbl(collections = c("H","C2","C5"), include_c7 = TRUE)  # poné FALSE si no querés C7

# 3) Correr por bucket (pasando restrict_subcats para usarlo si gs_cat/gs_subcat existen)
out_inflam <- run_gsea_bucket(rank_vec, msig_tbl, kw_inflam, "Inflamacion",
                              minGS = 25, maxGS = 500, p_cut = 0.05,
                              save_prefix = paste0("GSEA_", contraste_gsea, "_", metric),
                              restrict_subcats = restrict_subcats)

out_cardio <- run_gsea_bucket(rank_vec, msig_tbl, kw_cardio, "Cardio-Metabolico",
                              minGS = 25, maxGS = 500, p_cut = 0.05,
                              save_prefix = paste0("GSEA_", contraste_gsea, "_", metric),
                              restrict_subcats = restrict_subcats)

out_bone   <- run_gsea_bucket(rank_vec, msig_tbl, kw_bone, "Hueso-Osteoporosis",
                              minGS = 10, maxGS = 500, p_cut = 0.05,  # en hueso a veces los sets son más chicos
                              save_prefix = paste0("GSEA_", contraste_gsea, "_", metric),
                              restrict_subcats = restrict_subcats)

# (Opcional) inspección rápida de columnas disponibles para depurar:
# colnames(msig_tbl); if ("gs_cat" %in% colnames(msig_tbl)) unique(msig_tbl$gs_cat)
# if ("gs_subcat" %in% colnames(msig_tbl)) unique(msig_tbl$gs_subcat)





############################################################
# LECTURA Y RESUMEN DE GSEA (compatible clusterProfiler/fgsea)
# Requiere que ya existan: rank_vec, out_inflam, out_cardio, out_bone
# y el objeto fit2_allcontrasts (para validar el contraste).
############################################################

## === 0) Comprobaciones mínimas ===
stopifnot(exists("fit2_allcontrasts"))
stopifnot(exists("rank_vec"))
contraste_gsea <- "M12_vs_CTL"  # tu contraste (típico: M12 - CTL)
if (!("t" %in% names(fit2_allcontrasts))) {
  stop("No encuentro fit2_allcontrasts$t. Asegurate de tener el objeto de limma con t-values.")
}
stopifnot(contraste_gsea %in% colnames(fit2_allcontrasts$t))
cat("Contraste usado:", contraste_gsea, "\n")
cat("Regla de interpretación: NES > 0 ⇒ enriquecido en M12; NES < 0 ⇒ enriquecido en CTL.\n\n")

## === 1) Función robusta para convertir resultados a una tabla estándar ===
normalize_gsea_tbl <- function(x) {
  # Caso clusterProfiler::gseResult guardado en x$gsea
  if (!is.null(x$gsea)) {
    tab <- as.data.frame(x$gsea@result)
    out <- data.frame(
      pathway     = if ("ID" %in% names(tab)) tab$ID else tab$Description,
      description = if ("Description" %in% names(tab)) tab$Description else tab$ID,
      size        = if ("setSize" %in% names(tab)) tab$setSize else NA_integer_,
      NES         = tab$NES,
      padj        = if ("p.adjust" %in% names(tab)) tab$p.adjust else tab$qvalues,
      leadingEdge = if ("core_enrichment" %in% names(tab)) tab$core_enrichment else NA_character_,
      stringsAsFactors = FALSE
    )
    return(out)
  }
  # Caso fgsea: data.frame típico en x$tbl
  if (!is.null(x$tbl)) {
    tab <- x$tbl
    if (!"description" %in% names(tab)) tab$description <- tab$pathway
    size_col <- intersect(c("size","setSize"), names(tab))
    padj_col <- intersect(c("padj","p.adjust","qvalues"), names(tab))
    le_col   <- intersect(c("leadingEdge","core_enrichment"), names(tab))
    out <- data.frame(
      pathway     = tab$pathway,
      description = tab$description,
      size        = if (length(size_col)) tab[[size_col[1]]] else NA_integer_,
      NES         = tab$NES,
      padj        = if (length(padj_col)) tab[[padj_col[1]]] else NA_real_,
      leadingEdge = if (length(le_col))   tab[[le_col[1]]]   else NA_character_,
      stringsAsFactors = FALSE
    )
    return(out)
  }
  stop("Estructura no reconocida: el objeto no tiene $gsea ni $tbl.")
}

## === 2) Normalizar tus tres buckets ===
stopifnot(exists("out_inflam"), exists("out_cardio"), exists("out_bone"))
tbl_inflam <- normalize_gsea_tbl(out_inflam)
tbl_cardio <- normalize_gsea_tbl(out_cardio)
tbl_bone   <- normalize_gsea_tbl(out_bone)

## === 3) Helpers para ordenar/filtrar y mostrar ===
ord_top <- function(df, n = 10) {
  if (!nrow(df)) return(df)
  df[order(df$padj, df$NES), c("pathway","NES","padj")][1:min(nrow(df), n), ]
}

peek_sets <- function(tbl, ids) {
  hit <- tbl[tbl$pathway %in% ids | tbl$description %in% ids, c("pathway","NES","padj")]
  if (!nrow(hit)) return(hit)
  hit[order(hit$padj, hit$NES), ]
}

## === 4) Resumen rápido por bucket ===
cat("=== INFLAMACIÓN (top por padj) ===\n")
print(ord_top(tbl_inflam, 10)); cat("\n")

cat("=== CARDIO-METABÓLICO (top por padj) ===\n")
print(ord_top(tbl_cardio, 10)); cat("\n")

cat("=== HUESO-OSTEOPOROSIS (top por padj) ===\n")
print(ord_top(tbl_bone, 10)); cat("\n")

## === 5) Chequeo dirigido de Hallmarks cardio-inflamatorios clave ===
hallmarks_cv <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COAGULATION",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_FATTY_ACID_METABOLISM"
)

cat(">>> Hallmarks cardio-inflamatorios (buscados en INFLAMACIÓN):\n")
print(peek_sets(tbl_inflam, hallmarks_cv)); cat("\n")

cat(">>> Hallmarks cardio-inflamatorios (buscados en CARDIO-METABÓLICO):\n")
print(peek_sets(tbl_cardio, hallmarks_cv)); cat("\n")

## === 6) Leading edge del pathway de tu gráfico (ajustá si querés otro) ===
pw_target <- "GSE34156_UNTREATED_VS_24H_NOD2_AND_TLR1_TLR2_LIGAND_TREATED_MONOCYTE_DN"

fila <- tbl_inflam[tbl_inflam$pathway == pw_target | tbl_inflam$description == pw_target, ]
where <- "INFLAMACIÓN"
if (nrow(fila) == 0) {
  fila <- tbl_cardio[tbl_cardio$pathway == pw_target | tbl_cardio$description == pw_target, ]
  where <- "CARDIO-METABÓLICO"
}
if (nrow(fila) == 0) {
  fila <- tbl_bone[tbl_bone$pathway == pw_target | tbl_bone$description == pw_target, ]
  where <- "HUESO"
}

if (nrow(fila) > 0) {
  cat("Pathway objetivo encontrado en bucket:", where, "\n")
  print(fila[, c("pathway","NES","padj")])
  
  # Parseo robusto del leading edge (puede venir como 'g1/g2/...' o como lista o vector)
  le_raw <- fila$leadingEdge[1]
  le_genes <- character(0)
  if (is.list(le_raw)) {
    le_genes <- unique(unlist(le_raw))
  } else if (!is.na(le_raw)) {
    s <- as.character(le_raw)
    if (grepl("/", s, fixed = TRUE)) {
      le_genes <- unlist(strsplit(s, "/", fixed = TRUE))
    } else if (grepl(",", s, fixed = TRUE)) {
      le_genes <- trimws(unlist(strsplit(s, ",", fixed = TRUE)))
    } else {
      le_genes <- unique(unlist(strsplit(s, "[ \t;]+")))
    }
  }
  le_genes <- sort(unique(le_genes[le_genes != ""]))
  cat("\nLeading edge (", length(le_genes), " genes):\n", paste(le_genes, collapse = ", "), "\n\n", sep = "")
  
  # (Opcional) recuperar signos de cada gen en tu ranking para validar dirección
  common <- intersect(names(rank_vec), le_genes)
  if (length(common)) {
    signos <- rank_vec[common]
    cat("Resumen de signos en ranking (t de ", contraste_gsea, "):\n", sep = "")
    print(summary(signos))
    cat("Top 10 más negativos del leading edge:\n")
    print(sort(signos, decreasing = FALSE)[1:min(10, length(signos))])
  }
} else {
  cat("No encontré el pathway objetivo en los resultados normalizados.\n")
}

## === 7) (Opcional) Graficar el enriquecimiento si usaste clusterProfiler ===
suppressWarnings({
  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
      !is.null(out_inflam$gsea) &&
      pw_target %in% out_inflam$gsea@result$ID) {
    p <- clusterProfiler::gseaplot2(out_inflam$gsea, geneSetID = pw_target, title = pw_target)
    print(p)
  }
})

############################################################
# Cómo leer los resultados:
# - NES > 0 ⇒ vía enriquecida en M12 (lado positivo del ranking M12-CTL).
# - NES < 0 ⇒ vía enriquecida en CTL (lado negativo del ranking).
# - 'padj' es el FDR: cuanto más bajo, más significativo.
# - El leading edge es el núcleo de genes que más contribuye a la señal.
############################################################



