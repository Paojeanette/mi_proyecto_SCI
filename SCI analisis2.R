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

# 3) Dejo solo el símbolo del gen (quitar lo que está antes de ">" que es el identificador de transcripto RefSeq mRNA)
#    Ej: "NM_000016>ACADM" -> "ACADM"
colnames(df) <- sub(".*>", "", colnames(df))

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

#Vector de símbolos tal como aparecen en la matriz
syms <- rownames(matriz_gene)


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
map[map$SYMBOL == "ACADM", ]


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

# Colapsar duplicados por símbolo base (sin sufijos .1, .2, ...)
#Creo un ID base por gen quitando sufijos numéricos del final del nombre de fila.
#Ejemplos: RPL8.1 → RPL8, TMUB1.2 → TMUB1.
#Resultado: tenia RPL8, RPL8.1, RPL8.2, ahora queda una sola fila RPL8 con el promedio de esas filas en cada muestra.
id_base <- sub("\\.\\d+$", "", rownames(counts_data))
counts_data <- limma::avereps(counts_data, ID = id_base, fun = mean)
stopifnot(!any(duplicated(rownames(counts_data))))



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

head(counts_data)
gene_variability <-  apply(counts_data,1,function(x) var(x))

top_genes <- sort(gene_variability, decreasing = TRUE)[1:50]
print(top_genes)

# Selecciono de counts_data las filas(genes) con mayor variabilidaad
toplot <- counts_data[names(top_genes), ]

#sincronizar diccionario con los genes presentes en counts_data
map <- map[match(rownames(counts_data), map$SYMBOL), , drop = FALSE]
rownames(map) <- map$SYMBOL

# Verifico qué nombres devuelve map (es como annotation_data (fData))
colnames(map)
# [1] "SYMBOL" "ENTREZID" "ENSEMBL" "GENENAME"

# rownames en 'map' 
rownames(map) <- map$SYMBOL

# Renombro filas de toplot con los nombres de gen oficiales
rownames(toplot) <- map[names(top_genes), "SYMBOL"]

head(toplot[, 1:5])  # primeras filas y columnas

# VER LOS DUPLICADOS, TENGO QUE VOLVER A VERLO PARA CORREGIRLO
# YA SEA DESDE COUNT_DATA O DESDE MAP O AMBOS !!!

#Como es una matriz, busco asi los nombres de los genes en la matriz

#Vector de expresión de RPL8 
counts_data["RPL8", ]   # si el nombre es exactamente "RPL8"

counts_data["TMUB1",]

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

# Ejemplo: ver resultados para un contraste puntual
head(topTable(fit2_allcontrasts, coef = "Acute_vs_CTL", number = 10))


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



# 1) Volcano para contraste puntual (ej: "Acute_vs_CTL")
contraste <- "M12_vs_CTL"   # elegir el contraste que quieras
lfc_thr   <- 1                 
padj_thr  <- 0.05

res1 <- topTable(fit2_allcontrasts, coef = contraste, number = Inf, sort.by = "P")
res1$gen <- rownames(res1)
res1$neglog10 <- -log10(pmax(res1$adj.P.Val, .Machine$double.xmin))

# Clasifico puntos
res1$clase <- "NS"
res1$clase[res1$logFC >=  lfc_thr & res1$adj.P.Val <= padj_thr] <- "UP"
res1$clase[res1$logFC <= -lfc_thr & res1$adj.P.Val <= padj_thr] <- "DOWN"

# Elijo genes a rotular
to_label <- res1[res1$clase != "NS", ]
to_label <- to_label[order(to_label$neglog10, decreasing = TRUE), ]
to_label <- head(to_label, 15)   # top 15

p_volcano1 <- ggplot(res1, aes(x = logFC, y = neglog10, color = clase)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = to_label,
    aes(label = gen),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = paste("Volcano —", contraste),
    x = "log2 Fold Change",
    y = expression(-log[10]("(adj. p-value)")),
    color = "Regulacion"
  ) +
  theme_minimal()

print(p_volcano1)
ggsave(paste0("volcano_", contraste, ".png"), p_volcano1, width = 7, height = 5, dpi = 300)


# 2) Volcano para contraste colapsado (SCI_all vs CTL)
volc <- within(tt, {
  neglog10 <- -log10(pmax(adj.P.Val, .Machine$double.xmin))
  gene     <- rownames(tt)
})

p_volcano2 <- ggplot(volc, aes(x = logFC, y = neglog10, color = direction)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = subset(volc, direction != "NS")[1:min(15, sum(volc$direction != "NS")), ],
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    title = "Volcano — SCI_all vs CTL",
    x = "log2 Fold Change",
    y = expression(-log[10]("(adj. p-value)")),
    color = "Regulación"
  ) +
  theme_minimal()

print(p_volcano2)
ggsave("volcano_SCIall_vs_CTL.png", p_volcano2, width = 7, height = 5, dpi = 300)

