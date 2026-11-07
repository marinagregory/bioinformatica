##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

tabla <- read.table("archivo-entrada.txt",
                    header = TRUE,          # la primera fila son los nombres de columnas
                    sep = "\t",             # el archivo está separado por tabuladores
                    check.names = FALSE,    # respeta nombres como "Gene Name" o "1W5-TPM"
                    stringsAsFactors = FALSE)

head(tabla)   # primeras 6 filas
tail(tabla)   # últimas 6 filas
dim(tabla)    # filas x columnas: debe dar 619 x 7


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

# Columnas 2-4 = Control   (C2-TPM, C1-TPM, C3-TPM)
# Columnas 5-7 = Condición (1W5-TPM, 1W6-TPM, 1W5-COV)
control   <- tabla[, 2:4]
condicion <- tabla[, 5:7]

tabla$Media_Control   <- apply(control,   1, mean)   # el 1 indica "por filas" (cada gen)
tabla$DE_Control      <- apply(control,   1, sd)
tabla$Media_Condicion <- apply(condicion, 1, mean)
tabla$DE_Condicion    <- apply(condicion, 1, sd)

head(tabla)   # comprobamos que hay 4 columnas nuevas (11 en total)


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

# Colores que vamos a usar:
col_control   <- "steelblue"
col_condicion <- "tomato"

# Los TPM tienen valores muy dispares, así que se usa log2(media + 1)
# para que el histograma sea legible.
log_ctrl <- log2(tabla$Media_Control   + 1)
log_cond <- log2(tabla$Media_Condicion + 1)
rango    <- range(c(log_ctrl, log_cond))

par(mfrow = c(1, 2))   # 1 fila, 2 gráficos
hist(log_ctrl, col = col_control, xlim = rango,
     main = "Media del Control", xlab = "log2(Media TPM + 1)", ylab = "Número de genes")
hist(log_cond, col = col_condicion, xlim = rango,
     main = "Media de la Condición", xlab = "log2(Media TPM + 1)", ylab = "Número de genes")
dev.copy(png, "histogramas.png", width = 1000, height = 500)   # guarda la imagen
dev.off()
par(mfrow = c(1, 1))   # restablece la ventana gráfica


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

boxplot(log_ctrl, log_cond,
        names = c("Control", "Condición"),
        col   = c(col_control, col_condicion),
        main  = "Expresión media por gen: Control vs Condición",
        ylab  = "log2(Media TPM + 1)")
legend("topleft",
       legend = c("Control", "Condición"),
       fill   = c(col_control, col_condicion),
       title  = "Grupo")
dev.copy(png, "boxplot.png", width = 700, height = 600)
dev.off()


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

# Test de Shapiro-Wilk: H0 = los datos siguen una distribución normal.
# Si p-value < 0.05 se rechaza la normalidad.
shap_control   <- shapiro.test(tabla$Media_Control)
shap_condicion <- shapiro.test(tabla$Media_Condicion)
shap_control
shap_condicion

normal <- (shap_control$p.value > 0.05) && (shap_condicion$p.value > 0.05)

if (normal) {
  cat("Los datos son normales -> se usa t-test\n")
  tabla$Pvalue <- sapply(1:nrow(tabla), function(i) {
    tryCatch(t.test(as.numeric(control[i, ]), as.numeric(condicion[i, ]))$p.value,
             error = function(e) NA)
  })
} else {
  cat("Los datos NO son normales -> se usa el test no paramétrico de Wilcoxon\n")
  tabla$Pvalue <- sapply(1:nrow(tabla), function(i) {
    suppressWarnings(
      wilcox.test(as.numeric(control[i, ]), as.numeric(condicion[i, ]))$p.value
    )
  })
}

head(tabla)
summary(tabla$Pvalue)


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

DE <- tabla[!is.na(tabla$Pvalue) & tabla$Pvalue < 0.05, ]

cat("Número de genes diferencialmente expresados:", nrow(DE), "\n")
cat("p-value mínimo obtenido:", min(tabla$Pvalue, na.rm = TRUE), "\n")

if (nrow(DE) > 0) {
  sobreexpresados <- DE[DE$Media_Condicion > DE$Media_Control, ]
  reprimidos      <- DE[DE$Media_Condicion < DE$Media_Control, ]
  cat("Sobreexpresados:", nrow(sobreexpresados), "\n")
  cat("Reprimidos:",      nrow(reprimidos),      "\n")
  print(DE[, c(1, 8, 10, 12)])   # listado: gen, media control, media condición, p-value
} else {
  cat("No sale ningún gen con p < 0.05.\n")
}

# JUSTIFICACIÓN (si no sale ninguno):
# Los datos no son normales, así que se usa Wilcoxon, que es un test basado en rangos. Con solo
# 3 réplicas por grupo hay únicamente 20 ordenaciones posibles, de modo que el p-value
# mínimo posible del test exacto es 2/20 = 0.1. Con muchos ceros R usa una
# aproximación normal, y aun en la separación perfecta entre grupos el p-value queda
# en ~0.064. Por tanto ningún gen puede bajar de 0.05, sino que con n = 3 el test no 
# tiene potencia suficiente para detectarlas.


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

if (nrow(DE) > 0) {
  genes_hm <- DE
  titulo   <- "Heatmap de genes diferencialmente expresados"
} else {
  # Como no hay genes con p < 0.05, se muestran los 50 genes con menor p-value
  # y mayor diferencia de medias solo para ilustrar el heatmap. NO son significativos.
  orden    <- order(tabla$Pvalue, -abs(tabla$Media_Condicion - tabla$Media_Control))
  genes_hm <- tabla[orden[1:50], ]
  titulo   <- "Heatmap: 50 genes con menor p-value (ninguno significativo)"
}

mat <- as.matrix(genes_hm[, 2:7])
rownames(mat) <- genes_hm[, 1]
mat <- log2(mat + 1)
mat <- mat[apply(mat, 1, sd) > 0, ]   # quita genes constantes (no se pueden escalar)

heatmap(mat,
        scale         = "row",   # cada gen se escala para comparar patrones, no magnitudes
        col           = colorRampPalette(c("blue", "white", "red"))(100),
        ColSideColors = c(rep(col_control, 3), rep(col_condicion, 3)),
        main          = titulo,
        margins       = c(8, 8),
        cexRow        = 0.6)
dev.copy(png, "heatmap.png", width = 800, height = 900)
dev.off()

### Guardar la tabla final con todas las columnas nuevas
write.table(tabla, "tabla_resultados.txt", sep = "\t", quote = FALSE, row.names = FALSE)
