# =============================================================================
#  OE1 · Caracterización y proyección (2042) de coberturas de suelo
#  Preparación de drivers, correlación (Anexo 1) y validación 2022 (Cuadro 19)
# -----------------------------------------------------------------------------
#  Autor   : Diego Salazar Canales
#  Tesis   : Magíster en Geografía, Universidad de Chile (2026)
#  Postulación: Premio MapBiomas Chile, 1ª edición
#
#  Qué hace este script (pasos realizados en R):
#    1. Calcula el percentil 95 (P95) de cada driver.
#    2. Normaliza los drivers a 0–1 con la expresión usada en la tesis
#       (equivalente en R de la expresión aplicada en ArcGIS Pro).
#    3. Alinea los rásteres a la grilla de referencia (30 m) para TerrSet.
#    4. Calcula la matriz de correlación de Pearson entre drivers (Anexo 1).
#    5. Calcula Hits, Misses, False Alarms, Correct rejections y el
#       Figure of Merit de la validación 2022 (Cuadro 19).
#
#  Pasos realizados FUERA de R (documentados en el README):
#    - Descarga y reclasificación de MapBiomas Chile Colección 1 (ArcGIS Pro).
#    - Distancias euclidianas de los drivers de proximidad (ArcGIS Pro).
#    - Matrices de Márkov, submodelos, V de Cramér, MLP, simulación 2022 y
#      proyección 2042 (TerrSet, Land Change Modeler).
#    - Kappa global y por clase (TerrSet, módulo VALIDATE / CROSSTAB).
#
#  Requisitos: R >= 4.3; paquetes terra, psych, reshape2, ggplot2, ggcorrplot
# =============================================================================

library(terra)
library(psych)
library(reshape2)
library(ggplot2)
library(ggcorrplot)

# -----------------------------------------------------------------------------
# 0) CONFIGURACIÓN
# -----------------------------------------------------------------------------
# Todas las rutas son relativas a la carpeta OE1/ del repositorio.
dir_dist   <- "datos/drivers_distancia"    # distancias euclidianas brutas (m) y topografía
dir_norm   <- "datos/drivers_normalizados" # salida de la normalización P95
dir_final  <- "datos/drivers_finales"      # drivers alineados que ingresan al LCM
dir_res    <- "resultados"
ruta_ref   <- file.path(dir_final, "DEM_resample_final.tif")  # grilla de referencia (30 m)

dir.create(dir_norm,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_final, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_res,   showWarnings = FALSE, recursive = TRUE)

# Drivers finales del modelo (Cuadro 2): archivo final -> etiqueta
drivers <- c(
  "DEM_resample_final.tif"        = "Elevación",
  "SLOPE_resample_final.tif"      = "Pendiente",
  "dist_nucleos02_p95_final.tif"  = "Núcleos urbanos (2002)",
  "dist_vial02_p95_final.tif"     = "Red vial (2002)",
  "dist_LU94_p95_final.tif"       = "Límite urbano PRMS (1994)",
  "dist_erosion_p95_final.tif"    = "Área histórica erosionada",
  "dist_APB_p95_final.tif"        = "Áreas protegidas",
  "forestal02_resample2.tif"      = "Plantaciones existentes (2002)",
  "agri02_resample2.tif"          = "Agricultura-pasturas (2002)"
)
# Nota: núcleos urbanos, red vial, plantaciones y agricultura-pasturas (2002)
# se derivan de MapBiomas Chile Colección 1 (año 2002). La red vial corresponde
# a la red de caminos principales 2019 (rutas clase 1 y 3) restringida a la
# clase Infraestructura de MapBiomas 2002.


# -----------------------------------------------------------------------------
# 1) PERCENTIL 95 DE CADA DRIVER
# -----------------------------------------------------------------------------
# Para rásteres muy grandes, global(quantile) puede ser lento; en ese caso
# puede estimarse sobre una muestra con spatSample().
calcular_p95 <- function(r) {
  as.numeric(global(r, fun = quantile, probs = 0.95, na.rm = TRUE)[1, 1])
}


# -----------------------------------------------------------------------------
# 2) NORMALIZACIÓN 0–1 POR P95
# -----------------------------------------------------------------------------
# Expresiones aplicadas en ArcGIS Pro (Raster Calculator):
#   Drivers de proximidad:  Con(x > P95, P95, x) / P95
#       0 = sobre el elemento · 1 = a una distancia >= P95
#   Elevación y pendiente (invertidas):  1 - Con(x > P95, P95, x) / P95
#       1 = menor elevación / terreno plano · 0 = valores >= P95
normalizar_p95 <- function(r, p95, invertir = FALSE) {
  x <- clamp(r, upper = p95, values = TRUE) / p95
  if (invertir) 1 - x else x
}

# Uso: para cada driver x (distancias en m; DEM en m; pendiente en grados)
#   normalizar_p95(x, calcular_p95(x))                  # proximidad
#   normalizar_p95(x, calcular_p95(x), invertir = TRUE) # elevación y pendiente

# -----------------------------------------------------------------------------
# 3) ALINEACIÓN A LA GRILLA DE REFERENCIA (requisito de TerrSet/LCM)
# -----------------------------------------------------------------------------
# TerrSet exige que todos los rásteres compartan exactamente extensión,
# resolución y origen. Se remuestrea a la grilla del DEM (30 m) por vecino más
# cercano para no alterar los valores normalizados.
alinear_a_referencia <- function(archivo, ref, salida) {
  r <- resample(rast(archivo), ref, method = "near")
  writeRaster(r, salida, overwrite = TRUE, datatype = "FLT4S",
              NAflag = -9999, gdal = "COMPRESS=LZW")
  compareGeom(ref, rast(salida), stopOnError = FALSE)
}

# Uso: alinear_a_referencia("driver_p95.tif", rast(ruta_ref), "driver_p95_final.tif")

# -----------------------------------------------------------------------------
# 4) CORRELACIÓN DE PEARSON ENTRE DRIVERS (Anexo 1)
# -----------------------------------------------------------------------------
# Se correlacionan los 9 drivers finales (normalizados y alineados) sobre una
# muestra aleatoria de 100.000 píxeles. Los pares con |r| >= 0,7 se mantuvieron
# en el modelo, dado que el MLP tolera variables correlacionadas (Anexo 1).
stack_drivers <- rast(file.path(dir_final, names(drivers)))
names(stack_drivers) <- unname(drivers)

set.seed(42)
muestra <- spatSample(stack_drivers, size = 1e5, method = "random",
                      na.rm = TRUE, values = TRUE)
cat("Píxeles muestreados:", nrow(muestra), "\n")

M  <- cor(muestra, method = "pearson", use = "pairwise.complete.obs")
Pm <- psych::corr.test(as.matrix(muestra), use = "pairwise",
                       method = "pearson", adjust = "none")$p

write.csv(round(M, 6),  file.path(dir_res, "OE1_correlacion_pearson_drivers.csv"))
write.csv(round(Pm, 6), file.path(dir_res, "OE1_correlacion_pvalues_drivers.csv"))

# Pares con |r| >= 0,7 (documentación de la multicolinealidad)
M_long <- melt(M, varnames = c("Var1", "Var2"), value.name = "r")
M_long <- M_long[as.character(M_long$Var1) < as.character(M_long$Var2), ]
pares_altos <- M_long[abs(M_long$r) >= 0.7, ]
write.csv(pares_altos[order(-abs(pares_altos$r)), ],
          file.path(dir_res, "OE1_pares_correlacion_mayor_0.7.csv"), row.names = FALSE)

# Figura del Anexo 1
fig_corr <- ggcorrplot(M, p.mat = Pm, type = "upper", hc.order = TRUE,
                       outline.color = "white", lab = TRUE, lab_size = 3,
                       tl.cex = 10, tl.srt = 45, show.diag = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        plot.margin = margin(10, 20, 10, 10))

ggsave(file.path(dir_res, "OE1_Anexo1_correlacion_drivers.png"), fig_corr,
       width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 5) VALIDACIÓN 2022: HITS, MISSES, FALSE ALARMS Y FIGURE OF MERIT (Cuadro 19)
# -----------------------------------------------------------------------------
# Insumo: ráster de validación de la corrida final (F99) generado en TerrSet
# a partir de la tabulación cruzada 2012 (base) × 2022 simulado × 2022 observado.
# Códigos: 0 = Correct rejection · 1 = Hit · 2 = Miss · 3 = False alarm
ruta_val <- "datos/validacion/VAL_PRED22_F99_FoM.rst"

r_val <- rast(ruta_val)
levels(r_val) <- NULL   # trabajar con los códigos numéricos, no con etiquetas

frec <- as.data.frame(freq(r_val))[, c("value", "count")]
names(frec) <- c("codigo", "pixeles")

# Si el .rdc de TerrSet declara 0 como valor de fondo, terra lo lee como NA.
# En ese caso, los rechazos correctos se obtienen como el total de celdas
# menos Hits, Misses y False alarms (equivale a "Values < min. specified").
if (!0 %in% frec$codigo) {
  frec <- rbind(data.frame(codigo = 0, pixeles = ncell(r_val) - sum(frec$pixeles)), frec)
  message("Código 0 leído como NoData: rechazos correctos calculados por diferencia.")
}

etiquetas <- c("0" = "Correct rejection", "1" = "Hit", "2" = "Miss", "3" = "False alarm")
frec$categoria <- unname(etiquetas[as.character(frec$codigo)])

# Superficie con el tamaño real de celda (29,998 m en la grilla de TerrSet)
area_celda_km2 <- prod(res(r_val)) / 1e6
frec$area_km2   <- frec$pixeles * area_celda_km2
frec$porcentaje <- 100 * frec$pixeles / sum(frec$pixeles)

H  <- frec$pixeles[frec$categoria == "Hit"]
M_ <- frec$pixeles[frec$categoria == "Miss"]
FA <- frec$pixeles[frec$categoria == "False alarm"]
FoM <- H / (H + Mi + FA)   # Figure of Merit (Leta et al., 2021; Shoyama, 2021)

cuadro19 <- rbind(
  frec[order(frec$codigo), c("categoria", "pixeles", "area_km2", "porcentaje")],
  data.frame(categoria = "Total", pixeles = sum(frec$pixeles),
             area_km2 = sum(frec$area_km2), porcentaje = 100)
)
print(cuadro19, row.names = FALSE)
cat(sprintf("\nFigure of Merit: %.4f\n", FoM))

write.csv(cuadro19, file.path(dir_res, "OE1_Cuadro19_validacion_2022.csv"), row.names = FALSE)

