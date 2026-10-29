# =============================================================================
#  OE2 · Distribución potencial de Jubaea chilensis, Quillaja saponaria y
#        Vachellia caven mediante ensamble de modelos (biomod2)
# -----------------------------------------------------------------------------
#  Autor      : Diego Salazar Canales
#  Tesis      : Magíster en Geografía, Universidad de Chile (2026)
#  Postulación: Premio MapBiomas Chile, 1ª edición
#
#  Flujo (en el orden de la metodología, apartado 4.2.2):
#    1. Descarga de ocurrencias ............ (manual, documentada)
#    2. Depuración de ocurrencias .......... se ejecuta
#    3. Preprocesamiento de predictores .... opcional (regenera rásteres)
#    4. Correlaciones (dos etapas) ......... se ejecuta
#    5. Modelamiento biomod2 ............... opcional (≈ 1 día de cómputo)
#
#  Requisitos: R >= 4.3; terra, sf, dplyr, readr, caret, pheatmap,
#              biomod2 4.3-4-5 (+ randomForest, maxnet; Java y maxent.jar)
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(readr)
library(caret)
library(pheatmap)

# -----------------------------------------------------------------------------
# 0) CONFIGURACIÓN
# -----------------------------------------------------------------------------
# Interruptores de los bloques pesados
ejecutar_preprocesamiento <- FALSE   # TRUE regenera los predictores a 250 m
ejecutar_modelacion       <- FALSE   # TRUE vuelve a correr biomod2 (≈ 1 día)

# Parámetros de la metodología
fecha_minima       <- as.Date("2000-01-01")  # filtro temporal de registros
dist_thinning_m    <- 500                    # distancia mínima entre registros (m)
umbral_cribado     <- 0.75                   # cribado automático (caret)
umbral_final       <- 0.70                   # criterio de selección final |r|
ratio_pa           <- 2                      # pseudoausencias por presencia
n_sets_pa          <- 2                      # sets independientes (2 × 1:1)
algoritmos         <- c("GLM", "RF", "MAXENT")
n_replicas         <- 10
prop_entrenamiento <- 0.75
semilla            <- 123

# Rutas (relativas a la carpeta OE2/)
dir_gbif      <- "datos/gbif"                    # descargas GBIF Chile (.csv, tab)
ruta_pobladas <- "datos/areas_pobladas/AREAS_POBLADAS_WGS194_19S.shp"  # BCN, mapas vectoriales
ruta_area     <- "datos/area_descarga/AREA_DESCARGA.shp"
ruta_ref250   <- "datos/referencia/ref_raster250m_19s.tif"
ruta_ref500   <- "datos/referencia/ref_raster500m_19s.tif"
dir_brutos    <- "datos/predictores_brutos"      # worldclim/, topografia/, heterogeneidad/, soilgrids/
dir_250       <- "datos/predictores_250m"        # salida del bloque 3 (misma estructura)
dir_finales   <- "datos/predictores_finales"     # 8 capas seleccionadas por especie
dir_res       <- "resultados"
ruta_maxent   <- "maxent.jar"

dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)

especies <- data.frame(
  id     = c("Jchilensis", "Qsaponaria", "Vcaven"),
  nombre = c("Jubaea chilensis", "Quillaja saponaria", "Vachellia caven"),
  gbif   = file.path(dir_gbif, c("J_chilensis.csv", "Q_saponaria.csv", "V_caven.csv")),
  stringsAsFactors = FALSE
)

# Predictores finales por especie (códigos del Cuadro 3; resultado del bloque 4).
# Las capas correspondientes se encuentran en datos/predictores_finales/<especie>/
# y son las que usan los bloques 4 y 5.
predictores_finales <- list(
  Jchilensis = c("V1", "V17", "V21", "V27", "V33", "V34", "V35", "V38"),
  Qsaponaria = c("V3", "V14", "V23", "V26", "V28", "V30", "V35", "V38"),
  Vcaven     = c("V5", "V15", "V21", "V26", "V28", "V30", "V35", "V38")
)
sets_predictores <- c("worldclim", "soilgrids", "topografia", "heterogeneidad")


# =============================================================================
# 1) DESCARGA DE OCURRENCIAS (manual)
# =============================================================================
# Los registros de cada especie se descargaron del portal GBIF Chile
# (https://gbifchile.mma.gob.cl/) para todo el territorio nacional, en formato
# de texto delimitado por tabulaciones. Los archivos originales se guardan en
# datos/gbif/ sin modificaciones.

# =============================================================================
# 2) DEPURACIÓN DE OCURRENCIAS
# =============================================================================
ref250 <- rast(ruta_ref250)
areas_pobladas <- st_read(ruta_pobladas, quiet = TRUE) |>
  st_transform(32719) |>
  st_make_valid()
areas_pobladas <- areas_pobladas[st_geometry_type(areas_pobladas) %in% c("POLYGON", "MULTIPOLYGON"), ]

# 2.1 Depuración tabular: duplicados, registros anteriores a 2000 o sin fecha
#     válida, sin identificación a nivel de especie y sin coordenadas
depurar_gbif <- function(ruta) {
  df <- read_delim(ruta, delim = "\t", escape_double = FALSE, trim_ws = TRUE,
                   show_col_types = FALSE)
  n_gbif <- nrow(df)
  df <- df[!duplicated(df), ]
  df$eventDate <- as.Date(df$eventDate)
  df <- df[!is.na(df$eventDate) & df$eventDate > fecha_minima, ]
  df <- df[!is.na(df$species), ]
  df <- df[!is.na(df$decimalLongitude) & !is.na(df$decimalLatitude),
           c("scientificName", "decimalLatitude", "decimalLongitude", "eventDate")]
  df$ID <- seq_len(nrow(df))
  attr(df, "n_gbif") <- n_gbif
  df
}

# 2.2 Filtro espacial: fuera de áreas pobladas, dentro de la grilla de 250 m
#     y un único registro por celda
filtrar_espacial <- function(df) {
  pts <- st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326, remove = FALSE) |> st_transform(32719)
  n_inicial <- nrow(pts)
  
  fuera <- lengths(st_intersects(pts, areas_pobladas)) == 0
  pts <- pts[fuera, ]
  n_fuera_pobladas <- nrow(pts)
  
  v <- vect(pts)
  pts$celda <- cellFromXY(ref250, crds(v))
  pts$valor <- terra::extract(ref250, v)[, 2]
  pts <- pts[!is.na(pts$celda) & !is.na(pts$valor), ]
  pts <- pts |> arrange(ID) |> group_by(celda) |> slice(1) |> ungroup()
  
  attr(pts, "conteos") <- c(antes_espacial = n_inicial,
                            fuera_pobladas = n_fuera_pobladas,
                            una_por_celda = nrow(pts))
  pts
}

# 2.3 Spatial thinning por distancia mínima: se recorren los registros en orden
#     y se descartan los que quedan a menos de dist_min de uno ya conservado
thinning_distancia <- function(pts, dist_min) {
  pts <- pts[order(pts$ID), ]
  n <- nrow(pts)
  conservar <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!conservar[i]) next
    d <- as.numeric(st_distance(pts[i, ], pts))
    conservar[which(d < dist_min & seq_len(n) > i)] <- FALSE
  }
  pts[conservar, ]
}

ocurrencias <- list()
resumen_registros <- data.frame()

for (i in seq_len(nrow(especies))) {
  id <- especies$id[i]
  df <- depurar_gbif(especies$gbif[i])
  pts <- filtrar_espacial(df)
  cont <- attr(pts, "conteos")
  pts <- thinning_distancia(pts, dist_thinning_m)
  
  ocurrencias[[id]] <- pts
  write_csv(st_drop_geometry(pts), file.path(dir_res, paste0("OE2_ocurrencias_finales_", id, ".csv")))
  
  resumen_registros <- rbind(resumen_registros, data.frame(
    especie = especies$nombre[i],
    registros_GBIF = attr(df, "n_gbif"),
    desde_2000_con_especie = nrow(df),
    fuera_areas_pobladas = cont[["fuera_pobladas"]],
    un_registro_por_celda_250m = cont[["una_por_celda"]],
    thinning_final = nrow(pts)
  ))
}
print(resumen_registros, row.names = FALSE)
write_csv(resumen_registros, file.path(dir_res, "OE2_resumen_registros_por_etapa.csv"))

# =============================================================================
# 3) PREPROCESAMIENTO DE PREDICTORES
# =============================================================================
# Variables de 1 km (WorldClim 2.1, EarthEnv topografía y heterogeneidad):
# recorte → reproyección a WGS 84 / UTM 19S → remuestreo bilineal a 500 m →
# remuestreo bilineal a 250 m sobre la grilla de SoilGrids (paquete terra).
# SoilGrids (250 m): solo se reproyecta y alinea a la grilla de referencia.
#
# Pasos previos (fuera de R):
#  - SoilGrids se descargó desde Google Earth Engine (8 variables × 6 intervalos
#    de profundidad, 0–200 cm) y se promedió entre profundidades.
#  - El índice de humedad topográfica (V34) se derivó del DEM de EarthEnv.
preparar_predictor <- function(ruta, set, area, ref500, ref250) {
  r <- rast(ruta)
  if (set == "soilgrids") {
    if (!same.crs(r, ref250)) r <- project(r, ref250, method = "bilinear")
    r <- resample(r, ref250, method = "bilinear")
  } else {
    a <- if (same.crs(area, r)) area else project(area, crs(r))
    r <- mask(crop(r, a), a)
    r <- project(r, crs(ref250), method = "bilinear")
    r <- mask(resample(r, ref500, method = "bilinear"), ref500)
    r <- resample(r, ref250, method = "bilinear")
  }
  crop(mask(r, ref250), ref250, snap = "near")
}

if (ejecutar_preprocesamiento) {
  area   <- vect(ruta_area)
  ref500 <- rast(ruta_ref500)
  for (set in sets_predictores) {
    dir.create(file.path(dir_250, set), showWarnings = FALSE, recursive = TRUE)
    for (f in list.files(file.path(dir_brutos, set), "\\.tif$", full.names = TRUE)) {
      r <- preparar_predictor(f, set, area, ref500, ref250)
      writeRaster(r, file.path(dir_250, set, basename(f)), overwrite = TRUE)
    }
  }
}

# =============================================================================
# 4) CORRELACIONES EN DOS ETAPAS
# =============================================================================
# La correlación de Pearson se calcula con los valores de los predictores
# extraídos en los puntos de ocurrencia de cada especie.
#   Etapa 1: dentro de cada set; cribado automático |r| > 0,75 (caret) y
#            selección final de dos predictores por set (|r| < 0,7 y
#            relevancia ecológica; selección manual).
#   Etapa 2: correlación global entre los predictores preseleccionados; de cada
#            par con |r| >= 0,7 se eliminó uno según su relevancia ecológica
#            (selección manual). Las capas finales (máx. 8 por especie) se
#            copian a datos/predictores_finales/<especie>/ y aquí se documenta
#            su matriz de correlación.
extraer_valores <- function(pts, carpeta) {
  stk <- rast(list.files(carpeta, "\\.tif$", full.names = TRUE))
  v <- project(vect(pts), crs(stk))
  vals <- terra::extract(stk, v)[, -1, drop = FALSE]
  vals <- vals[, apply(vals, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  na.omit(vals)
}

guardar_heatmap <- function(M, archivo) {
  png(archivo, width = 900, height = 900, res = 150)
  pheatmap(M, display_numbers = TRUE, number_format = "%.2f",
           clustering_method = "complete", fontsize = 6, fontsize_number = 5,
           color = colorRampPalette(c("#2C7FB8", "#F7F7F7", "#D7191C"))(100),
           breaks = seq(-1, 1, length.out = 101))
  dev.off()
}

cribado <- data.frame()
for (id in especies$id) {
  # Etapa 1: por set ambiental
  for (set in sets_predictores) {
    vals <- extraer_valores(ocurrencias[[id]], file.path(dir_250, set))
    M <- cor(vals, method = "pearson")
    eliminadas <- colnames(M)[findCorrelation(M, cutoff = umbral_cribado)]
    write.csv(round(M, 4), file.path(dir_res, paste0("OE2_corr_", id, "_", set, ".csv")))
    guardar_heatmap(M, file.path(dir_res, paste0("OE2_corr_", id, "_", set, ".png")))
    cribado <- rbind(cribado, data.frame(
      especie = id, set = set, n_variables = ncol(M),
      eliminadas_cribado = paste(eliminadas, collapse = ", ")))
  }
  
  # Etapa 2: predictores finales
  vals <- extraer_valores(ocurrencias[[id]], file.path(dir_finales, id))
  M <- cor(vals, method = "pearson")
  write.csv(round(M, 4), file.path(dir_res, paste0("OE2_corr_", id, "_final.csv")))
  guardar_heatmap(M, file.path(dir_res, paste0("OE2_corr_", id, "_final.png")))
  cat(sprintf("%s: |r| máximo entre predictores finales = %.2f\n",
              id, max(abs(M[upper.tri(M)]))))
}
write_csv(cribado, file.path(dir_res, "OE2_cribado_correlacion.csv"))

# =============================================================================
# 5) MODELAMIENTO CON biomod2 4.3-4-5
# =============================================================================
if (ejecutar_modelacion) {
  library(biomod2)
  stopifnot(packageVersion("biomod2") == "4.3.4.5")
  options(java.parameters = "-Xmx4g")
  Sys.setenv(PATH_TO_MAXENT_JAR = normalizePath(ruta_maxent))
  set.seed(semilla)
  
  for (i in seq_len(nrow(especies))) {
    id  <- especies$id[i]
    stk <- rast(list.files(file.path(dir_finales, id), "\\.tif$", full.names = TRUE))
    stk <- resample(stk, ref250, method = "bilinear")
    mascara <- app(stk, function(x) as.numeric(all(!is.na(x))))
    mascara[mascara == 0] <- NA
    
    # 5.1 Presencias con valores ambientales completos
    pres <- project(vect(ocurrencias[[id]]), crs(stk))
    ok <- complete.cases(terra::extract(stk, pres)[, -1])
    xy_pres <- as.data.frame(crds(pres[ok, ]))
    
    # 5.2 Pseudoausencias aleatorias: 2 sets independientes de 1:1 (total 2:1),
    #     en celdas válidas no ocupadas por presencias
    celdas_validas <- which(!is.na(values(mascara, mat = FALSE)))
    celdas_libres  <- setdiff(celdas_validas, unique(cellFromXY(mascara, xy_pres)))
    n_por_set <- nrow(xy_pres) * ratio_pa / n_sets_pa
    celdas_pa <- c()
    for (k in seq_len(n_sets_pa)) {
      nuevas <- sample(setdiff(celdas_libres, celdas_pa), n_por_set)
      celdas_pa <- c(celdas_pa, nuevas)
    }
    xy_pa <- as.data.frame(xyFromCell(mascara, celdas_pa))
    names(xy_pa) <- names(xy_pres)
    
    # 5.3 Formateo
    datos <- BIOMOD_FormatingData(
      resp.name = id,
      resp.var  = c(rep(1, nrow(xy_pres)), rep(0, nrow(xy_pa))),
      resp.xy   = rbind(xy_pres, xy_pa),
      expl.var  = stk,
      dir.name  = file.path("datos/biomod2", id),
      data.type = "binary", PA.nb.rep = 0, na.rm = TRUE,
      filter.raster = FALSE, seed.val = semilla)
    
    # 5.4 Modelos individuales: GLM, RF y Maxent × 10 réplicas (75/25),
    #     con la configuración por defecto de biomod2
    modelos <- BIOMOD_Modeling(
      bm.format = datos, modeling.id = paste0(id, "_Modeling"),
      models = algoritmos, CV.strategy = "random",
      CV.nb.rep = n_replicas, CV.perc = prop_entrenamiento,
      var.import = 3, metric.eval = c("AUCroc", "TSS"),
      scale.models = FALSE, OPT.strategy = "default",
      nb.cpu = 1, seed.val = semilla)
    
    # 5.5 Ensamble por réplica: EMwmean de los 3 algoritmos, ponderado por su
    #     AUC de calibración
    ensamble <- BIOMOD_EnsembleModeling(
      bm.mod = modelos, models.chosen = "all", em.by = "PA+run",
      em.algo = c("EMmean", "EMwmean"),   # EMmean se calculó, pero solo se usó EMwmean
      metric.select = "AUCroc", metric.select.thresh = 0,
      metric.select.dataset = "calibration",
      metric.eval = c("AUCroc", "TSS"), var.import = 3,
      EMwmean.decay = "proportional", nb.cpu = 1, seed.val = semilla)
    
    # 5.6 Evaluación (AUC y TSS sobre los datos de validación de cada réplica)
    #     e importancia de variables
    write_csv(get_evaluations(modelos),
              file.path(dir_res, paste0("OE2_evaluacion_individuales_", id, ".csv")))
    write_csv(get_evaluations(ensamble),
              file.path(dir_res, paste0("OE2_evaluacion_ensamble_", id, ".csv")))
    write_csv(get_variables_importance(ensamble),
              file.path(dir_res, paste0("OE2_importancia_variables_", id, ".csv")))
    
    # 5.7 Proyección de los 10 ensambles EMwmean y promedio simple entre réplicas
    em_wmean <- grep("EMwmean", get_built_models(ensamble), value = TRUE)
    proy <- BIOMOD_EnsembleForecasting(
      bm.em = ensamble, proj.name = paste0(id, "_EMwmean_RUNS"),
      new.env = mask(stk, mascara), models.chosen = em_wmean,
      metric.binary = NULL, metric.filter = NULL, na.rm = TRUE, nb.cpu = 1)
    
    preds <- get_predictions(proy)
    capas_runs <- grep("EMwmeanByAUCroc_allData_RUN([1-9]|10)_mergedAlgo$", names(preds), value = TRUE)
    idoneidad <- app(preds[[capas_runs]], mean, na.rm = TRUE)
    # biomod2 entrega la idoneidad en escala 0–1000; se lleva a 0–1
    if (global(idoneidad, "max", na.rm = TRUE)[1, 1] > 1) idoneidad <- idoneidad / 1000
    
    writeRaster(idoneidad, file.path(dir_res, paste0("OE2_idoneidad_", id, ".tif")), overwrite = TRUE)
    
    # 5.8 Reclasificación en 5 categorías de igual amplitud y superficie por categoría
    clases <- classify(idoneidad, rcl = matrix(c(0, .2, 1, .2, .4, 2, .4, .6, 3,
                                                 .6, .8, 4, .8, 1.01, 5), ncol = 3, byrow = TRUE),
                       include.lowest = TRUE)
    writeRaster(clases, file.path(dir_res, paste0("OE2_idoneidad_clases_", id, ".tif")), overwrite = TRUE)
    
    sup <- as.data.frame(freq(clases))
    sup$area_km2 <- sup$count * prod(res(clases)) / 1e6
    write_csv(sup, file.path(dir_res, paste0("OE2_superficie_categorias_", id, ".csv")))
  }
}

