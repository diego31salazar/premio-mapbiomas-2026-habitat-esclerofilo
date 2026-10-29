# =============================================================================
#  OE3 · Estructura espacial del hábitat potencial de Jubaea chilensis,
#        Quillaja saponaria y Vachellia caven (2002–2042)
# -----------------------------------------------------------------------------
#  Autor      : Diego Salazar Canales
#  Tesis      : Magíster en Geografía, Universidad de Chile (2026)
#  Postulación: Premio MapBiomas Chile, 1ª edición
#
#  Flujo (en el orden de la metodología, apartado 4.2.3):
#    1. Delimitación del hábitat potencial ..... ArcGIS Pro (documentada)
#    2. Lectura y validación de los 36 rásteres
#    3. Métricas de paisaje CA y NP ............ -> Cuadros 13 y 14, Figuras 26 y 27
#    4. Modelos GLS con AR1 .................... -> Cuadros 15 y 16, Anexos 5 y 6
#    5. Índice PFI adaptado (Rivas et al. 2022)  -> Figuras 28–30
#
#  Entrada: 36 rásteres binarios de hábitat potencial (1 = hábitat, NA = resto),
#  250 m, WGS 84 / UTM 19S, nombrados <especie>_<año>_<nivel>.tif
#  (p. ej. Jchilensis_2042_NH3.tif): 3 especies × 4 años × 3 niveles.
#
#  Requisitos: R 4.6.1; terra 1.9-50, sf, landscapemetrics, nlme, dplyr, tidyr,
#              ggplot2, readr
# =============================================================================

library(terra)
library(sf)
library(landscapemetrics)
library(nlme)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# -----------------------------------------------------------------------------
# 0) CONFIGURACIÓN
# -----------------------------------------------------------------------------
# Interruptor
calcular_pfi <- FALSE   # TRUE: calcula el PFI por parche (≈ 1–2 h)

# Parámetros de la metodología
especies     <- c("Jchilensis", "Qsaponaria", "Vcaven")
anios        <- c(2002, 2012, 2022, 2042)
niveles      <- c("NH1", "NH2", "NH3")   # NH1 0,8–1 · NH2 0,6–1 · NH3 0,4–1
vecindad     <- 8                        # regla de vecindad de los parches
peso_espacio <- 4/5                      # PFI: componente 1 − Ap/Ai
peso_forma   <- 1/5                      # PFI: componente FRAC_patch / 2

# Rutas (relativas a la carpeta OE3/)
dir_habitat  <- "datos/habitat_potencial"            # 36 rásteres binarios
dir_res      <- "resultados"
dir_pfi      <- file.path(dir_res, "PFI")

dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)

etiquetas_sp <- c(Jchilensis = "J. chilensis",
                  Qsaponaria = "Q. saponaria",
                  Vcaven     = "V. caven")


# =============================================================================
# 1) DELIMITACIÓN DEL HÁBITAT POTENCIAL (ArcGIS Pro 3.7.1, documentada)
# =============================================================================
# Los 36 rásteres de entrada se construyeron fuera de R:
#   a) Niveles anidados de hábitat a partir del ráster de idoneidad del
#      ensamble de cada especie (OE2): NH1 = 0,8–1; NH2 = 0,6–1; NH3 = 0,4–1.
#   b) Máscara de bosque + matorral para 2002, 2012 y 2022 (MapBiomas Chile,
#      Colección 1) y 2042 (proyección del OE1).
#   c) Intersección de la máscara vectorial con una grilla de 250 m equivalente
#      al ráster de idoneidad; se conservaron las celdas con superficie de
#      intersección >= 62.499 m² (celdas completamente cubiertas) y se
#      rasterizaron.
#   d) Extract by Mask de cada nivel de hábitat con la máscara del año.
# Resultado: <especie>_<año>_<nivel>.tif en datos/habitat_potencial/.


# =============================================================================
# 2) LECTURA Y VALIDACIÓN DE LOS RÁSTERES
# =============================================================================
archivos <- list.files(dir_habitat, pattern = "\\.tif$", full.names = TRUE)

tabla_rasteres <- data.frame(archivo = archivos, stringsAsFactors = FALSE) %>%
  mutate(
    nombre        = tools::file_path_sans_ext(basename(archivo)),
    especie       = sapply(strsplit(nombre, "_"), `[`, 1),
    anio          = as.integer(sapply(strsplit(nombre, "_"), `[`, 2)),
    nivel_habitat = sapply(strsplit(nombre, "_"), `[`, 3)
  ) %>%
  filter(especie %in% especies, anio %in% anios, nivel_habitat %in% niveles) %>%
  arrange(especie, anio, nivel_habitat)

stopifnot("Se esperan 36 rásteres (3 especies × 4 años × 3 niveles)" =
            nrow(tabla_rasteres) == 36)

rasteres <- lapply(tabla_rasteres$archivo, rast)

# Misma geometría (resolución, extensión y CRS) y valores 1 / NA
stopifnot(all(sapply(rasteres, function(r) compareGeom(r, rasteres[[1]]))))
binario_ok <- sapply(rasteres, function(r) {
  v <- unique(values(r))
  all(v[!is.na(v)] == 1)
})
stopifnot("Algún ráster contiene valores distintos de 1 y NA" = all(binario_ok))

cat("Validación OK: 36 rásteres binarios,",
    paste(res(rasteres[[1]]), collapse = " × "), "m\n")


# =============================================================================
# 3) MÉTRICAS DE PAISAJE: CANTIDAD DE HÁBITAT (CA) Y NÚMERO DE PARCHES (NP)
# =============================================================================
# CA (lsm_c_ca, en ha; se expresa en km²) y NP (lsm_c_np) con landscapemetrics,
# regla de 8 vecinos. Las celdas NA se recodifican a 0 (matriz) y se reporta
# la clase 1 (hábitat).
metricas_raster <- function(r) {
  r_lsm <- ifel(is.na(r), 0, r)
  lsm <- calculate_lsm(r_lsm, what = c("lsm_c_ca", "lsm_c_np"),
                       directions = vecindad) %>%
    filter(class == 1)
  
  data.frame(
    area_km2  = lsm$value[lsm$metric == "ca"] / 100,
    n_parches = lsm$value[lsm$metric == "np"]
  )
}

tabla_metricas <- bind_cols(
  tabla_rasteres %>% select(especie, anio, nivel_habitat),
  bind_rows(lapply(rasteres, metricas_raster))
)

write_csv(tabla_metricas, file.path(dir_res, "OE3_tabla_metricas.csv"))

# Cuadros 13 y 14 (formato de la tesis: especie × nivel en filas, años en columnas)
cuadro_13 <- tabla_metricas %>%
  select(especie, nivel_habitat, anio, area_km2) %>%
  pivot_wider(names_from = anio, values_from = area_km2) %>%
  arrange(especie, nivel_habitat)

cuadro_14 <- tabla_metricas %>%
  select(especie, nivel_habitat, anio, n_parches) %>%
  pivot_wider(names_from = anio, values_from = n_parches) %>%
  arrange(especie, nivel_habitat)

write_csv(cuadro_13, file.path(dir_res, "OE3_cuadro13_cantidad_habitat_km2.csv"))
write_csv(cuadro_14, file.path(dir_res, "OE3_cuadro14_numero_parches.csv"))

cat("\nCuadro 13. Cantidad de hábitat potencial (km²)\n")
print(as.data.frame(cuadro_13), digits = 7)
cat("\nCuadro 14. Número de parches potenciales\n")
print(as.data.frame(cuadro_14))

# Figuras 26 y 27: trayectorias observadas por especie y nivel de hábitat
datos_fig <- tabla_metricas %>%
  mutate(anio          = factor(anio, levels = anios),
         nivel_habitat = factor(nivel_habitat, levels = niveles),
         especie_lab   = etiquetas_sp[especie])

graficar_trayectoria <- function(variable, eje_y) {
  ggplot(datos_fig, aes(anio, .data[[variable]],
                        color = nivel_habitat, group = nivel_habitat)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ especie_lab, scales = "free_y") +
    labs(x = "Año", y = eje_y, color = "Nivel de hábitat") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "italic", size = 12))
}

ggsave(file.path(dir_res, "OE3_figura26_trayectorias_CA.png"),
       graficar_trayectoria("area_km2", "Cantidad de hábitat potencial (km²)"),
       width = 11, height = 6, dpi = 300)
ggsave(file.path(dir_res, "OE3_figura27_trayectorias_NP.png"),
       graficar_trayectoria("n_parches", "Número de parches potenciales"),
       width = 11, height = 6, dpi = 300)


# =============================================================================
# 4) MÍNIMOS CUADRADOS GENERALIZADOS (GLS) CON ESTRUCTURA AR1
# =============================================================================
# Diseño: 36 observaciones = 3 especies × 3 niveles × 4 años.
#  - Año y nivel de hábitat como factores ordenados -> contrastes polinomiales
#    ortogonales (año: L, Q, C; nivel: L, Q).
#  - Especie con contraste de tratamiento (referencia: J. chilensis).
#  - Efectos aditivos: respuesta ~ año + nivel + especie.
#  - AR1 sobre los pasos sucesivos de la serie (2002, 2012, 2022, 2042 = 0–3)
#    dentro de cada serie especie × nivel de hábitat.
#  - Estimación por REML; IC 95 % con nlme::intervals().
datos_modelo <- tabla_metricas %>%
  mutate(
    especie        = factor(especie, levels = especies),
    anio_factor    = factor(anio, levels = anios, ordered = TRUE),
    nivel_habitat  = factor(nivel_habitat, levels = niveles, ordered = TRUE),
    tiempo         = match(anio, anios) - 1,
    serie_temporal = interaction(especie, nivel_habitat, drop = TRUE)
  ) %>%
  arrange(especie, nivel_habitat, anio)

ajustar_gls <- function(formula) {
  gls(formula, data = datos_modelo,
      correlation = corAR1(form = ~ tiempo | serie_temporal),
      method = "REML")
}

modelos <- list(
  CA = ajustar_gls(area_km2  ~ anio_factor + nivel_habitat + especie),
  NP = ajustar_gls(n_parches ~ anio_factor + nivel_habitat + especie)
)

# Cuadros 15 y 16: estimaciones, IC 95 % y p-valor
nombres_efectos <- c(
  "anio_factor.L"     = "Año (L)",
  "anio_factor.Q"     = "Año (Q)",
  "anio_factor^3"     = "Año (C)",
  "nivel_habitat.L"   = "Nivel hábitat (L)",
  "nivel_habitat.Q"   = "Nivel hábitat (Q)",
  "especieQsaponaria" = "Quillaja saponaria",
  "especieVcaven"     = "Vachellia caven"
)

tabla_coeficientes <- function(modelo) {
  tt <- summary(modelo)$tTable
  ic <- intervals(modelo, which = "coef")$coef
  data.frame(
    parametro  = rownames(tt),
    efecto     = unname(nombres_efectos[rownames(tt)]),
    estimacion = tt[, "Value"],
    ic_inf     = ic[, "lower"],
    ic_sup     = ic[, "upper"],
    p_valor    = tt[, "p-value"],
    row.names  = NULL
  ) %>%
    filter(!is.na(efecto))
}

cuadro_15 <- tabla_coeficientes(modelos$CA)
cuadro_16 <- tabla_coeficientes(modelos$NP)

write_csv(cuadro_15, file.path(dir_res, "OE3_cuadro15_GLS_cantidad_habitat.csv"))
write_csv(cuadro_16, file.path(dir_res, "OE3_cuadro16_GLS_numero_parches.csv"))

cat("\nCuadro 15. GLS cantidad de hábitat\n");  print(cuadro_15[, -1], digits = 5)
cat("\nCuadro 16. GLS número de parches\n");    print(cuadro_16[, -1], digits = 5)

# Anexos 5 y 6: diagnóstico de residuos normalizados (Cuadros 22 y 23)
diagnostico <- function(modelo) {
  e <- residuals(modelo, type = "normalized")
  sw <- shapiro.test(e)
  data.frame(
    shapiro_W = unname(sw$statistic),
    shapiro_p = sw$p.value,
    acf_lag1  = acf(e, plot = FALSE)$acf[2],
    phi_AR1   = unname(coef(modelo$modelStruct$corStruct, unconstrained = FALSE))
  )
}

tabla_diagnostico <- bind_rows(lapply(modelos, diagnostico), .id = "modelo")
write_csv(tabla_diagnostico, file.path(dir_res, "OE3_anexos5y6_diagnostico_GLS.csv"))
cat("\nDiagnóstico de residuos (Cuadros 22 y 23)\n")
print(tabla_diagnostico, digits = 3)

# Gráficos de residuos de los Anexos 5 (CA) y 6 (NP)
for (m in names(modelos)) {
  df <- data.frame(ajustados = fitted(modelos[[m]]),
                   residuos  = residuals(modelos[[m]], type = "normalized"))
  
  p_res <- ggplot(df, aes(ajustados, residuos)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(x = "Valores ajustados", y = "Residuos normalizados") +
    theme_minimal(base_size = 13)
  
  p_qq <- ggplot(df, aes(sample = residuos)) +
    stat_qq(size = 2.3) +
    stat_qq_line(color = "red", linewidth = 1) +
    labs(x = "Cuantiles teóricos", y = "Cuantiles observados") +
    theme_minimal(base_size = 13)
  
  ggsave(file.path(dir_res, paste0("OE3_anexo_residuos_vs_ajustados_", m, ".png")),
         p_res, width = 8, height = 5, dpi = 300)
  ggsave(file.path(dir_res, paste0("OE3_anexo_qqplot_", m, ".png")),
         p_qq, width = 8, height = 5, dpi = 300)
}


# =============================================================================
# 5) ÍNDICE DE FRAGMENTACIÓN POR PARCHE (PFI adaptado de Rivas et al., 2022)
# =============================================================================
#   PFI = 4/5 · (1 − Ap/Ai) + 1/5 · (FRAC_patch / 2)
#   Ap         = área del parche
#   Ai         = área de influencia: polígono de Voronoi generado a partir de un
#                punto interior de cada parche (st_point_on_surface) y recortado
#                a la extensión del ráster; la razón Ap/Ai se acota a [0, 1]
#   FRAC_patch = 2 · ln(0,25 · P) / ln(A)  (dimensión fractal del parche,
#                misma fórmula que lsm_p_frac de landscapemetrics; P en m, A en m²)
# Los rásteres continuos del PFI para 2042 son la base de las Figuras 28–30
# (composición cartográfica en ArcGIS Pro).
if (calcular_pfi) {
  
  dir.create(dir_pfi, showWarnings = FALSE, recursive = TRUE)
  
  calcular_pfi_raster <- function(r, especie, anio, nivel) {
    
    cat("PFI:", especie, anio, nivel, "\n")
    
    # Parches (8 vecinos) y polígonos
    r_id <- patches(r, directions = vecindad)
    names(r_id) <- "patch_id"
    parches <- st_as_sf(as.polygons(r_id, dissolve = TRUE, na.rm = TRUE)) %>%
      st_make_valid() %>%
      filter(!is.na(patch_id))
    
    # Área de influencia (Voronoi recortado a la extensión del ráster)
    limite <- st_as_sfc(st_bbox(c(xmin = xmin(r), ymin = ymin(r),
                                  xmax = xmax(r), ymax = ymax(r)),
                                crs = st_crs(parches)))
    puntos <- st_point_on_surface(parches) %>% select(patch_id)
    
    if (nrow(puntos) == 1) {
      voronoi <- st_sf(patch_id = puntos$patch_id, geometry = limite)
    } else {
      voronoi <- st_voronoi(st_union(st_geometry(puntos))) %>%
        st_collection_extract("POLYGON") %>%
        st_sf(geometry = ., crs = st_crs(parches)) %>%
        st_make_valid() %>%
        st_intersection(st_sf(geometry = limite)) %>%
        st_make_valid() %>%
        st_join(puntos, join = st_intersects, left = FALSE) %>%
        group_by(patch_id) %>%
        summarise(geometry = st_union(geometry), .groups = "drop")
    }
    
    ai <- data.frame(patch_id = voronoi$patch_id,
                     Ai_m2    = as.numeric(st_area(voronoi)))
    
    # Ap, perímetro, FRAC y PFI
    tabla <- data.frame(
      patch_id = parches$patch_id,
      Ap_m2    = as.numeric(st_area(parches)),
      perim_m  = as.numeric(st_length(st_boundary(parches)))
    ) %>%
      left_join(ai, by = "patch_id") %>%
      mutate(
        Ap_Ai      = pmin(pmax(Ap_m2 / Ai_m2, 0), 1),
        FRAC_patch = 2 * log(0.25 * perim_m) / log(Ap_m2),
        FRAC_patch = ifelse(is.finite(FRAC_patch), FRAC_patch, NA),
        PFI        = peso_espacio * (1 - Ap_Ai) + peso_forma * (FRAC_patch / 2),
        especie = especie, anio = anio, nivel_habitat = nivel
      )
    
    # Ráster continuo del PFI
    r_pfi <- subst(r_id, from = tabla$patch_id, to = tabla$PFI, others = NA)
    writeRaster(r_pfi,
                file.path(dir_pfi, paste0(especie, "_", anio, "_", nivel, "_PFI.tif")),
                overwrite = TRUE, datatype = "FLT4S")
    
    tabla
  }
  
  tabla_pfi <- bind_rows(Map(calcular_pfi_raster, rasteres,
                             tabla_rasteres$especie,
                             tabla_rasteres$anio,
                             tabla_rasteres$nivel_habitat))
  
  write_csv(tabla_pfi, file.path(dir_pfi, "OE3_PFI_por_parche.csv"))
}
