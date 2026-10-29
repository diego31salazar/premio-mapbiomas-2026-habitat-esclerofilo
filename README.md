# Estructura espacial del hábitat potencial de especies del bosque y matorral esclerófilo de Chile central

Código de la tesis de Magíster en Geografía (Universidad de Chile, 2026):

> **Análisis de la evolución de la estructura espacial de hábitats para las especies *Jubaea chilensis*, *Quillaja saponaria* y *Vachellia caven* bajo un enfoque que integra los cambios de uso y cobertura de suelo y la idoneidad ambiental, en las regiones de Valparaíso, Metropolitana de Santiago y del Libertador Bernardo O'Higgins.**
>
> Autor: Diego Salazar Canales · Profesor guía: Pablo Sarricolea Espinoza

Trabajo postulado al **Premio MapBiomas Chile, 1ª edición (2026)**.

---

## Resumen del enfoque

El estudio integra tres productos espaciales para evaluar cómo cambia la estructura del hábitat potencial de tres especies leñosas entre 2002 y 2042:

| Objetivo | Pregunta | Método principal |
|---|---|---|
| **OE1** | ¿Cómo cambian las coberturas de suelo y cómo podrían configurarse al 2042? | Cadenas de Márkov + Perceptrón Multicapa (Land Change Modeler, TerrSet) sobre MapBiomas Chile Colección 1 |
| **OE2** | ¿Dónde se encuentran las condiciones ambientales idóneas para cada especie? | Ensamble de modelos de distribución de especies (GLM, Random Forest y Maxent) con `biomod2` |
| **OE3** | ¿Cómo evoluciona la cantidad y fragmentación del hábitat potencial? | Métricas de paisaje (`landscapemetrics`), Mínimos Cuadrados Generalizados con AR1 (`nlme`) e índice de fragmentación por parche adaptado de Rivas et al. (2022) |

El hábitat potencial se define como el cruce entre la idoneidad ambiental (OE2) y las coberturas de bosque y matorral de cada año (OE1), en tres niveles anidados: NH1 (0,8–1), NH2 (0,6–1) y NH3 (0,4–1).

---

## Contenido del repositorio

```
├── OE1/OE1_coberturas_drivers_validacion.R
├── OE2/OE2_distribucion_especies.R
├── OE3/OE3_estructura_espacial_habitat.R
├── README.md
└── LICENSE
```

| Script | Apartado de la metodología | Qué hace |
|---|---|---|
| `OE1_coberturas_drivers_validacion.R` | 4.2.1.1 y 4.2.1.3 | Normalización de *drivers* por percentil 95, alineación a la grilla de referencia, correlación de Pearson entre *drivers* (Anexo 1) y cálculo de Hits, Misses, False Alarms y Figure of Merit de la validación 2022 (Anexo 2) |
| `OE2_distribucion_especies.R` | 4.2.2 | Depuración de ocurrencias GBIF y *spatial thinning*, remuestreo de predictores a 250 m, correlaciones en dos etapas, modelamiento con `biomod2`, ensamble EMwmean, evaluación (AUC, TSS) y reclasificación de la idoneidad |
| `OE3_estructura_espacial_habitat.R` | 4.2.3.2 y 4.2.3.3 | Cantidad de hábitat (CA) y número de parches (NP), modelos de Mínimos Cuadrados Generalizados con estructura AR1 (Anexos 5 y 6) e índice de fragmentación adaptado (PFI) |

Cada script está organizado en bloques numerados que siguen el orden de la metodología. Los bloques de mayor costo computacional (preprocesamiento de predictores, modelamiento con `biomod2` y cálculo del PFI) se activan mediante interruptores al inicio de cada script.

---

## Pasos realizados fuera de R

| Paso | Software | Detalle |
|---|---|---|
| Reclasificación de MapBiomas Chile (2002, 2012, 2022) en seis clases | ArcGIS Pro 3.7.1 | Bosque, matorral, plantación forestal, agricultura y ganadería, urbanización y otros |
| *Drivers* de proximidad | ArcGIS Pro 3.7.1 | Distancia euclidiana a 30 m; normalización con Raster Calculator (expresiones equivalentes en el bloque 2 del script del OE1) |
| Matrices de Márkov, submodelos de transición, entrenamiento MLP, simulación 2022 y proyección 2042 | TerrSet, Land Change Modeler | Submodelos urbano, agrícola-ganadero y forestal entrenados con el periodo 2002–2012 |
| Kappa global y por clase | TerrSet (VALIDATE / CROSSTAB) | Comparación entre el mapa observado y el simulado de 2022 |
| Descarga de SoilGrids | Google Earth Engine | Ocho variables, promedio de seis intervalos de profundidad (0–200 cm), 250 m |
| Delimitación del hábitat potencial | ArcGIS Pro 3.7.1 | Niveles NH1–NH3; máscaras de bosque y matorral por año; celdas de 250 m con intersección ≥ 62.499 m²; Extract by Mask |

---

## Datos de entrada

Los datos no se incluyen en el repositorio. Todos provienen de fuentes públicas y deben citarse según los términos de uso de cada proveedor:

| Datos | Fuente |
|---|---|
| Coberturas de suelo 2002, 2012 y 2022 (30 m) | Proyecto MapBiomas Chile – Colección 1.0 de la Serie Anual de Mapas de Cobertura y Uso del Suelo de Chile. Consultada a inicios de 2025 mediante enlaces de descarga directa por año, vigentes en esa fecha; acceso actual a través de [chile.mapbiomas.org](https://chile.mapbiomas.org) y Google Earth Engine. Datos bajo licencia CC-BY. |
| Registros de ocurrencia | [GBIF Chile](https://gbifchile.mma.gob.cl) |
| Variables bioclimáticas (30 s, ~1 km) | [WorldClim 2.1](https://www.worldclim.org) |
| Variables edáficas (250 m) | [SoilGrids](https://soilgrids.org), ISRIC |
| Topografía y heterogeneidad de hábitat (1 km) | [EarthEnv](https://www.earthenv.org) |
| Áreas pobladas | Biblioteca del Congreso Nacional, mapas vectoriales |

---

**Cita de MapBiomas Chile:**

> Proyecto MapBiomas Chile. (s. f.). *Colección 1.0 de la Serie Anual de Mapas de Cobertura y Uso del Suelo de Chile* [Conjunto de datos]. Consultado a inicios de 2025 mediante enlaces de descarga directa por año, vigentes en esa fecha. Acceso actual a través de https://chile.mapbiomas.org.

---

## Software y versiones

- R 4.6.1
- `terra` 1.9-50 · `biomod2` 4.3-4-5 · `landscapemetrics` · `nlme` · `sf` · `caret` · `dplyr` · `tidyr` · `ggplot2`
- ArcGIS Pro 3.7.1
- TerrSet (Land Change Modeler)
- Google Earth Engine

---

## Cómo citar

Salazar Canales, D. (2026). *Análisis de la evolución de la estructura espacial de hábitats para las especies Jubaea chilensis, Quillaja saponaria y Vachellia caven bajo un enfoque que integra los cambios de uso y cobertura de suelo y la idoneidad ambiental* [Tesis de Magíster en Geografía, Universidad de Chile]. Código disponible en: https://github.com/diego31salazar/premio-mapbiomas-2026-habitat-esclerofilo

---

## Nota

Los scripts consolidan y documentan el procesamiento realizado durante la tesis. Para su orden y documentación se contó con apoyo de herramientas de inteligencia artificial; las decisiones metodológicas, el análisis y la interpretación de los resultados son del autor.

## Licencia

Código distribuido bajo licencia [MIT](LICENSE).
