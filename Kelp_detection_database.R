# Macrocystis pyrifera database from remote sensing imaginery
# layers download using Mora et al. 2020. Remote sensing Kelp Difference (KD) index and Floating Algae Index (FAI)

##### libraries ####################### ####
library(tidyverse)
library(purrr)
library(raster)
library(terra)
library(sf)
library(furrr)
setwd("C:/Users/Acer/OneDrive - uc.cl/Documentos/GitHub/Macrocystis_pyrifera/")
#######################################

##### ecoregional polygons ############
#In this section, the shape files that will be used for remote sensing through Sentinel-2 images are generated.

#setting the spatial domain
marspec <- load_layers(c("MS_biogeo05_dist_shore_5m"), datadir = "./Capas_climaticas/FRE_selected/MARSPEC")
ext <- extent(-76, -66, -41.84, -18) #Study area
marspec <- crop(marspec, ext)
dist_pol <- rasterToPolygons(marspec$MS_biogeo05_dist_shore_5m, fun=function(x){x <= 8}, dissolve = T)
dist_pol <- buffer(dist_pol, 2000)
plot(dist_pol)
dist_pol <- aggregate(dist_pol, fact= 1, dissolve= T)
dist_pol <- st_as_sf(dist_pol)
plot(dist_pol)
eco_zee <- st_intersection(ecoregions, dist_pol)
#plot(eco_zee)
eco_zee$ECOREGION <- factor(eco_zee$ECOREGION, levels = c("Humboldtian","Central Chile",  "Araucanian"))


# Dividir los polígonos

eco_zee_Hum <- filter(eco_zee, ECOREGION== "Humboldtian")
eco_zee_CC <- filter(eco_zee, ECOREGION== "Central Chile")
eco_zee_Arc <- filter(eco_zee, ECOREGION== "Araucanian")
# Crear línea horizontal en -21.5° S (como LINESTRING)
bbox_Hum <- st_bbox(eco_zee_Hum)
bbox_CC <- st_bbox(eco_zee_CC)
bbox_Arc <- st_bbox(eco_zee_Arc)

linea_corte_Hum <- st_linestring(
  rbind(
    c(bbox_Hum["xmin"], -21.5),  # Punto izquierdo (oeste)
    c(bbox_Hum["xmax"], -21.5)   # Punto derecho (este)
  )
) %>% 
  st_sfc(crs = st_crs(eco_zee))  # Mismo CRS que el original
linea_corte_CC <- st_linestring(
  rbind(
    c(bbox_CC["xmin"], -29),  # Punto izquierdo (oeste)
    c(bbox_CC["xmax"], -29)   # Punto derecho (este)
  )
) %>% 
  st_sfc(crs = st_crs(eco_zee))  # Mismo CRS que el original
linea_corte_Arc <- st_linestring(
  rbind(
    c(bbox_Arc["xmin"], -37.5),  # Punto izquierdo (oeste)
    c(bbox_Arc["xmax"], -37.5)   # Punto derecho (este)
  )
) %>% 
  st_sfc(crs = st_crs(eco_zee))  # Mismo CRS que el original

# Dividir el polígono con la línea
eco_zee_Hum_dividido <- eco_zee_Hum %>% 
  lwgeom::st_split(linea_corte_Hum) %>%  # Aplica el corte
  st_collection_extract("POLYGON") %>%  # Extrae solo polígonos
  st_as_sf() %>%  # Convierte a sf
  mutate(
    # Calcular latitud (Y) del centroide de cada subpolígono
    centroide_y = st_coordinates(st_centroid(.))[, 2],  
    # Asignar etiqueta según posición relativa a -21.5° S
    subregion = ifelse(centroide_y > -21.5, "Humboldtian_N", "Humboldtian_S"),
    .before = 1  # Opcional: mover la nueva columna al inicio
  ) %>% 
  dplyr::select(-centroide_y) %>% 
  dplyr::select(ECOREGION, subregion)# Eliminar columna temporal (opcional)
eco_zee_CC_dividido <- eco_zee_CC %>% 
  lwgeom::st_split(linea_corte_CC) %>%  # Aplica el corte
  st_collection_extract("POLYGON") %>%  # Extrae solo polígonos
  st_as_sf() %>%  # Convierte a sf
  mutate(
    # Calcular latitud (Y) del centroide de cada subpolígono
    centroide_y = st_coordinates(st_centroid(.))[, 2],  
    # Asignar etiqueta según posición relativa a -21.5° S
    subregion = ifelse(centroide_y > -29, "Central_Chile_N", "Central_Chile_S"),
    .before = 1  # Opcional: mover la nueva columna al inicio
  ) %>% 
  dplyr::select(-centroide_y) %>% 
  dplyr::select(ECOREGION, subregion) %>% 
  group_by(ECOREGION, subregion) %>%  # Agrupar por estas columnas
  summarise(geometry = st_union(geometry)) %>%  # Unir geometrías
  ungroup()  # Opcional: remover agrupamiento
eco_zee_Arc_dividido <- eco_zee_Arc %>% 
  lwgeom::st_split(linea_corte_Arc) %>%  # Aplica el corte
  st_collection_extract("POLYGON") %>%  # Extrae solo polígonos
  st_as_sf() %>%  # Convierte a sf
  mutate(
    # Calcular latitud (Y) del centroide de cada subpolígono
    centroide_y = st_coordinates(st_centroid(.))[, 2],  
    # Asignar etiqueta según posición relativa a -21.5° S
    subregion = ifelse(centroide_y > -37.5, "Araucanian_N", "Araucanian_S"),
    .before = 1  # Opcional: mover la nueva columna al inicio
  ) %>% 
  dplyr::select(-centroide_y) %>% 
  dplyr::select(ECOREGION, subregion) %>% 
  group_by(ECOREGION, subregion) %>%  # Agrupar por estas columnas
  summarise(geometry = st_union(geometry)) %>%  # Unir geometrías
  ungroup()  # Opcional: remover agrupamiento

ggplot()+geom_sf(data=eco_zee_Arc_dividido, aes(fill= subregion))+geom_sf(data= linea_corte_Arc,aes(geometry= geometry))
GEE_geom <- bind_rows(eco_zee_Arc_dividido, eco_zee_CC_dividido, eco_zee_Hum_dividido)
GEE_geom$subregion <- as.factor(GEE_geom$subregion)
# Guardar cada ecorregión como un shapefile separado
for (ecorregion in levels(GEE_geom$subregion)) {
  # Filtrar la ecorregión actual
  temp <- GEE_geom %>% filter(subregion == ecorregion)
  
  # Crear un nombre de archivo válido (sin espacios)
  filename <- gsub(" ", "_", tolower(ecorregion))
  filepath <- paste0("./ecorregiones_shp/", filename, ".shp")
  
  # Guardar el shapefile
  st_write(temp, filepath, delete_dsn = TRUE)
  
  # Mensaje de confirmación
  message(paste("Shapefile guardado:", filepath))
}

# Calcular bbox para cada ecorregión
bbox_list <- lapply(split(eco_zee, eco_zee$ECOREGION), st_bbox)
bbox_list
# Definir los límites para cada ecorregión (extraídos de bbox_list)
xlims <- list(
  "Humboldtian" = c(-71.00, -70.00),
  "Central Chile" = c(-71.83, -70.42),
  "Araucanian" = c(-74.17, -71.58)
)

ylims <- list(
  "Humboldtian" = c(-25.12, -18.00),
  "Central Chile" = c(-33.13, -25.11),
  "Araucanian" = c(-41.83, -33.12)
)

library(cowplot)

plot_list <- lapply(levels(eco_zee$ECOREGION), function(ecor) {
  datos <- subset(eco_zee, ECOREGION == ecor)
  bbox <- st_bbox(datos)
  
  ggplot() +
    geom_sf(data = datos, aes(fill = ECOREGION), show.legend = F) +
    geom_sf(data = chile, aes(geometry = geometry)) +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"]),
      expand = FALSE
    ) +
    ggtitle(ecor) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0, size= 6))
})

# Combinar los gráficos
plot_grid(plotlist = plot_list, nrow = 1)

## The generated shape files are loaded into GEE to obtain images of the 
#area of interest and process the information to discretize the floating canopy signal.

############################################################

#### Stack temporal series and enframe it ################ ####
# The raster images obtained through remote recognition of the 
# floating canopy are read and processed to obtain the area of 
# each floating canopy patch in space and time.

# 1. Función mejorada para cargar rasters
load_rasters <- function(path) {
  # Listar archivos .tif
  tif_files <- list.files(path, pattern = "\\.tif$", full.names = TRUE)
  
  # Verificar si hay archivos
  if (length(tif_files) == 0) {
    warning(paste("No se encontraron archivos .tif en:", path))
    return(NULL)
  }
  
  # Leer rasters y convertirlos a lista nombrada
  raster_list <- map(tif_files, ~raster(.x)) %>% 
    set_names(tools::file_path_sans_ext(basename(tif_files)))
  
  return(raster_list)
}

# 2. Función para procesar todos los directorios
process_all_directories <- function(base_dir, dir_pattern = "Ecoregions") {
  # Encontrar todos los directorios relevantes
  dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE) %>% 
    str_subset(dir_pattern)
  
  # Cargar todos los rasters en una lista nombrada
  all_rasters <- map(dirs, load_rasters) %>% 
    set_names(basename(dirs))
  
  # Eliminar elementos nulos (directorios vacíos)
  all_rasters <- compact(all_rasters)
  
  return(all_rasters)
}

# 3. Función para convertir raster a sf y calcular áreas
raster_to_sf_area <- function(r) {
  # Convertir a polígonos (solo valores = 1)
  polys <- rasterToPolygons(r, fun = function(x) x == 1, dissolve = F)
  
  # Si no hay polígonos, retornar NULL
  if (is.null(polys) || nrow(polys) == 0) {
    return(NULL)
  }
  
  # Convertir a sf y procesar
  sf_obj <- st_as_sf(polys) %>% 
    st_transform(crs = "EPSG:32719") %>% 
    #st_buffer(20) %>%        # Buffer de 20m
    st_union() %>%           # Unir polígonos
    st_cast("POLYGON") %>%   # Separar multipolígonos
    st_as_sf() %>%           # Convertir a sf
    mutate(area = st_area(.)) # Calcular área
  
  return(sf_obj)
}


# 4. Función para extraer metadatos de los nombres 
extract_metadata <- function(label) {
  # Extraer componentes del nombre del archivo
  parts <- str_split(label, "_")[[1]]
  
  # Formato: lat1_lat2_season_year
    meta <- tibble(
      Ecoregion = as.character(parts[1]),
      subdivision = as.character(parts[2]),
      season = parts[3],  # Mantener como character por ahora
      year = as.numeric(parts[4])
    )
  
  # Normalizar nombres de temporadas y agregar trimestre
  meta <- meta %>%
    mutate(
      season = case_when(
        str_starts(season, "Sum") ~ "Sum",
        str_starts(season, "Fall|Aut") ~ "Fall",
        str_starts(season, "Win") ~ "Win",
        str_starts(season, "Spr") ~ "Spr",
        TRUE ~ season
      ),
      season = factor(season, levels = c("Sum", "Fall", "Win", "Spr")),
      quarter = case_when(
        season == "Sum" ~ "Q1",
        season == "Fall" ~ "Q2",
        season == "Win" ~ "Q3",
        season == "Spr" ~ "Q4"
      )
    )
  
  return(meta)
}

# Directorio base
base_dir <- "./GitHub/Macrocystis_pyrifera/Remote_sensing_database/"  #Set according to where the data is hosted

# Cargar todos los rasters
all_rasters <- process_all_directories(base_dir)
# Convertir a tibble anidada
raster_nested <- imap_dfr(all_rasters, function(raster_list, dir_name) {
  imap_dfr(raster_list, function(r, label) {
    # Extraer metadatos
    meta <- extract_metadata(label)
    
    # Crear tibble con metadatos y raster
    tibble(
      directory = dir_name,
      label = label,
      raster = list(r),
      season = meta$season,
      year= meta$year,
      Ecoregion = meta$Ecoregion,
      subdivision = meta$subdivision,
      Eco_sub= paste0(meta$Ecoregion,"_",meta$subdivision)
    )
  })
})
#check if there is any raster layer missing
(missing.layers <- complete(raster_nested, season, year, Eco_sub) %>% filter(is.na(label)))


# Procesar rasters para obtener polígonos y áreas (puede ser lento)
start.time <- Sys.time()
raster_processed <- raster_nested %>% filter(Ecoregion== "Hum") %>% 
  mutate(
    polygons = future_map(raster, raster_to_sf_area, .progress = TRUE) # Usar future_map para paralelizar
  )
end.time <- Sys.time()
print(round(end.time - start.time,2))
#spatial feature with sst information for the cluster between 18°S and 20°S in each season and 
sf_unnest <- raster_processed %>% dplyr::select(season, year, polygons) %>% unnest(polygons)
saveRDS(sf_unnest, "MP_st_Hum")

hum <- readRDS("MP_st_Hum")
cc <- readRDS("MP_st_CC")
arc <- readRDS("MP_st_Arc")
data <- bind_rows(hum, cc, arc)
saveRDS(data, "Mp_HCS_st")
############################################################################

### filtered data ##### 
#Data filtering through direct scanning of signals obtained from satellite images

data <- readRDS("./Mp_HCS_st")
#(missing.layers <- complete(data, season, year) %>% filter(is.na(area)))
data$area <- as.numeric(data$area)
data <- filter(data, area> 1000)
data <- rename(data, geometry= x)

sf <- st_as_sf(data, crs= "EPSG:32719" ) # 'EPSG:4326'
sf <- st_intersection(sf, ecozee_utm[,"ECOREGION"])
sf <- sf  %>% mutate(quarter= case_when(season %in% "Sum" ~ "Q1",
                                        season %in% "Fall" ~ "Q2",
                                        season %in% "Win" ~ "Q3",
                                        season %in% "Spr" ~ "Q4")) %>% 
  mutate(yearqtr= paste(year, quarter, sep= " ")) %>% 
  mutate(centroid= st_centroid(geometry)) %>% 
  mutate(yearqtr= zoo::as.yearqtr(yearqtr, format = "%Y Q%q"))


sf$lat <- st_coordinates(st_transform(sf$centroid, crs = "EPSG:4326"))[,"Y"]
sf$lon <- st_coordinates(st_transform(sf$centroid, crs = "EPSG:4326"))[,"X"]
sf$X <- st_coordinates(sf$centroid)[,"X"]/1000
sf$Y <- st_coordinates(sf$centroid)[,"Y"]/1000

sf <- st_transform(sf, crs= "EPSG:4326")
df <- st_drop_geometry(sf)


sf2 <- sf %>%
  filter(!(lat < -30.70 & lat > -30.75 & lon > -71.7)) %>% # aguas interiores de rio Limarí
  filter(!(lat < -33.61 & lat > -33.65 & lon > -71.7)) %>% # aguas interiores de rio Maipo
  filter(!(lat < -33.904 & lat > -33.95 & lon > -71.85)) %>% # aguas interiores de rio Rapel
  filter(!(lat < -34.47 & lat > -34.48 & lon > -72.15)) %>% # aguas interiores Nilahue
  filter(!(lat < -34.63 & lat > -34.65 & lon > -72.1)) %>% # aguas interiores de Humedal Bucalemu
  filter(!(lat < -34.7 & lat > -34.72 & lon > -72.08)) %>% #Desembocadura Rio Boyeruca
  filter(!(lat < -34.75 & lat > -34.76 & lon > -72.14)) %>% #Desembocadura Estero Llico
  filter(!(lat < -34.82 & lat > -34.84 & lon > -72.1)) %>% #Laguna Petrel
  filter(!(lat < -34.95 & lat > -35.02 & lon > -72.2)) %>% #Desembocadura Rio Mataquito
  filter(!(lat < -35.32 & lat > -35.33 & lon > -72.42)) %>% #Celulosa Arauco 
  filter(!(lat < -36.38 & lat > -36.40 & lon > -72.89)) %>% #Desembocadura Rio Itata 
  filter(!(lat < -36.74 & lat > -36.75 & lon > -73.105 & lon < -73)) %>% #Humedal Rocuant Andalien 
  filter(!(lat < -37.745 & lat > -38.24 & lon > -73.66)) %>% #Estero Pangue, Lebu , Cañete, rio Lleulleu, Rio Quidico
  filter(!(lat < -38.9 & lat > -39.295 & lon > -73.49)) %>% #Desembocadura rio Tolten, sur Lago Budi
  filter(!(lat < -39.43 & lat > -39.46 & lon > -73.23)) %>% # Desembocadura rio Mehuin
  filter(!(lat < -39.65 & lat > -39.67 & lon > -73.3)) %>% # Registro en Rio Cruces
  filter(!(lat < -39.8 & lat > -39.9 & lon > -73.385)) %>% # aguas interiores de rio Valdivia
  filter(!(lat < -39.9 & lat > -39.96 & lon > -73.4)) %>% # aguas interiores de rio Valdivia
  filter(!(lat < -40 & lat > -40.02 & lon > -73.7)) %>% # aguas interiores de rio Valdivia
  filter(!(lat < -36.80 & lat > -36.82 & lon > -73.17)) %>% # aguas interiores de rio Bio-Bio
  filter(!(lat < -37.22 & lat > -37.25 & lon > -73.45)) %>% # Humedal Tubul Raqui
  filter(!(lat < -37.60 & lat > -37.625 & lon > -73.65)) %>% # Rio Lebu
  filter(!(lat < -39.36 & lat > -39.395 & lon > -73.23)) %>% # Rio Quele
  filter(!(lat < -39.43 & lat > -39.45 & lon > -73.15)) %>% # Rio Lingue
  filter(!(lat < -40.68 & lat > -40.70 & lon > -73.75)) %>% # Rio Huellehue
  filter(!(lat < -40.77 & lat > -40.775 & lon > -73.84)) %>% # Rio Cholguaco
  filter(!(lat < -41.21 & lat > -41.22 & lon > -73.9)) %>% # Rio San Juanito
  filter(!(lat < -41.27 & lat > -41.29 & lon > -73.9)) %>% # Rio Llico
  filter(!(lat < -41.485 & lat > -41.51 & lon > -73.8)) %>% # Rio Quenuir
  filter(!(lat < -40.08 & lat > -40.09 & lon > -73.66)) %>% # Rio Colium
  filter(!(lat < -40.146 & lat > -40.155 & lon > -73.67)) %>% # Estero Hueicolla
  filter(!(lat < -40.235 & lat > -40.255 & lon > -73.73)) %>% # Rio Bueno
  filter(!(lat < -40.53 & lat > -40.54 & lon > -73.8)) %>% # Rio Bueno
  filter(!(yearqtr=="2017 Q2" & lat < -39.8 & lat > -39.94 & lon > -73.6 & lon < -73.25 )) %>% # missrecognition cloud interference
  filter(!(yearqtr=="2017 Q3" & lat < -40.8 & lat > -40.94 & lon > -73.88 & lon < -73.8 )) %>%  # missrecognition cloud interference
  filter(!(lat < -18 & lat > -18.26)) %>% #Offshore signal 
  filter(!(lat < -18.42 & lat > -18.47 & yearqtr == "2021 Q3")) %>% #Offshore signal 
  filter(!(lat < -18.5 & lat > -18.7 & lon < -70.37)) %>% #Offshore signal 
  filter(!(lat < -19.31 & lat > -19.32 & lon < -70.2)) %>% #Offshore signal 
  filter(!(lat < -19.66 & lat > -19.82 & lon < -70.18 & yearqtr== "2021 Q3")) %>% #Offshore signal 
  filter(!(lat < -20 & lat > -20.1 & lon < -70.175)) %>% #Offshore signal 
  filter(!(lat < -20.35 & lat > -20.376 & lon < -70.21)) %>% #Offshore signal 
  filter(!(lat < -22.23 & lat > -22.28 & lon < -70.26)) %>% #Offshore signal 
  filter(!(lat < -22.32 & lat > -22.33 & lon < -70.3)) %>% #Offshore signal 
  filter(!(lat < -22.37 & lat > -22.41 & lon < -70.2 & yearqtr== "2016 Q1")) %>% #Offshore signal 
  filter(!(lat < -22.42 & lat > -22.45 & lon < -70.261)) %>% #Offshore signal
  filter(!(lat < -22.51 & lat > -22.53 & yearqtr== "2022 Q3")) %>% #Cloude interference signal
  filter(!(lat < -22.56 & lat > -23.19 & yearqtr== "2018 Q1")) %>% #Cloude interference signal
  filter(!(lat < -22.56 & lat > -22.62 & yearqtr== "2025 Q1")) %>% #Cloude interference signal
  filter(!(lat < -23.18 & lat > -23.19 & lon < -70.6)) %>% #Cloude interference signal
  filter(!(lat < -23.3 & lat > -23.32 & lon < -70.6)) %>% #Cloude interference signal
  filter(!(lat < -23.95 & lat > -24.01 & lon < -70.5)) %>% #Cloude interference signal
  filter(!(lat < -24.54 & lat > -24.6)) %>% #Cloude interference signal
  filter(!(lat < -24.39 & lat > -25.16 & lon < -70.51)) %>% #Cloude interference signal
  filter(!(lat < -25 & lat > -25.01 & lon < -70.48)) %>% #Cloude interference signal
  filter(!(lat < -25.12 & lat > -25.5 & lon < -70.47)) %>% #Cloude interference signal
  filter(!(lat < -25.91 & lat > -26.26 & lon < -70.65)) %>% #Cloude interference signal
  filter(!(lat < -26.75 & lat > -26.765 & lon < -70.7)) %>% #offshore signal
  filter(!(lat < -26.83 & lat > -26.96 & lon < -70.86)) %>% #offshore signal
  filter(!(lat < -27 & lat > -27.36 & lon < -70.98)) %>% #Cloude interference signal
  filter(!(lat < -27.36 & lat > -27.57 & lon < -70.92)) %>% #Cloude interference signal
  filter(!(lat < -27.6 & lat > -27.85 & lon < -71.1)) %>% #Cloude interference signal
  filter(!(lat < -27.85 & lat > -28.05 & lon < -71.17)) %>% #Cloude interference signal
  filter(!(lat < -28.05 & lat > -28.21 & lon < -71.2)) %>% #Cloude interference signal
  filter(!(lat < -28.10 & lat > -28.25 & yearqtr == "2016 Q3")) %>% #off shore signal
  filter(!(lat < -28.48 & lat > -28.68 & lon < -71.3)) %>% #Cloude interference signal
  filter(!(lat < -28.73 & lat > -28.90 & lon < -71.3)) %>% #Cloude interference signal
  filter(!(lat < -28.9 & lat > -29.01 & lon < -71.52)) %>% #Cloude interference signal
  filter(!(lat < -29.25 & lat > -29.26 & lon > -71.5 & yearqtr== "2024 Q2")) %>% #Missrecognition
  filter(!(lat < -29 & lat > -29.21 & lon < -71.5)) %>% #Cloude interference signal/offshore signal
  filter(!(lat < -29.6 & lat > -29.9 & lon < -71.37)) %>% #Cloude interference signal/offshore signal
  filter(!(lat < -30.1 & lat > -30.25 & lon < -71.54)) %>% #offshore signal
  filter(!(lat < -30.6 & lat > -30.65 & lon < -71.54)) %>% #offshore signal
  filter(!(lat < -31.05 & lat > -31.1 & lon < -71.7)) %>% #offshore signal
  filter(!(lat < -31.17 & lat > -31.18 & lon < -71.6)) %>% #offshore signal
  filter(!(lat < -31.34 & lat > -31.38 & lon < -71.6)) %>% #offshore signal
  filter(!(lat < -31.4 & lat > -31.53 & lon < -71.6)) %>% #offshore signal
  filter(!(lat < -31.5 & lat > -31.61 & lon < -71.55)) %>% #offshore signal
  filter(!(lat < -31.6 & lat > -31.7 & lon < -71.55)) %>% #offshore signal
  filter(!(lat < -31.7 & lat > -31.72 & lon < -71.55)) %>% #offshore signal
  filter(!(lat < -31.8 & lat > -32.1 & lon < -71.53)) %>% #offshore signal
  filter(!(lat < -32.4 & lat > -32.5 & lon < -71.45)) %>% #offshore signal
  filter(!(lat < -32.52 & lat > -32.60 & lon < -71.45)) %>% #offshore signal 
  filter(!(lat < -32.66 & lat > -32.768 & lon < -71.5)) %>% #offshore signal  
  filter(!(lat < -32.9 & lat > -32.95 & lon < -71.5)) %>% #offshore signal
  filter(!(lat < -33.06 & lat > -33.085 & lon < -71.5)) %>% #offshore signal 
  filter(!(lat < -33.7 & lat > -33.75 & lon < -71.5)) %>% #offshore signal 
  filter(!(lat < -33.83 & lat > -33.85 & lon < -71.5)) %>% #offshore signal 
  filter(!(lat < -34.67 & lat > -34.68 & lon < -72)) %>% #offshore signal 
  filter(!(lat < -34.75 & lat > -34.765 & lon < -72.14)) %>% #offshore signal 
  filter(!(lat < -34.91 & lat > -34.92 & lon < -72.14)) %>% #offshore signal 
  filter(!(lat < -35.31 & lat > -35.54 & lon < -72.5)) %>% #offshore signal
  filter(!(lat < -35.4 & lat > -35.41 & lon < -72.5)) %>% #offshore signal
  filter(!(lat < -35.54 & lat > -35.56 & lon < -72.5)) %>% #offshore signal
  filter(!(lat < -35.63 & lat > -35.7 & lon < -72.5)) %>% #offshore signal
  filter(!(lat < -35.925 & lat > -35.94 & lon < -72.73)) %>% #offshore signal
  filter(!(lat < -36.16 & lat > -36.18 & lon < -72.782)) %>% #offshore signal
  filter(!(lat < -36.34 & lat > -36.35 & lon < -72.782)) %>% #offshore signal
  filter(!(lat < -37.24 & lat > -37.25 & lon < -73.5)) %>% #offshore signal
  filter(!(lat < -38.6 & lat > -38.95)) %>%  # interference with stuarine systems where spectral signal cluld be related to wetland vegetations
  filter(!(lat < -39.76 & lat > -39.79 & lon < -73.4)) %>%  # interference with stuarine systems where spectral signal cluld be related to wetland vegetations
  filter(!(lat < -39.97 & lat > -39.98 & lon < -73.7)) %>% # interference with stuarine systems where spectral signal cluld be related to wetland vegetations
  filter(!(lat < -39.69 & lat > -39.67 & lon < -73.41)) %>% #offshore signal
  filter(!(lat < -39.935 & lat > -39.96 & yearqtr == "2021 Q3")) %>%  # Cloud interference
  filter(!(lat < -39.97 & lat > -39.995 & lon < -73.6)) %>%  #offshore signal/Cloud interference
  filter(!(lat < -40.15 & lat > -40.16 & yearqtr == "2022 Q3" & lon > -73.675)) %>% # Cloud interference
  filter(!(lat < -40.25 & lat > -40.43 & lon < -73.76)) %>%  #offshore signal/Cloud interference
  filter(!(lat < -40.36 & lat > -40.415 & lon < -73.77)) %>%  #offshore signal/Cloud interference
  filter(!(lat < -40.45 & lat > -40.47 & lon < -73.78)) %>%  #offshore signal/Cloud interference
  filter(!(lat < -40.57 & lat > -41.03 & lon < -73.75 & yearqtr== "2020 Q3")) %>%  #offshore signal/Cloud interference
  filter(!(lat < -40.63 & lat > -40.73 & lon < -73.8)) %>% #offshore signal/Cloud interference
  filter(!(lat < -40.87 & lat > -40.88 & lon < -73.88)) %>% #offshore signal/Cloud interference
  filter(!(lat < -40.89 & lat > -40.93 & lon < -73.86)) %>% #offshore signal/Cloud interference
  filter(!(lat < -41.22 & lat > -41.48 & lon < -73.75 & yearqtr== "2020 Q3")) %>% #offshore signal/Cloud interference
  filter(!(lat < -41.21 & lat > -41.29 & lon < -73.75 & yearqtr == "2021 Q3")) %>% #offshore signal/Cloud interference
  filter(!(lat < -41.43 & lat > -41.46 & lon < -73.75 & yearqtr== "2022 Q3")) %>% #offshore signal/Cloud interference
  filter(!(lat < -41.47 & lat > -41.48 & lon < -73.8)) %>% #offshore signal/Cloud interference
  filter(!(lat < -41.48 & lat > -41.52 & lon < -73.79 & yearqtr== "2018 Q2")) #offshore signal/Cloud interference


loc <- sf2 %>% st_drop_geometry() %>%
  mutate(lat2= round(lat, 2)) %>%
  dplyr::select(yearqtr, lat, lat2, lon, area) %>%
  group_by(lat2) %>%
  filter(area== max(area)) %>%
  ungroup()


#### verificacion a posteriori 
# Extraer los valores únicos de yearqtr y ordenarlos
yearqtr_unique <- as.character(sort(unique(df$yearqtr)))
#yearqtr_unique <- as.character(sort(yearqtr_unique))  # Ordenar los trimestres

# Crear un mapeo entre índices y etiquetas yearqtr
yearqtr_mapping <- setNames(seq_along(yearqtr_unique), yearqtr_unique)

# Convertir loc a objeto sf
puntos_sf <- st_as_sf(loc, coords = c("lon", "lat"), crs = "EPSG:4326")
puntos_sf$lat <- st_coordinates(puntos_sf)[,"Y"]
puntos_sf$lon <- st_coordinates(puntos_sf)[,"X"]
puntos_sf$id <- seq(1, nrow(puntos_sf), by = 1)  # Asignar un ID único a cada punto

# Interfaz de usuario (UI)
ui <- fluidPage(
  titlePanel("Mapa Interactivo con Línea Temporal (Trimestres)"),
  sidebarLayout(
    sidebarPanel(
      sliderTextInput(  # Usar sliderTextInput en lugar de sliderInput
        inputId = "yearqtr_index",  # ID del control deslizante
        label = "Selecciona un trimestre:",  # Etiqueta
        choices = yearqtr_unique,  # Etiquetas personalizadas
        selected = yearqtr_unique[1],  # Valor inicial (primer trimestre)
        animate = FALSE  # Animación automática
      )
    ),
    mainPanel(
      leafletOutput("mapa")  # Salida del mapa
    )
  )
)


# Lógica del servidor (Server)
server <- function(input, output, session) {
  # Variable reactiva para almacenar el punto de interés
  punto_interes_sf <- reactiveVal(NULL)
  
  # Buffer alrededor del punto de interés
  buffer_zoom <- reactive({
    req(punto_interes_sf())
    st_buffer(punto_interes_sf(), dist = 0.05)  # Crear buffer alrededor del punto
  })
  
  # Observar clics en los puntos para definir el punto de interés
  observeEvent(input$mapa_marker_click, {
    click <- input$mapa_marker_click
    punto <- st_point(c(click$lng, click$lat))
    punto_interes_sf(st_sfc(punto, crs = "EPSG:4326"))
  })
  
  # Obtener la etiqueta yearqtr correspondiente al valor seleccionado
  yearqtr_label <- reactive({
    input$yearqtr_index  # Usar el valor seleccionado en el sliderTextInput
  })
  
  # Filtrar datos según el trimestre seleccionado
  datos_filtrados <- reactive({
    req(buffer_zoom())
    
    #sf2_z <-  st_transform(sf2, crs = 4326)  # Asegurar que esté en WGS84
    sf2_zoom <- st_intersection(sf2, buffer_zoom())  # Recortar sf2 con el buffer
    sf2_zoom %>% filter(yearqtr == yearqtr_label())
  })
  
  # Renderizar el mapa inicial
  output$mapa <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%  # Capa base del mapa
      addScaleBar(position = "bottomleft", options = scaleBarOptions(metric = TRUE)) %>%  # Agregar escala en km
      addCircleMarkers(
        data = puntos_sf,
        radius = 5,
        color = "blue",
        fillOpacity = 1,
        layerId = ~id,  # Asignar un ID único a cada marcador
        popup = ~paste("Lat: ", lat, "Lon: ", lon)  # Popup básico
      )
  })
  
  # Actualizar el mapa con el punto de interés y el buffer
  observe({
    req(punto_interes_sf())
    
    leafletProxy("mapa") %>%
      clearGroup("punto_buffer") %>%  # Limpiar capas anteriores
      setView(lng = st_coordinates(punto_interes_sf())[1],
              lat = st_coordinates(punto_interes_sf())[2],
              zoom = 15)
  })
  
  # Mostrar los datos filtrados en el mapa
  observe({
    req(datos_filtrados())
    
    leafletProxy("mapa") %>%
      clearGroup("datos_filtrados") %>%  # Limpiar capas anteriores
      addPolygons(
        data = datos_filtrados(),
        fillColor = "yellow",
        fillOpacity = 0.5,
        color = "darkgreen",
        weight = 2,
        group = "datos_filtrados",
        popup = ~paste("Trimestre:", yearqtr, "<br>Área:", area)
      )
  })
}

# Ejecutar la aplicación
shinyApp(ui, server)


#Once the database has been reviewed, it is saved
saveRDS(sf2, "./MP_data")

##################################################


