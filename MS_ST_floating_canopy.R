# Spatio-temporal modeling of the extent of the floating canopy 
# of Macrocystis pyrifera along the Humboldt Current,
# in seasonal intervals from summer 2016 to summer 2025.

## library and wd #### ####
library(sdmTMB)
library(INLA)
library(dplyr)
library(tidyr)
#library(purrr)
library(sf)
library(ggplot2)
library(gstat)
library(raster)
library(sdmpredictors)
library(sp) 
library(geoR) 
library(leaflet)
library(shiny)
library(shinyWidgets)


# Establecer el número de decimales a mostrar
options(digits = 4)

setwd("(...)/GitHub/Macrocystis_pyrifera/") # The local working directory should be added here.
###########################

####### Polygons and spatial domain ############# ####
chile <- read_sf("./Dependences/Shapes_chile/Chile.shp")
sf_use_s2(FALSE)
chile <- chile %>% st_as_sf() %>% st_make_valid() %>% st_crop(xmin= -76, ymin= -41.56, xmax= -69, ymax= -18) #-18
chile_utm <- chile %>% st_transform("+proj=utm +zone=19 + south +ellps=WGS84 +datum=WGS84 +units=km")
ecoregions <- read_sf("./Dependences/MEOW/meow_ecos.shp") 
ecoregions <- ecoregions %>% st_as_sf() %>% st_make_valid() %>% st_crop(xmin= -76, ymin= -41.56, xmax= -69, ymax= -18)
ecoregions <- filter(ecoregions, ECOREGION %in% c("Humboldtian", "Central Chile", "Araucanian"))

#setting the coastal spatial domain
marspec <- load_layers(c("MS_biogeo05_dist_shore_5m"), datadir = "./Dependences/MARSPEC")
ext <- extent(-76, -66, -41.56, -18) #Study area
marspec <- crop(marspec, ext)
dist_pol <- rasterToPolygons(marspec$MS_biogeo05_dist_shore_5m, fun=function(x){x <= 100}, dissolve = T)
dist_pol <- buffer(dist_pol, 2000) # to include all the geometry of the coast which is sometime not well represented in the marspec database
dist_pol <- aggregate(dist_pol, fact= 1, dissolve= T)
dist_pol <- st_as_sf(dist_pol)


eco_zee <- st_intersection(ecoregions, dist_pol) #to crop the marine ecoregions to the area of interest.

eco_zee$ECOREGION <- factor(eco_zee$ECOREGION, levels = c("Humboldtian","Central Chile",  "Araucanian"))
eco_zee_utm <- eco_zee %>% st_transform("+proj=utm +zone=19 + south +ellps=WGS84 +datum=WGS84 +units=km")

ggplot()+
  geom_sf(data=eco_zee, aes(fill= ECOREGION))+
  geom_sf(data= chile,aes(geometry= geometry))+
  theme(axis.text.x = element_text(angle= 45))


######################################################

## data base ######## ####
sf <- readRDS("./MP_data")
df <- st_drop_geometry(sf)
summary(df)
#########################

##### histograms ######## ####
# 25% of the data belong to an unique area value, which is due to single pixels.
h1 <- ggplot(data = df[df$area < 10000,]) + 
  geom_histogram(aes(x = area),
                 fill = "#1E88E5",          # Modern blue fill color
                 color = "white",           # White borders
                 alpha = 0.85,              # Slight transparency
                 bins = 100) +               # Adjust number of bins as needed
  labs(x = expression("Area ("*m^2*")"),    # Proper m² notation
       y = "Count") +
  annotate(geom= "text",x= 10000, y= 1500, label= "A)")+
  theme_minimal(base_size = 12) +           # Clean theme with larger font
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank())
h2 <-  ggplot(data = df[df$area > 10000 & df$area < 30000,]) + 
  geom_histogram(aes(x = area),
                 fill = "#1E88E5",          # Modern blue fill color
                 color = "white",           # White borders
                 alpha = 0.85,              # Slight transparency
                 bins = 100) +               # Adjust number of bins as needed
  labs(x = expression("Area ("*m^2*")"),    # Proper m² notation
       y = "Count") +
  annotate(geom= "text",x= 30000, y= 40, label= "B)")+
  theme_minimal(base_size = 12) +           # Clean theme with larger font
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank())
h3 <-  ggplot(data = df[df$area > 30000,]) + 
  geom_histogram(aes(x = area),
                 fill = "#1E88E5",          # Modern blue fill color
                 color = "white",           # White borders
                 alpha = 0.85,              # Slight transparency
                 bins = 100) +               # Adjust number of bins as needed
  labs(x = expression("Area ("*m^2*")"),    # Proper m² notation
       y = "Count") +
  annotate(geom= "text",x= 180000, y= 25, label= "C)")+
  theme_minimal(base_size = 12) +           # Clean theme with larger font
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank())

library(patchwork)
h1+h2+h3 + plot_annotation(title = "Floating canopy extension distribution",
                           subtitle = "Frequency histogram",
                           theme = theme(plot.title = element_text(size = 14, face = "bold")))

#########################

#### data grid and group ####### ####

NGroup <- seq(from= 1, to= length(unique(df$yearqtr)), by= 1)
fecha_group <- data.frame(Group= NGroup, yearqtr= sort(unique(df$yearqtr)))
head(fecha_group)
df <- left_join(df, fecha_group, by= "yearqtr")

# 
#Grouped data base of a 5 km grid 
lim_lat <- range(df$Y)    # latitude range
lim_lon <- range(df$X)    # longitude range
lat_axis <- seq(lim_lat[1],   # latitude axis grid
                lim_lat[2],
                by= 5) #~5 Km
lon_axis <- seq(lim_lon[1],   # latitude axis grid
                lim_lon[2],
                by= 5)

# Completa la grilla con los datos existentes y llena NAs con ceros

df_grid <- df    
dists <- abs(outer(df$Y, lat_axis, "-"))
dists_x <- abs(outer(df$X, lon_axis, "-"))
df_grid$Y <- lat_axis[apply(dists, 1, which.min)]
df_grid$X <- lon_axis[apply(dists_x, 1, which.min)]
df_grid$lon <- round(df$lon, digits = 1)

df_grouped <- group_by(df_grid, X, Y, Group, ECOREGION) %>%
  summarise(z =sum(area), lon= mean(lon), lat= mean(lat)) %>%
  mutate(lz= log(z)) %>%
  left_join(fecha_group, by= "Group")
df <- df %>% mutate(lz= log(area))   # To also have the log of each individual observation
################################

#### Interactive map for database visualization ##### ####
loc <- df %>%
  mutate(lat2= round(lat, 2)) %>%
  dplyr::select(yearqtr, lat, lat2, lon, area) %>%
  group_by(lat2) %>%
  filter(area== max(area)) %>%
  ungroup()

# Datos de los puntos (longitud, latitud, nombre)
upw_points <- data.frame(
  lon = c(-70.5, -71.6, -71.7, -73.5),  # Longitudes (Oeste)
  lat = c(-23.5, -29.25, -32.7, -37.1),  # Latitudes (Sur)
  name = c("Pta de Antofagasta", "Pta de Choros", "Valparaiso", "Pta Lavapíe")
)
# View(loc)
# Crear objeto sf a partir de las coordenadas
upw_sf <- st_as_sf(upw_points, 
                   coords = c("lon", "lat"),
                   crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
upw_sf <- st_intersection(upw_sf, eco_zee[,"ECOREGION"])
upw_points<- left_join(upw_points, st_drop_geometry(upw_sf), by="name")

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
    sf2_zoom <- st_intersection(sf, buffer_zoom())  # Recortar sf2 con el buffer
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


###########################

#### empirical spatial and temporal average ##### ####
spat_av <- group_by(df_grouped, Y, lat) %>% # group by lon-lat
  summarise(mu_emp = mean(lz)) # mean for each lon-lat
T_av <- group_by(df_grouped, yearqtr, Group) %>%
  summarise(meanA = mean(lz))

ecorreg_pos <- df_grouped %>%
  group_by(ECOREGION) %>%
  summarise(lat_mean = mean(lat))

hov<- ggplot(df_grouped) + # take data
    geom_tile(aes(x = yearqtr, y = lat, color= lz, linewidth= lz)) +   # plot
    scale_color_viridis_c(expression(paste("Log(area ", m^2, ")"))) +
    labs(linewidth= element_blank())+
    xlab("Year") + # add y label
    ylab("Latitude") + # add x label
    theme_bw()+
    annotate("text", x= zoo::as.yearqtr("2025 Q1"), y= -18, label= "C)", size= 4)+
    guides(linewidth = "none" # Elimina la leyenda del linewidth
    ) + 
    # Líneas de límite de ecorregiones (mejoradas)
    geom_hline(yintercept = c(-25, -33),
               color = "red" ,  # Colores más profesionales (RColorBrewer)
               linetype = "longdash",
               linewidth = 1.2,
               alpha = 0.7)+
    theme(
      legend.box.margin = margin(t = 0, r = 0, b = 0, l = 15),
      axis.text = element_text(size= 12),
      legend.text.position = "right",
      legend.box = "horizontal",
      legend.spacing =  unit(2, "cm"),  # Espacio horizontal extra  
      legend.key.size = unit(0.8, 'cm'),
      legend.text = element_text(size=10),
      legend.title = element_text(size=10))

(hov_label <- hov + 
    geom_text(
      data = ecorreg_pos,
      aes(
        x = max(df_grouped$yearqtr) + 0.5,  # Posición a la derecha del último dato
        y = lat_mean,
        label = ECOREGION,
        color = NULL  # Anulamos el mapeo de color
      ),
      hjust = 0.5,       # Alineación horizontal (0 = izquierda)
      vjust = -0.5,      # Alineación vertical (0.5 = centro)
      size = 4,
      fontface = "bold",
      check_overlap = TRUE,
      angle= -90
    ) +
    coord_cartesian(
      clip = "off",     # Permite que las etiquetas salgan del área del gráfico
      xlim = range(df_grouped$yearqtr)  # Mantenemos el límite original en X
    ) +
    theme(
      plot.margin = margin(r = 10)  # Aumentamos el margen derecho para las etiquetas
    ))

sa <- ggplot(spat_av) +
    geom_point(aes(lat, mu_emp), alpha= 0.02) +
    #Líneas de límite de ecorregiones (mejoradas)
    geom_vline(xintercept = c(-25, -33),
               color = "red" ,  # Colores más profesionales (RColorBrewer)
               linetype = "longdash",
               linewidth = 1.2,
               alpha = 0.7)+
    geom_smooth(aes(lat, mu_emp), method = "gam",  
                formula = y ~ s(x), 
                colour = "black")+
    xlab("Latitude") +
    annotate("text", x= -19, y= 12.5, label= "A)", size= 4)+
    ylab(expression(paste("Log(area ", m^2, ")"))) + theme_bw() + 
    theme(axis.text = element_text(size= 12))
(sa_label <- sa + geom_text(
  data = ecorreg_pos,
  aes(
    y = max(spat_av$mu_emp),  # Posición a la derecha del último dato
    x = lat_mean,
    label = ECOREGION,
    color = NULL  # Anulamos el mapeo de color
  ),
  hjust = 0.5,       # Alineación horizontal (0 = izquierda)
  vjust = -2,      # Alineación vertical (0.5 = centro)
  size = 4,
  fontface = "bold") + 
    coord_cartesian(
      clip = "off") +
    theme(
      plot.margin = margin(t = 20)  # Aumentamos el margen derecho para las etiquetas
    )  # Permite que las etiquetas salgan del área del gráfico
)
(ta <- ggplot() +
    geom_point(data = T_av,aes(x = yearqtr, y = meanA)) +
    geom_smooth(data= T_av, aes(x=yearqtr, meanA), colour= "black") +
    geom_line(data = T_av, aes(x = yearqtr, y = meanA)) +
    annotate("text", x= zoo::as.yearqtr("2025 Q1"), y= 9.3, label= "B)", size= 4)+
    xlab("Year") + ylab(expression(paste("Log(area ", m^2, ")"))) +
    theme_bw()+
    theme(axis.text = element_text(size= 12)))
library(patchwork)
(sa_label/ta)|hov_label
#################################################

##### sdmTMB null modeling ###################### ####
fit0 <- sdmTMB(lz ~ 1,
               data = df_grouped,
               spatial= "off",
               family = gaussian(link= "identity")
)
summary(fit0)    
BIC(fit0)       # 16729
preds0 <- fit0$family$linkinv(predict(fit0)$est)
# Pseudo-R²
observed <- df_grouped$lz
(pseudo_R2 <- 1 - (var(observed - preds0) / var(observed))) #0


fit0 <- sdmTMB(lz ~ 0+ECOREGION,
               data = df_grouped,
               #mesh = mesh,
               spatial= "off",
               family = gaussian(link= "identity")
)
summary(fit0)    
BIC(fit0)       # 16583

preds0 <- fit0$family$linkinv(predict(fit0)$est)
# Pseudo-R²
observed <- df_grouped$lz
(pseudo_R2 <- 1 - (var(observed - preds0) / var(observed))) #0.032

res0 <- residuals(fit0, type= "mle-mvn")
df_grouped$E0 <- res0

##### Variogram 
# Variogram 
MyData <- data.frame(E0 = res0,
                     Xkm = df_grouped$X, 
                     Ykm = df_grouped$Y,
                     lz= df_grouped$lz,
                     Group= df_grouped$Group)

# Create a "gstat" object 
coordinates(MyData) <- ~ Xkm + Ykm  #effectively convert the data into a spatial data frame
proj4string(MyData) <- "+proj=utm +zone=19 +south +ellps=WGS84 +units=km"
TheGStat_res <- gstat::gstat(id="Variogram residuals", formula=E0 ~ 1, data=MyData)

Vario0 <- variogram(object = E0 ~ 1 ,  
                    data = MyData,  
                    cressie = T,  
                    cutoff = 500, width = 5) 
Vario0_local <- variogram(object = E0 ~ 1,  
                          data = MyData,  
                          cressie = F,  
                          cutoff = 50, width = 2) 
#plot(variogram(TheGStat_res, map=T, cutoff=200, width=10))

(Vs0 <- ggplot(data = Vario0,  
               aes(x = dist, y = gamma)) +
    geom_point() + 
    geom_smooth(method = "gam",  
                formula = y ~ s(x,  
                                bs = "cs"), 
                colour = "black") +
    annotate("text", y= 1.25, x= 500, label= "B)", size= 5)+
    #ylim(0,1) + 
    theme(axis.text = element_text(size = 12),
          text= element_text(size= 12),
          axis.text.x = element_text(angle= 45, vjust= 0.5)) +  
    xlab("Distance (km)") +
    ylab("Variogram"))

(Vs0_local <- ggplot(data = Vario0_local,  
                     aes(x = dist, y = gamma)) +
    geom_point() + 
    geom_smooth(method = "gam",  
                formula = y ~ s(x,  
                                bs = "cs"), 
                colour = "black") +  
    annotate("text", y= 1.25, x= 50, label= "A)", size= 5)+
    #ylim(0,1) + 
    theme(axis.text = element_text(size = 12),
          text= element_text(size= 12)) +  
    xlab("Distance (km)") +
    ylab("Variogram"))


library(sp) 
library(geoR) 
MyData$Zeros  <- rep(0, nrow(MyData)) 
Var_df       <- as.data.frame(na.exclude(MyData[,c("E0", "Group","Zeros")]))
breaks <- seq(1, 16, length = 100) 
V0 <- geoR::variog(coords = Var_df[,c("Group", "Zeros")],  
                   data = Var_df[,"E0"], 
                   option= "bin",
                   breaks = breaks) 
# Convertir resultados a dataframe para ggplot
Vario_t0 <- data.frame(
  dist = V0$u,  # Distancias temporales
  gamma = V0$v  # Semivarianza
)
# Graficar usando ggplot
(Vt0 <- ggplot(data = Vario_t0, aes(x = dist, y = gamma)) +
    geom_line() +
    geom_point(size= 3)+
    scale_x_continuous(breaks = 0:16)+
    theme_bw() +
    ylim(c(0.9, 1.02))+
    xlab("Temporal Interval") +
    annotate("text", y= 1.015, x= 16, label= "C)", size= 5)+
    ylab("Variogram")+
    theme(#text = element_text(size = 15),
      axis.text = element_text(size= 12),
      text = element_text(size= 12),
      panel.grid.minor.x= element_blank()))
acf(MyData$E0, lag.max= 15)

plot(V0,  
     type = "b", 
     xlab = "Temporal interval", 
     ylab = "Semi-variance",
     #ylim=c(0.085,0.095),
     col = 1, pch = 16, cex = 1.5) 

#################################################

######## Model fit sdmTMB ###################### #####

##### fit spatial model and spatiotemporal random field ar1

mesh <- make_mesh(
  df_grouped,
  xy_cols = c("X", "Y"),
  cutoff = 1) # min. distance between knots in X-Y units
plot(mesh) 
mesh$mesh$n

fit1 <- sdmTMB(lz ~ 0+ECOREGION,
               data = df_grouped,
               mesh = mesh,
               spatial= "on",
               family = gaussian(link = "identity"))
summary(fit1)    
BIC(fit1) #14177
tidy(fit1, effects = "ran_pars", conf.int = T)

preds1 <- fit1$family$linkinv(predict(fit1)$est)
# Pseudo-R²
observed <-df_grouped$lz
(pseudo_R2 <- 1 - (var(observed - preds1) / var(observed))) # 0.529

## Residuals
res1 <- residuals(fit1, type= "mle-mvn")
MyData$E1 <- res1
df_grouped$E1 <- res1
summary(MyData)


Vario1 <- variogram(object = E1 ~ 1,  
                    data = MyData,  
                    cressie = TRUE,  
                    cutoff = 500, width = 5) 
Vario1_local <- variogram(object = E1 ~ 1,  
                          data = MyData,  
                          cressie = TRUE,  
                          cutoff = 50, width = 2) 

(Vs1 <- ggplot(data = Vario1,  
               aes(x = dist, y = gamma)) +
    geom_point() + 
    geom_smooth(method = "gam",  
                formula = y ~ s(x,  
                                bs = "cs"), 
                colour = "black") +  
    annotate("text", y= 1.1, x= 450, label= "B)", size= 5)+
    #ylim(0,1) + 
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle= 45, vjust= 0.5),
          text= element_text(size= 12)) +  
    xlab("Distance (km)") +
    ylab("Variogram"))

(Vs1_local <- ggplot(data = Vario1_local,  
                     aes(x = dist, y = gamma)) +
    geom_point() + 
    geom_smooth(method = "gam",  
                formula = y ~ s(x,  
                                bs = "cs"), 
                colour = "black") +  
    annotate("text", y= 1.1, x= 45, label= "A)", size= 5)+
    #ylim(0,1) + 
    theme(axis.text = element_text(size = 12),
          text= element_text(size= 12)) +   
    xlab("Distance (km)") +
    ylab("Variogram"))

#MyData$Zeros <- 0
Var_df       <- as.data.frame(na.exclude(MyData[,c("E1", "Group","Zeros")]))
breaks <- seq(1, 16, by = 1) 
V1 <- geoR::variog(coords = Var_df[,c("Group", "Zeros")],  
                   data = Var_df[,"E1"], 
                   breaks = breaks) 
# Convertir resultados a dataframe para ggplot
Vario_t1 <- data.frame(
  dist = as.integer(V1$u),  # Distancias temporales
  gamma = V1$v  # Semivarianza
)

# Graficar usando ggplot
(Vt1 <- ggplot(data = Vario_t1, aes(x = dist, y = gamma)) +
    geom_line() +
    geom_point(size= 3)+
    scale_x_continuous(breaks = 0:16)+
    theme_bw() +
    ylim(c(0.8, 1.1))+
    xlab("Temporal Interval") +
    annotate("text", y= 1.1, x= 15, label= "C)", size= 5)+
    ylab("Variogram")+
    theme(#text = element_text(size = 15),
      axis.text = element_text(size= 12),
      text = element_text(size= 12),
      panel.grid.minor.x= element_blank()))


### Espaciotemporal ##

fit2 <- sdmTMB(lz ~ 0 + ECOREGION,
               data = df_grouped,
               mesh = mesh,
               spatial= "on",
               family = gaussian(link = "identity"),
               time= "Group",
               share_range = F,
               spatiotemporal = "ar1",
               priors =  sdmTMBpriors(matern_st = pc_matern(range_gt = 100, sigma_lt = 1, range_prob = 0.01)),
               control = sdmTMBcontrol(newton_loops = 1))
summary(fit2) # 
tidy(fit2, effects = "ran_pars", conf.int = T)
BIC(fit2)     # 13276
sanity(fit2)
preds2 <- fit2$family$linkinv(predict(fit2)$est)
# Pseudo-R²
observed <- df_grouped$lz
(pseudo_R2 <- 1 - (var(observed - preds2) / var(observed)))

saveRDS(fit2, "./Spatiotemporal/stmodel_full_eco_shareF")


res2 <- residuals(fit2, type= "mle-mvn")

# Variogram
MyData$SR <- res2
df_grouped$SR <- res2

Vario3 <- variogram(object = SR ~ 1,  
                    data = MyData,  
                    cressie = TRUE,  
                    cutoff = 500, width = 5) 

Vario3_local <- variogram(object = SR ~ 1,  
                          data = MyData,  
                          cressie = TRUE,  
                          cutoff = 50, width = 2) 


(Vs3 <- ggplot(data = Vario3,  
               aes(x = dist, y = gamma)) +
    geom_point() + 
    geom_smooth(method = "gam",  
                formula = y ~ s(x,  
                                bs = "cs"), 
                colour = "black") +  
    annotate("text", y= 1.1, x= 500, label= "B)", size= 5)+
    #ylim(0,1) + 
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle= 45, vjust= 0.5),
          text= element_text(size= 12)) + 
    xlab("Distance (km)") +
    ylab("Variogram"))

(Vs3_local <- ggplot(data = Vario3_local,  
                     aes(x = dist, y = gamma)) +
    geom_point() + 
    geom_smooth(method = "gam",  
                formula = y ~ s(x,  
                                bs = "cs"), 
                colour = "black") +  
    annotate("text", y= 1.05, x= 50, label= "A)", size= 5)+
    #ylim(0,1) + 
    theme(axis.text = element_text(size = 12),
          text= element_text(size= 12)) +   
    xlab("Distance (km)") +
    ylab("Variogram"))

#MyData$Zeros  <- rep(0, nrow(MyData)) 
Var_df       <- as.data.frame(na.exclude(MyData[,c("SR", "Group","Zeros")]))
breaks <- seq(1, 16, by= 1) 
V3 <- geoR::variog(coords = Var_df[,c("Group", "Zeros")],  
                   data = Var_df[,"SR"], 
                   breaks = breaks) 

# Convertir resultados a dataframe para ggplot
Vario_t3 <- data.frame(
  dist = as.integer(V3$u),  # Distancias temporales
  gamma = V3$v  # Semivarianza
)
(Vt3 <- ggplot(data = Vario_t3, aes(x = dist, y = gamma)) +
    geom_line() +
    geom_point(size= 3)+
    scale_x_continuous(breaks = 1:15)+
    theme_bw() +
    ylim(c(0.8, 1.1))+
    xlab("Temporal Interval") +
    ylab("Variogram")+
    theme(panel.grid.minor.x= element_blank()))

(Vt3 <- ggplot(data = Vario_t3, aes(x = dist, y = gamma)) +
    geom_line() +
    geom_point(size= 3)+
    scale_x_continuous(breaks = 0:16)+
    theme_bw() +
    ylim(c(0.8, 1.1))+
    xlab("Temporal Interval") +
    annotate("text", y= 1.1, x= 15, label= "C)", size= 5)+
    ylab("Semi-variance")+
    theme(#text = element_text(size = 15),
      axis.text = element_text(size= 12),
      text = element_text(size= 12),
      panel.grid.minor.x= element_blank()))


acf(MyData$E1, lag.max= 15)
acf(MyData$SR, lag.max= 15)

################################################

#### spacetime variograms ###################### #### 
time_lat_cov <- df_grouped %>% 
  ungroup %>%
  expand(Y, X, Group) %>%
  left_join(df_grouped, by = c("Y","X", "Group")) %>%
  #dplyr::select(-ECOREGION) %>%
  dplyr::select(-c("yearqtr")) %>%
  left_join(fecha_group, by= c("Group")) %>%
  #replace_na(list(SR = 0)) %>%
  arrange(Y, Group) %>%
  mutate(lon = ifelse(is.na(lon), zoo::na.approx(lon), lon),
         lat = ifelse(is.na(lat), zoo::na.approx(lat), lat))

summary(time_lat_cov)

# Ordena time_lat_Hov por fecha
time_lat_cov <- time_lat_cov[order(time_lat_cov$yearqtr),]
spat_df <- filter(time_lat_cov, Group ==22) %>% # lon/lat coords of stations
  dplyr::select(lat, lon) %>% # select lon/lat only
  arrange(lon, lat) # sort ascending by lon/lat
summary(spat_df)
spat_part <- SpatialPoints(coords = spat_df[, c("lon", "lat")], proj4string=CRS("+proj=longlat +datum=WGS84"))
temp_part <- sort(unique(as.Date(as.POSIXct(zoo::as.yearqtr(time_lat_cov$yearqtr)))))
STObj3 <- spacetime::STFDF(sp = spat_part,
                           time = temp_part,
                           data = time_lat_cov)

vvE0 <- variogramST(formula = E0~ 1, # fixed effect component
                    data = STObj3, # July data
                    width = 5, # spatial bin (5 km)
                    cutoff = 500, # consider pts < 500 km apart
                    tlags = 0:16, 
                    covariogram= F) 
vvE1 <- variogramST(formula = E1~ 1, # fixed effect component
                    data = STObj3, # July data
                    width = 5, # spatial bin (5 km)
                    cutoff = 500, # consider pts < 500 km apart
                    tlags = 0:16, 
                    covariogram= F) 

vvSR <- variogramST(formula = SR~ 1, # fixed effect component
                    data = STObj3, # July data
                    width = 5, # spatial bin (5 km)
                    cutoff = 500, # consider pts < 500 km apart
                    tlags = 0:16, 
                    covariogram= F) 

plot(vvE1)
vvdfE0 <- as.data.frame(vvE0)
vvdfE0 <- filter(vvdfE0, dist > 0 & !is.na(id) & id != "lag0")
vvdfE1 <- as.data.frame(vvE1)
vvdfE1 <- filter(vvdfE1, dist > 0 & !is.na(id) & id != "lag0")
vvdfSR <- as.data.frame(vvSR)
vvdfSR <- filter(vvdfSR, dist > 0 & !is.na(id) & id != "lag0")


(st_varE0 <- ggplot(vvdfE0, aes(x = dist, y = reorder(id, timelag), fill = gamma)) +
    geom_tile(color = NA, width = 5.5, height = 1) +  # Gráfico de calor
    scale_fill_viridis_c(name = "Variogram") +  # Escala continua de colores
    labs(
      x = "Distance (km)",
      y = "Time lag (Seasonal interval)",
      title = "Space-Time Variogram Null model"
    ) +
    theme_minimal()+
    annotate("text", y= "lag16", x= 480, label= "D)", color= "white", size= 5)+
    scale_x_continuous(breaks = seq(0, 500, by = 50))+  
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.1,
      title.vjust = 3))+
    theme(legend.margin = margin(l = -10, unit = "pt"),
          axis.text = element_text(size= 12),
          panel.grid = element_blank(),
          text = element_text(size= 12),
          panel.grid.minor.x= element_blank())) 
(st_varE1 <- ggplot(vvdfE1, aes(x = dist, y = reorder(id, timelag), fill = gamma)) +
    geom_tile(color = NA, width = 5.5, height = 1) +  # Gráfico de calor
    scale_fill_viridis_c(name = "Variogram") +  # Escala continua de colores
    labs(
      x = "Distance (km)",
      y = "Time lag (Seasonal interval)",
      title = "Space-Time Variogram Spatial Model"
    ) +
    theme_minimal()+
    annotate("text", y= "lag16", x= 480, label= "D)", color= "white", size= 5)+
    scale_x_continuous(breaks = seq(0, 500, by = 50))+  
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.1,
      title.vjust = 3))+
    theme(legend.margin = margin(l = -10, unit = "pt"),
          axis.text = element_text(size= 12),
          panel.grid = element_blank(),
          text = element_text(size= 12),
          panel.grid.minor.x= element_blank())) 

(st_varSR <- ggplot(vvdfSR, aes(x = dist, y = reorder(id, timelag), fill = gamma)) +
    geom_tile(color = NA, width = 5.8, height = 1) +  # Gráfico de calor
    scale_fill_viridis_c(name = "Variogram") +  # Escala continua de colores
    labs(
      x = "Distance (km)",
      y = "Time lag (Seasonal interval)",
      title = "Space-Time Variogram Spatio-temporal model"
    ) +
    theme_minimal()+
    annotate("text", y= "lag16", x= 480, label= "D)",  color= "white", size= 5)+
    scale_x_continuous(breaks = seq(0, 500, by = 50))+  
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.1,
      title.vjust = 3))+
    theme(legend.margin = margin(l = -10, unit = "pt"),
          axis.text = element_text(size= 12),
          panel.grid = element_blank(),
          text = element_text(size= 12),
          panel.grid.minor.x= element_blank())) 


plot(st_var)

################################################

#### Residual exploration ###################### ####
p <- predict(fit2)
df_grouped$est <- p$est
library(qqplotr)
# qqplot
(qq.fit=ggplot(data=df_grouped, mapping=aes(sample = SR)) + 
    stat_qq_point()+
    theme_bw()+
    stat_qq_line()+
    stat_qq_band(alpha=0.3)+
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
    theme(axis.title=element_text(size=12), 
          axis.text = element_text(size=12)))

# Residuales vs fitted
(res.fit.fit= ggplot(data=df_grouped, aes(x=est,y=SR))+ 
    geom_point(size=2)+
    geom_hline(yintercept =0,linetype = 2, size=1.2)+
    theme_bw()+ 
    labs(x="Fitted",y="Residuals")+
    theme(axis.title=element_text(size=12), 
          axis.text = element_text(size=12),
          legend.position = "none"))
# Residuales vs yearqtr
(res.yearqtr= ggplot(data=df_grouped, aes(x=yearqtr,y=SR, group= yearqtr))+ 
    geom_boxplot()+
    geom_hline(yintercept =0,linetype = 2, size=1.2)+
    theme_bw()+ 
    labs(x="Year",y="Residuals")+
    theme(axis.title=element_text(size=12), 
          axis.text = element_text(size=12),
          legend.position = "none"))
# Residuales vs yearqtr
(res.lat= ggplot(data=df_grouped, aes(x=lat,y=SR))+ 
    geom_point(size= 2)+
    geom_hline(yintercept =0,linetype = 2, size=1.2)+
    theme_bw()+ 
    labs(x="Latitude",y="Residuals")+
    theme(axis.title=element_text(size=12), 
          axis.text = element_text(size=12),
          legend.position = "none"))
library(patchwork)
qq.fit|res.fit.fit|res.yearqtr|res.lat


ggplot() + 
  geom_point(data = df_grouped, aes(x= lon, y= lat, color = SR, size= abs(SR)), stroke= 2, fill= NA, shape = 21, show.legend = T) + 
  #geom_sf(data= sf, aes(geometry=geometry, color = epsilon_st, fill= epsilon_st), show.legend = T) +
  #scale_fill_gradient2(name = ) +
  scale_color_gradient2(name = "Residuals") +
  #labs(title= "Spatiotemporal Random Field")+
  ggnewscale::new_scale_color() +
  geom_sf(data= eco_zee, aes(color= ECOREGION), linewidth=0.5, fill= NA, show.legend = T)+
  geom_sf(data= chile, fill= "white", color= "grey2")+
  #scale_color_brewer()+
  ylab("") + xlab("") + 
  labs(color= "Ecoregion")+
  guides(color = guide_legend(override.aes = list(alpha = 0.01)),
         size= "none") + 
  theme_bw()+
  facet_wrap(~yearqtr, nrow= 3)+
  theme(axis.text.x = element_blank())

################################################

#### plot conditional effect Ecoregions ######## ####

eco_eff <- ggeffects::ggeffect(fit2, terms = "ECOREGION") %>% as.data.frame()
datos_observados <- df_grouped %>%
  dplyr::select(ECOREGION, lz) %>%           # Seleccionar solo ecorregión y la variable modelada
  rename(x = ECOREGION, observed = lz)  # Renombrar para coincidir con eco_eff
datos_completos <- left_join(eco_eff, datos_observados, by = "x")

(cond <- ggplot() +
    # Datos observados (violin plot)
    geom_violin(
      data = datos_observados,
      aes(x = x, y = observed, fill = x), show.legend = F,
      trim = T
    ) +
    geom_jitter(
      data = datos_observados,
      aes(x = x, y = observed, fill= x), alpha= 0.05, show.legend = F
    )+ 
    labs(fill= "Ecoregion")+
    # Predicciones del modelo (puntos rojos)
    geom_point(
      data = eco_eff,
      aes(x = x, y = predicted),
      color = "black",
      size = 4,
      shape = 18
    ) +
    # Intervalos de confianza (barras de error)
    geom_errorbar(
      data = eco_eff,
      aes(x = x, ymin = conf.low, ymax = conf.high),
      width = 0.2,
      color = "black"
    ) +
    labs(
      x = "Ecoregion",
      y = expression(paste("Log(area ", m^2, ")"))
    ) +
    #annotate("text", y= 13, x= 3.5, label= "A)", size= 4)+
    theme_bw()+
    theme(axis.text = element_text(size= 12),
          axis.title = element_text(size= 12)))
################################################

#### Spatial and Spatio-temporal effects ####### ####
#grid
#1 - spatial domain
library(sdmpredictors)
marspec <- load_layers(c("MS_biogeo05_dist_shore_5m"), datadir = "./Capas_climaticas/FRE_selected/MARSPEC")
ext <- raster::extent(-78, -68, -42, -18) #Study area # -33 para el limite con la ecoregion Araucanian
marspec <- raster::crop(marspec, ext)
dist_pol <- rasterToPolygons(marspec$MS_biogeo05_dist_shore_5m, fun=function(x){x <= 50}, dissolve = T)
dist_pol <- aggregate(dist_pol, fact= 1, dissolve= T)
m <- dist_pol %>%
  st_as_sf() %>%
  st_cast("POLYGON")

m <- m %>% st_transform("+proj=utm +zone=19 + south +ellps=WGS84 +datum=WGS84 +units=km")
ggplot(m) + geom_sf() + theme_bw() + coord_sf(datum = st_crs(m))
bnd <- as(m, "Spatial")
bnd.df <- fortify(bnd)
bnd.df <- data.frame(X= bnd.df$long, Y= bnd.df$lat)
w.proj <- inla.mesh.projector(mesh$mesh,
                              xlim = range(bnd.df[,1]),
                              ylim = range(bnd.df[,2]),
                              dims = c(20, 500))  #from boudaries of the map polygon 
grid <- expand.grid(X = w.proj$x, 
                    Y =w.proj$y)
grid_sf <- st_as_sf(grid, coords = c("X", "Y"), crs="+proj=utm +zone=19 + south +ellps=WGS84 +datum=WGS84 +units=km")
grid_sf <- st_intersection(grid_sf, m)
grid_sf <- st_intersection(grid_sf, eco_zee_utm)
grid_sf$X <- st_coordinates(grid_sf)[,1]
grid_sf$Y <- st_coordinates(grid_sf)[,2]

grid_sf$ECOREGION <- factor(grid_sf$ECOREGION, levels = c("Humboldtian","Central Chile",  "Araucanian"))
grid <- st_drop_geometry(grid_sf)
grid <- replicate_df(grid, "Group", unique(df$Group))
summary(grid)

#predict on grid
p <- predict(fit2, newdata = grid)
summary(p)
sf <- st_as_sf(p, coords = c("X", "Y"), crs="+proj=utm +zone=19 + south +ellps=WGS84 +datum=WGS84 +units=km")
sf$X <- st_coordinates(sf)[,1]
sf$Y <- st_coordinates(sf)[,2]
sf_latlon <- st_transform(sf, crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
sf$lon <- st_coordinates(sf_latlon)[,1]
sf$lat <- st_coordinates(sf_latlon)[,2]
sf <- left_join(sf, fecha_group, by= c("Group"))


# Datos de los puntos (longitud, latitud, nombre)
upw_points <- data.frame(
  lon = c(-70.5, -71.6, -71.7, -73.5),  # Longitudes (Oeste)
  lat = c(-23.5, -29.25, -32.7, -37.1),  # Latitudes (Sur)
  name = c("Pta de Antofagasta", "Pta de Choros", "Valparaíso", "Pta Lavapíe")
)

# Crear objeto sf a partir de las coordenadas
upw_sf <- st_as_sf(upw_points, 
                   coords = c("lon", "lat"),
                   crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Crear una versión modificada de 'name' con saltos de línea
upw_sf$name_wrapped <- stringr::str_wrap(upw_sf$name, width = 10)  # Ajusta 'width' según necesidad

#plot spatiotemporal effect: epsilon_st
ggplot() + 
  geom_sf(data= sf, aes(geometry=geometry, color = epsilon_st, fill= epsilon_st), show.legend = T) +
  scale_fill_gradient2(name = expression(epsilon[st])) +
  scale_color_gradient2(name = expression(epsilon[st])) +
  labs(title= "Spatiotemporal Random Field")+
  ggnewscale::new_scale_color() +
  geom_sf(data= eco_zee, aes(color= ECOREGION), linewidth=0.5, fill= NA, show.legend = T)+
  geom_sf(data= chile, fill= "white", color= "grey2")+
  ylab("") + xlab("") + 
  labs(color= "Ecoregion")+
  guides(color = guide_legend(override.aes = list(alpha = 0.01)), size= "none") + 
  theme_bw()+
  facet_wrap(~yearqtr, nrow= 3)+
  theme(axis.text.x = element_blank())

# plot spatial effects: omega_s
(map_o <- ggplot() + 
    geom_sf(data = sf[sf$Group==22,], aes(fill = omega_s, color= omega_s, geometry= geometry, size= abs(omega_s)), show.legend = T)+
    geom_sf(data= chile_utm, fill= "white")+
    scale_fill_gradient2(name= expression(omega[s]))+
    scale_color_gradient2(name= expression(omega[s]))+
    ggnewscale::new_scale_color() +
    geom_sf(data= eco_zee, aes(color= ECOREGION), fill= NA)+
    geom_sf(data=upw_sf, aes(geometry= geometry), color= "red", size= 2)+
    ggrepel::geom_label_repel(data= upw_sf, aes(geometry= geometry,label = name_wrapped),
                              stat = "sf_coordinates",  # Necesario para objetos sf
                              size = 3.5,
                              nudge_x= 600,
                              min.segment.length = 0,   # Siempre dibuja segmentos
                              segment.color = "gray50")+
    guides(size = "none", color = guide_legend(override.aes = list(alpha = 0.01))) + 
    ylab("") + xlab("") + guides(size= "none")+ labs(color= "Ecoregion") + 
    annotate("text", y= 8000, x= 760, label= "B)", size= 4)+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.text = element_text(size= 12),
          legend.text = element_text(size=12),
          legend.title = element_text(size= 12)))

(lat_o <-ggplot() + 
    geom_line(data = sf, aes(y = omega_s, x= lat), show.legend = T)+
    geom_point(data= upw_points, aes(y= 0, x= lat), color= "red", size= 3)+
    coord_flip()+
    geom_hline(yintercept= 0, 
               color = "red", linewidth=0.5) +
    xlab("Latitud") + ylab(expression(omega[s])) + 
    annotate("text", x= -18, y= 2.5, label= "A)", size= 4)+
    theme_bw()+
    theme(axis.text = element_text(size= 12),
          axis.title.x = element_text(size= 14)))
library(patchwork)
lat_o|map_o

#################################################

#### Ecoregional mean estimate ################## ####
div_index_lists <- list()

for(i in unique(grid$ECOREGION)){
  
  pred_grid_sub <- grid %>% filter(ECOREGION == i)
  
  pred_sim_sub <- predict(fit2, newdata = pred_grid_sub, nsim = 100, type= "link")
  
  index_sim <- get_index_sims(pred_sim_sub, 
                              area = rep(5, nrow(pred_sim_sub)), #given the 5X5 km grid
                              est_function = function(x) mean(x),
                              agg_function = function(x) mean(x)) %>% 
    mutate(Ecoregion = i)
  
  div_index_lists[[i]] <- index_sim
  
}

div_index_all <- bind_rows(div_index_lists)
div_index_all$Ecoregion <- factor(div_index_all$Ecoregion, levels = c("Humboldtian","Central Chile","Araucanian"))
div_index_all <- left_join(div_index_all, fecha_group, by= "Group")

filtered_data <- div_index_all %>%
  filter(grepl("Q1|Q3", yearqtr))  # Busca "Q1" O "Q3" en la columna
breaks <- unique(filtered_data$yearqtr)

summary(div_index_all)

ggplot(div_index_all, aes(yearqtr, est, color= Ecoregion)) +
  geom_line(linewidth= 0.8, alpha=0.8, show.legend = F) +
  geom_ribbon(alpha = 0.2, aes(ymin = lwr, ymax = upr, fill= Ecoregion, color= Ecoregion), show.legend = F)+
  scale_x_continuous(
    breaks = as.numeric(breaks),  
    labels = function(x) format(zoo::as.yearqtr(x), "%Y Q%q")  
  ) +
  facet_wrap(~Ecoregion, ncol= 1, scales= "free_y")+
  labs(y=expression(paste("Estimate Log(area  ",m^2,")")), x= NULL, title = "Ecoregional mean estimate", color= "Ecoregion", fill= "Ecoregion") +
  theme_minimal()+
  theme(panel.background = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle= 45, hjust= 1))

#################################################