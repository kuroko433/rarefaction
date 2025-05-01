
rm(list = ls())

setwd("C:/Users/halex/OneDrive/Documentos/Programacion/R scripts/R Diversidad/rarefaction/rarefaction")

library(velociraptr) ##manejar datos de pbdb
library(vegan) ##analisis ecologicos
library(dplyr) ##para manejar datos
library(ggplot2) ##para graficar
library(glue)
library(iNEXT) ##analisis ecologicos
##### datos de dromaeosaurios ###########
## vamos a descargar datos de raptores (familia de dinosaurios)
## a nivel de Genero
dromaeosaurios=downloadPBDB(Taxa="Dromaeosauridae",StartInterval="Jurassic",StopInterval="Cretaceous")

## quitamos ocurrencias que no llegan a genero
dromaeosaurios=subset(dromaeosaurios,!dromaeosaurios$genus=="")
dromaeosaurios<-dromaeosaurios[,-c(16:18)]

### veamos cuantos generos tenemos
collections=tapply(dromaeosaurios$collection_no,list(dromaeosaurios$collection_no,dromaeosaurios$genus),length)

##pasamos a matriz 1/0
collections[is.na(collections)]=0
collections[collections>1]=1

##numero de generos
ngenus<-ncol(collections)

##tenemos 54 generos de raptores
print(ngenus)

## veamos la completitud de nuestro ensamble (aproximacion tipo presencia-ausencia por coleccion)

specpool(collections)[[2]] ##basado en CHAO 2 abrian 47 generos totales para el Mioceno

print(glue("en total tenemos {specpool(collections)[1]} generos entre el jurasico y el cretacico,
           y una completitud de ensambles del {round((specpool(collections)[1]/specpool(collections)[2])*100,2)}%"))

print("lo anterior quiere decir que faltan más de la mitad de generos de raptores
      aun por encontrar")


###rarefaccion por coleccion 1

##calculamos diversidad rarefaccionada por coleccion (presencia-ausencia de generos)
dromaeo_rare<- specaccum(collections,method = "random", permutations = 500,
                      ci=0.95)

dromaeo_rare<-data.frame(rarefaction=dromaeo_rare$richness,
                         ci_lower=dromaeo_rare$richness-dromaeo_rare$sd,
                         ci_upper=dromaeo_rare$richness+dromaeo_rare$sd,
                         sites=dromaeo_rare$sites)

##### graficar resultado de vegan
rare_1<-ggplot() +
  # Graficar el intervalo de confianza como un área sombreada
  geom_ribbon(data = dromaeo_rare, aes(x = sites, ymin = ci_lower, ymax = ci_upper), 
              fill = "#698B22", alpha = 0.2) +
  # Graficar el efecto promedio como una línea
  geom_line(data = dromaeo_rare, aes(x = sites, y = rarefaction), 
            color = "#698B22", linewidth = 1) +
  # Etiquetas de los ejes y título
  labs(title = "",
       x = "Collections",
       y = "Rarefied Richness") +
  # Escala del eje X para incluir marcas en los extremos
  scale_x_continuous(breaks = c(0,seq(450, 0, by = -50),450)) +
  #Escala del eje Y para incluir marcas en los extremos
  scale_y_continuous(breaks = c(0,seq(120, 0, by = -40),120)) +
  ##agregar linea del valor de CHAO2
  geom_hline(yintercept = specpool(collections)[[2]] , linetype = "solid", color = "red",linewidth=1)+
  theme(
    panel.background = element_rect(fill = "white", colour = NA), # Fondo blanco sin bordes
    plot.background = element_rect(fill = "white", colour = NA),  # Fondo general blanco
    panel.grid = element_blank(),                                 # Eliminar cuadrícula
    axis.line = element_line(color = "black"),                    # Ejes X e Y como líneas negras
    axis.ticks = element_line(color = "black"),                   # Marcas de los ejes visibles
    axis.text = element_text(color = "black"),                    # Texto de los ejes en negro
    axis.title = element_text(color = "black")                    # Títulos de los ejes en negro
  )
rare_1

##rarefaccion por coleccion 2
dromaeo_list <- list(dromaeo = t(collections))


rare.2 <- iNEXT(dromaeo_list, 
             q = 0,  # Riqueza (q=0), Shannon (q=1), Simpson (q=2)
             datatype = "incidence_raw",  # Datos de presencia-ausencia
             endpoint = 800,  # Máximo número de unidades de muestreo, 
                              # puedes poner un numero mayor a los que tienes 
                              # para extrapolar al futuro (en que momento la curva se va a estabilizar?)
             knots = 100,  # Número de puntos para la curva
             se = TRUE,  # Calcular intervalos de confianza
             nboot = 500)  # Número de réplicas bootstrap

ggiNEXT(rare.2, type = 1)  # Curva basada en tamaño de muestra
ggiNEXT(rare.2, type = 3)  # Curva basada en cobertura
ggiNEXT(rare.2, type = 2)  # Curva de completitud muestral


