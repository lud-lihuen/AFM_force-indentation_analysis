# Procesar histograma de datos de modulo elastico de mapa QNM


# En software NanoScope usar la función Depth en el mapa de módulo elástico y exportar los datos XZ del histograma.
# Verificar que el X Axis este seteado en Absolute y no en Relative.


## Importar matriz ascii de resultados de modulo y frecuencias con Import Dataset > From Text (base) > Import
# Cambiar los nombres de los archivos a controlData y modifiedData (segun corresponda) al importar

# Metodo alternativo:
controlData <- read.delim("rutaArchivo")
modifiedData <- read.delim("rutaArchivo")


## Control group data

E_control <- controlData[1] # datos de módulo elástico
p_control <- controlData[2]/1 # porcentaje de cada dato (dividir por el numero de histogramas procesados)
frel_control <- p_control/100 # frecuencia relativa de cada dato
fabs_control <- ceiling(frel_control*nrow(controlData)) # frecuencia absoluta de cada dato


## Modified group data

E_modified <- modifiedData[1] # datos de módulo elástico
p_modified <- modifiedData[2]/1 # porcentaje de cada dato (dividir por el numero de histogramas procesados)
frel_modified <- p_modified/100 # frecuencia relativa de cada dato
fabs_modified <- ceiling(frel_modified*nrow(modifiedData)) # frecuencia absoluta de cada dato


## MEDIA PONDERADA

mean_control <- weighted.mean(E_control,frel_control)
mean_modified <- weighted.mean(E_modified,frel_modified)


## VARIANZA Y DESVIACION ESTANDAR

library(Hmisc) # inicializo package para usar estas funciones

var_control <- wtd.var(E_control,p_control)
deviation_control <- sqrt(wtd.var(E_control,p_control))

var_modified <- wtd.var(E_modified,p_modified)
deviation_modified <- sqrt(wtd.var(E_modified,p_modified))


## HISTOGRAM

# Preparo datos para histograma y boxplot (creo vector a partir de matriz con frecuencias):
Ev_control <- rep(controlData$kPa, fabs_control$X.)
Ev_modified <- rep(modifiedData$kPa, fabs_modified$X.)
# La longitud de este vector se corresponde con el numero de curvas de fuerza analizadas

# Grafico histograma:
hist(Ev_control, xlim = c(0,max(Ev_control)), breaks = 30, main = "Wild type cells", xlab = "Young's Modulus (kPa)", ylab = "Force curves")
hist(Ev_modified, xlim = c(0,max(Ev_modified)), breaks = 30, main = "Modified cells", xlab = "Young's Modulus (kPa)", ylab = "Force curves")

# Con xlim se puede cambiar el rango del eje x y con breaks la cantidad de barras
# Modificar main segun los tipos de muestras analizadas

# Guardar histograma: Export > Save as Image > Modificar Directory y File name


## BOXPLOT

E = data.frame(modulus = c(Ev_control, Ev_modified), cells = rep(c("Wild type cells", "Modified cells"), times = c(length(Ev_control),length(Ev_modified))))
boxplot(modulus ~ cells, data = E, ylab = "Young's Modulus (kPa)", xlab = "", ylim = c(0,250))

# Con ylim se puede cambiar el rango del eje 
# Modificar cells segun los tipos de muestras analizadas

# Guardar boxplot: Export > Save as Image > Modificar Directory y File name


## TEST DE STUDENT

t.test(Ev_control, Ev_modified)


## RESULTADOS

# Wild type cells:
# E = mean_control +/- deviation_control kPa

# Modified cells:
# E = mean_modified +/- deviation_modified kPa

# Test de Student:
# p-value = 