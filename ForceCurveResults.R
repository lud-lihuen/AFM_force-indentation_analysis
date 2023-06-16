# Procesar resultados del analisis de curvas de fuerza (nanoindentacion)


## Importar vectores ascii de resultados con Import Dataset > From Text (base) > Import


## CONTROL GROUP DATA

E_control = archivo$V1 #kPa (Reemplazar "archivo" por nombre del archivo de modulo de Young)
E_control_mean <- mean(E_control) #media
E_control_deviation <- sd(E_control) #standard deviation

A_control = archivo$V1 #kPa (Reemplazar "archivo" por nombre del archivo de trabajo de adhesion)
A_control_mean <- mean(E_control) #media
A_control_deviation <- sd(E_control) #standard deviation

## MODIFIED GROUP DATA

E_modified = archivo$V1 #kPa (Reemplazar "archivo" por nombre del archivo de modulo de Young)
E_modified_mean <- mean(E_modified) #media
E_modified_deviation <- sd(E_modified) #standard deviation

A_modified = archivo$V1 #kPa (Reemplazar "archivo" por nombre del archivo de trabajo de adhesion)
A_modified_mean <- mean(E_modified) #media
A_modified_deviation <- sd(E_modified) #standard deviation


## HISTOGRAM

hist(E_control, xlim = c(0,max(E_control)), breaks = 30, main = "Wild type cells", xlab = "Young's Modulus (kPa)", ylab = "Force curves")
hist(E_modified, xlim = c(0,max(E_modified)), breaks = 30, main = "Modified cells", xlab = "Young's Modulus (kPa)", ylab = "Force curves")

hist(A_control, xlim = c(0,max(A_control)), breaks = 30, main = "Wild type cells", xlab = "Work of Adhesion (J/m2)", ylab = "Force curves")
hist(A_modified, xlim = c(0,max(A_modified)), breaks = 30, main = "Modified cells", xlab = "Work of Adhesion (J/m2)", ylab = "Force curves")

# Con xlim se puede cambiar el rango del eje x y con breaks la cantidad de barras
# Modificar main segun los tipos de muestras analizadas

# Guardar histograma: Export > Save as Image > Modificar Directory y File name


## BOXPLOT

E = data.frame(modulus = c(E_control, E_modified), cells = rep(c("Wild type cells", "Modified cells"), times = c(length(E_control),length(E_modified))))
A = data.frame(adhesion = c(A_control, A_modified), cells = rep(c("Wild type cells", "Modified cells"), times = c(length(A_control),length(A_modified))))

boxplot(modulus ~ cells, data = E, ylab = "Young's Modulus (kPa)", xlab = "", ylim = c(0,250))
boxplot(adhesion ~ cells, data = A, ylab = "Work of Adhesion (J/m2)", xlab = "", ylim = c(0,250))

# Con ylim se puede cambiar el rango del eje 
# Modificar cells segun los tipos de muestras analizadas

# Guardar boxplot: Export > Save as Image > Modificar Directory y File name


## TEST DE STUDENT

t.test(E_control, E_modified) # Para modulo de Young
t.test(A_control, A_modified) # Para trabajo de adhesion


## RESULTADOS

# Wild type cells:
# E = E_control_mean +/- E_control_deviation kPa
# A = A_control_mean +/- A_control_deviation J/m2

# Modified cells:
# E = E_modified_mean +/- E_modified_deviation kPa
# A = A_modified_mean +/- A_modified_deviation J/m2

# Test de Student:
# p-value para E = 
# p-value para A = 