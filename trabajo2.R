## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
set.seed(77)

## ------------------------------------------------------------------------
E = function(u,v){(u^2*exp(v)-2*v^2*exp(-u))^2}
duE = function(u,v){2*(u^2*exp(v)-2*v^2*exp(-u))*(2*exp(v)*u+2*v^2*exp(-u))}
dvE = function(u,v){2*(u^2*exp(v)-2*v^2*exp(-u))*(u^2*exp(v)-4*v*exp(-u))}

## ------------------------------------------------------------------------
# Función que calcula el gradiente descendente de una función dada,
# junto a sus derivadas previamente calculadas. El algoritmo actualiza
# los pesos utilizando el valor obtenido del gradiente y la tasa de
# aprendizaje que se quiera utilizar en cada caso. En el momento en el
# que la diferencia de valores en la función original para pesos(t) y
# pesos(t-1) sea menor que el umbral propuesto, el algoritmo se detendrá.
# También se detendrá si se alcanza el máximo de iteraciones establecido.
# Para facilitar la resolución de los próximos ejercicios, se devuelve
# una lista con los pesos que dan lugar al mínimo calculado y las
# iteraciones necesarias para concurrir hasta él.
gradienteDescendente = function(wini, mu, threshold, max_iter, E, duE, dvE){
  
  wold = wini
  iter = 0
  seguir_iterando = T
  
  while (iter < max_iter & seguir_iterando){
    g=c(duE(wold[1],wold[2]),dvE(wold[1],wold[2]))
    v = -g
    wnew = wold + mu*v
    if (abs(E(wold[1],wold[2]) - E(wnew[1],wnew[2])) < threshold){
      seguir_iterando = F
    }
    wold = wnew
    iter = iter+1  
  }
  list(minimo=wnew, iter=iter)
}

resultado1a = gradienteDescendente(c(1,1), 0.1, 10^-4, 50, E, duE, dvE)
resultado1a$iter


## ------------------------------------------------------------------------
resultado1a$minimo

## ------------------------------------------------------------------------
fun1b = function(x,y){(x-2)^2 + 2*(y-2)^2 + 2*sin(2*pi*x)*sin(2*pi*y)}
dxfun1b = function(x,y){2*(-2+x+2*pi*cos(2*pi*x)*sin(2*pi*y))}
dyfun1b = function(x,y){4*(-2+y+pi*cos(2*pi*y)*sin(2*pi*x))}

# Función equivalente al gradiente descendente pero con la particularidad
# de que dibuja en una gráfica los puntos "mínimos" que va encontrando el
# algoritmo en función del número de iteraciones que se han evaluado.
gradienteDescendenteGrafica = function(wini, mu, threshold, max_iter, E, duE, dvE){
  
  wold = wini
  iter = 0
  seguir_iterando = T
  iteraciones_grafica = list()
  minimos_grafica = list()

  while (iter < max_iter & seguir_iterando){
    g=c(duE(wold[1],wold[2]),dvE(wold[1],wold[2]))
    v = -g
    wnew = wold + mu*v
    iteraciones_grafica = c(iteraciones_grafica, iter)
    minimos_grafica = c(minimos_grafica, E(wnew[1],wnew[2]))
    if (abs(E(wold[1],wold[2]) - E(wnew[1],wnew[2])) < threshold){
      seguir_iterando = F
    }
    wold = wnew
    iter = iter+1  
  }
  plot(iteraciones_grafica, minimos_grafica)
  list(minimo=wnew, iter=iter)
}

par(mfrow = c(1,2))
resultado1b = gradienteDescendenteGrafica(c(1,1), 0.01, 10^-4, 50, fun1b, dxfun1b, dyfun1b)
resultado1b$minimo
resultado1b$iter

resultado1b2 = gradienteDescendenteGrafica(c(1,1), 0.1, 10^-4, 50, fun1b, dxfun1b, dyfun1b)
resultado1b2$minimo
resultado1b2$iter

## ------------------------------------------------------------------------
valores_punto1 = gradienteDescendente(c(2.1,2.1), 0.1, 10^-4, 50, fun1b, dxfun1b, dyfun1b)
valores_punto2 = gradienteDescendente(c(3,3), 0.1, 10^-4, 50, fun1b, dxfun1b, dyfun1b)
valores_punto3 = gradienteDescendente(c(1.5,1.5), 0.1, 10^-4, 50, fun1b, dxfun1b, dyfun1b)
valores_punto4 = gradienteDescendente(c(1,1), 0.1, 10^-4, 50, fun1b, dxfun1b, dyfun1b)

valores_punto1$minimo
print(fun1b(valores_punto1$minimo[1], valores_punto1$minimo[2]))
valores_punto2$minimo
print(fun1b(valores_punto2$minimo[1], valores_punto2$minimo[2]))
valores_punto3$minimo
print(fun1b(valores_punto3$minimo[1], valores_punto3$minimo[2]))
valores_punto4$minimo
print(fun1b(valores_punto4$minimo[1], valores_punto4$minimo[2]))


## ------------------------------------------------------------------------
# Funciones recicladas de la práctica 1, para generar datos aleatorios uniformemente
# distribuidos, rectas, y asignar etiquetas. También se incluye la función que
# calcula una recta en función de unos pesos, que se usará más adelante.
simula_unif = function (N=2,dims=2, rango = c(0,1)){
 m = matrix(runif(N*dims, min=rango[1], max=rango[2]),
 nrow = N, ncol=dims, byrow=T)
 m
}

simula_recta = function (intervalo = c(-1,1), visible=F){

  ptos = simula_unif(2,2,intervalo) # se generan 2 puntos
   a = (ptos[1,2] - ptos[2,2]) / (ptos[1,1]-ptos[2,1]) # calculo de la pendiente
   b = ptos[1,2]-a*ptos[1,1]  # calculo del punto de corte

   if (visible) {  # pinta la recta y los 2 puntos
       if (dev.cur()==1) # no esta abierto el dispositivo lo abre con plot
           plot(1, type="n", xlim=intervalo, ylim=intervalo)
       points(ptos,col=3)  #pinta en verde los puntos
       abline(b,a,col=3)   # y la recta
   }
   c(a,b) # devuelve el par pendiente y punto de corte
}

asignarEtiquetasSegunRecta = function(recta,puntos){
	etiquetas = sign(puntos[,2] - (recta[1]*puntos[,1]+recta[2]))
}

obtenerRectaPesos = function(pesos){
  c(-pesos[1]/pesos[2], -pesos[3]/pesos[2])
}

datosej2 = simula_unif(100, 2, c(0,2))
rectaej2 = simula_recta(c(0,2),F)
etiquetasej2 = asignarEtiquetasSegunRecta(rectaej2, datosej2)

plot(datosej2, col=etiquetasej2+3)
abline(rectaej2[2], rectaej2[1])

## ------------------------------------------------------------------------
# Función que computa el gradiente descendente estocástico para un punto con su
# etiqueta y un vector de pesos. Esta función será llamada dentro del algoritmo
# de regresión logística cuando se deban ir actualizando los pesos.
GDEstocastico = function(punto, etiqueta, w){
  -(crossprod(etiqueta,punto)[1,] / (1 + exp(etiqueta*w%*%punto)))
}

# Función que calcula el módulo de un vector como la raíz de la sumatoria de sus
# componentes al cuadrado.
moduloVector = function(vector){
  sqrt(sum(vector^2))
}

# Función que calcula unos pesos que separan unos datos de unas etiquetas por
# regresión logística. El algoritmo itera el número de épocas que hagan falta,
# con los pesos actualizándose con cada dato y comprobándose al final de la misma
# si el módulo del nuevo vector de pesos varía más del umbral con respecto al
# de la época anterior. Este algoritmo ha sido extraído de la página 95 del libro
# "Learning from data" de Abu-Mostafa et al.
regresionLogistica = function(datos, etiquetas, threshold, mu){
  wnew = c(0,0,0)
  wold = c(0,0,0)
  datos = cbind(datos,1)
  epoca = 0
  modulo = threshold*2
  evaluado = 1
  barajados = sample(1:dim(datosej2)[1])

  while (modulo > threshold){
    wold = wnew
    
    while (evaluado <= length(barajados)){
      g = GDEstocastico(datos[barajados[evaluado],], etiquetas[barajados[evaluado]], wold) 
      evaluado = evaluado+1
      v = -g
      wnew = wnew + mu*v
    }
    
    evaluado = 1
    barajados = sample(1:dim(datosej2)[1])
    
    modulo = moduloVector(wnew - wold)
    epoca = epoca + 1
  }

  resultados = list(pesos=wnew, epocas=epoca)
}

## ------------------------------------------------------------------------
# Función que calcula el error para la regresión logística basándose en la fórmula
# de la página 98 del libro "Learning from data", a partir de los datos, las etiquetas
# y los pesos que se han generado. Finalmente devuelve la media del error de todos
# los datos.
errorRLog = function(datos,etiquetas,pesos){
  error = 0
  datos = cbind(datos,1)
  for (i in 1:dim(datos)[1]){
    error = error + log(1 + exp(-etiquetas[i]*(pesos%*%datos[i,])), exp(1))
  }  
  error/dim(datos)[1]
}

# Función que computa el experimento pedido en el ejercicio 2b, generando N muestras
# y manteniendo los pesos calculados con regresión logística y la recta separadora
# original para comprobar el error fuera de la muestra.
ejercicio2b = function(recta, pesos, muestras){
  datos = simula_unif(muestras, 2, c(0,2))
  etiquetas = asignarEtiquetasSegunRecta(recta, datos)
  
  Eout = errorRLog(datos, etiquetas, pesos)
  Eout
}

resultadoRLog = regresionLogistica(datosej2, etiquetasej2, 0.01, 0.01)
rectapesosRLog = obtenerRectaPesos(resultadoRLog$pesos)
resultadoRLog$epocas
resultadoRLog$pesos

par(mfrow=c(1,2))
plot(datosej2, col=etiquetasej2+3)
abline(rectaej2[2], rectaej2[1])
plot(datosej2, col=etiquetasej2+3)
abline(rectapesosRLog[2], rectapesosRLog[1])

Eout = ejercicio2b(rectapesosRLog, resultadoRLog$pesos, 1000)
Ein = errorRLog(datosej2, etiquetasej2, resultadoRLog$pesos)
Eout
Ein

## ------------------------------------------------------------------------
digit.train <- read.table("datos/zip.train",
                          quote="\"", comment.char="", stringsAsFactors=FALSE)

digitos48.train = digit.train[digit.train$V1==4 | digit.train$V1==8,]
digitos = digitos48.train[,1]  # etiquetas
ndigitos = nrow(digitos48.train)

# se retira la clase y se monta una matriz 3D: 599*16*16
grises = array(unlist(subset(digitos48.train,select=-V1)),c(ndigitos,16,16))
rm(digit.train) 
rm(digitos48.train)

# Para visualizar los 4 primeros
## ------------------------------------------------------------------------

par(mfrow=c(2,2)) 
for(i in 1:4){
  imagen = grises[i,,16:1] # se rota para verlo bien
  image(z=imagen)
}

digitos[1:4] # etiquetas correspondientes a las 4 imágenes

fsimetria <- function(A){
  A = abs(A-A[,ncol(A):1])
  -mean(A)
}

simetria = apply(grises, 1, fsimetria)
intensidadPromedio = apply(grises, 1, mean)
par(mfrow = c(1,2))
plot(x=intensidadPromedio, y=simetria, col=digitos+1)

digit.test <- read.table("datos/zip.test",
                          quote="\"", comment.char="", stringsAsFactors=FALSE)

digitos48.test = digit.test[digit.test$V1==4 | digit.test$V1==8,]
digitos.test = digitos48.test[,1]  # etiquetas
ndigitos.test = nrow(digitos48.test)

# se retira la clase y se monta una matriz 3D: 599*16*16
grises.test = array(unlist(subset(digitos48.test,select=-V1)),c(ndigitos.test,16,16))
rm(digit.test) 
rm(digitos48.test)

simetria.test = apply(grises.test, 1, fsimetria)
intensidadPromedio.test = apply(grises.test, 1, mean)
plot(x=intensidadPromedio.test, y=simetria.test, col=digitos.test+1)

## ------------------------------------------------------------------------
# Función que calcula el error de unas etiquetas de unos datos para un vector de pesos
# en base a cuántas etiquetas están mal situadas con respecto a la recta que formaría
# dicho vector.
errores = function(datos,label,pesos){
  (sum(sign(datos%*%pesos) != label))/length(label)
}

# Función que calcula los pesos por regresión lineal de unos datos y unas etiquetas
# aportadas. Se hace uso de la fórmula para calcular la matriz pseudoinversa de
# los datos que aparece en las diapositivas de teoría. Según esto, la matriz
# pseudoinversa de X se puede obtener como la multiplicación de la inversa de
# la multiplicación de X traspuesta por X, por la traspuesta de X. En fórmula
# matemática, pseudoX = (XT*X)-1 * XT.
# El resto de operaciones son las mismas que aparecen en las transparencias,
# empleando la función svd para extraer las matrices U, D y V de X, y sabiendo
# que (XT*X)-1 = V * (pseudoD)^2 * VT (en las transparencias no aparece el cuadrado
# por error, pero lo corregimos en clase) y que la pseudoD es la diagonal de 1/D.
# Finalmente, los pesos se obtienen de multiplicar la pseudoinversa de X por las
# etiquetas.
Regress_Lin = function(datos,label){
  
  datos=cbind(datos,1)

  descomposicion = svd(datos)
  pseudoinversaD = diag(1/descomposicion$d)
  inversadatosTdatos = (descomposicion$v)%*%(pseudoinversaD^2) %*%t(descomposicion$v)
  pseudoinversa = inversadatosTdatos%*%(t(datos))
  
  pesos = (pseudoinversa%*%label)
}

# Función que ajusta el algoritmo Perceptron para unos datos y unas etiquetas.
# Funciona similarmente al PLA original pero este guarda la mejor solución
# por la que ha pasado (la que arrojaba el error más pequeño) y sólo la actualiza
# cuando los nuevos pesos generados en cada iteración dan un mejor resultado.
# Devuelve en este caso el mejor resultado y no el último.
PLA_pocket = function(datos, label, max_iter, vini){
  w = vini
  cambio = T
  datos = cbind(datos,1)
  iteraciones = 0
  
  menor_error = errores(datos,label,w)
  mejor_w = w
  
  while(iteraciones < max_iter & cambio){
    cambio = F
    for(j in sample(1:dim(datos)[1])){
      if (sign(crossprod(datos[j,],w)) != label[j]){
        w = w + datos[j,]*label[j]
        cambio = T
        break
      }
    }
    
    if (cambio){
      error_actual = errores(datos,label,w)
      if (error_actual < menor_error){
        menor_error = error_actual
        mejor_w = w
      }
    }
    
    iteraciones = iteraciones+1
  }
  
  list(w=mejor_w, iter=iteraciones)
}

## ------------------------------------------------------------------------
datos48train = matrix(c(intensidadPromedio,simetria),nrow=length(intensidadPromedio),ncol=2)
digitos[digitos==4]=1
digitos[digitos==8]=-1

datos48test = matrix(c(intensidadPromedio.test,simetria.test),nrow=length(intensidadPromedio.test),ncol=2)
digitos.test[digitos.test==4]=1
digitos.test[digitos.test==8]=-1

pesos_train = Regress_Lin(datos48train, digitos)
resPLApocket = PLA_pocket(datos48train, digitos, 10000, pesos_train)
pesos_train = resPLApocket$w
recta_pesos_train = obtenerRectaPesos(pesos_train)

par(mfrow=c(1,2))
plot(datos48train, col=digitos+3)
abline(recta_pesos_train[2], recta_pesos_train[1])
plot(datos48test, col=digitos+3)
abline(recta_pesos_train[2], recta_pesos_train[1])

## ------------------------------------------------------------------------
Ein_3 = errores(cbind(datos48train, 1), digitos, pesos_train)
Etest_3 = errores(cbind(datos48test, 1), digitos.test, pesos_train)
Ein_3
Etest_3

## ------------------------------------------------------------------------
# Función que calcula la cota de error estimada en base a la fórmula de la
# cota de generalización Vapnik-Chervonenkis. La tolerancia determinará el
# porcentaje de confianza que podemos esperar para estos resultados.
# La dimensión de Vapnik-Chervonenkis del Perceptron 2D es 3.
calcularCotaRespectoError = function(error, dvc, tolerancia, N){
  error + sqrt(8/N * log((4*((2*N)^dvc + 1)/tolerancia), exp(1)))
}

cota_segun_ein = calcularCotaRespectoError(Ein_3, 3, 0.05, dim(datos48train)[1])
cota_segun_etest = calcularCotaRespectoError(Etest_3, 3, 0.05, dim(datos48test)[1])
cota_segun_ein
cota_segun_etest

## ------------------------------------------------------------------------
# Función aportada en la práctica 1, con una ligera modificación para que acepte
# la media que se le quiere dar a la distribución gaussiana.
simula_gaus = function(N=2,dim=2,sigma, mean){

  if (missing(sigma)) stop("Debe dar un vector de varianzas")
  sigma = sqrt(sigma)  # para la generación se usa sd, y no la varianza
  if(dim != length(sigma)) stop ("El numero de varianzas es distinto de la dimensión")

  simula_gauss1 = function() rnorm(dim, mean=mean, sd = sigma) # genera 1 muestra, con las desviaciones especificadas
  m = t(replicate(N,simula_gauss1())) # repite N veces, simula_gauss1 y se hace la traspuesta
  m
}

# Función que calcula las etiquetas de los datos de la gaussiana en función de los pesos wf,
# y metiendo un ruido como se indica en el enunciado.
obtenerEtiquetasGauss = function(datos, pesos, sigma){
  etiquetas = as.vector(as.vector(datos%*%t(pesos))+sigma*simula_gaus(dim(datos)[1],1,sigma,0))
}

datosgauss3 = simula_gaus(130, 3, c(1,1,1), 1)

pesosgauss3 = simula_gaus(4,1,1,0)
pesosgauss3[4] = pesosgauss3[4]+1

etiquetasgauss3 = obtenerEtiquetasGauss(cbind(datosgauss3,1), pesosgauss3, 0.5)

## ------------------------------------------------------------------------
# Función que calcula los pesos según la fórmula de weight decay como se indica
# en la sesión 7 de teoría, utilizando el parámetro de regularización que se
# especifica en el enunciado y empleando la técnica de cálculo de pseudoinversas
# con descomposición SVD.
weightDecay = function(datos, etiquetas, regularizacion){
  tam = dim(datos)[2]
  matriz = t(datos)%*%datos + regularizacion*diag(tam)

  descomposicion = svd(matriz)
  pseudoinversaD = diag(1/descomposicion$d)
  inversadatosTdatos = (descomposicion$v)%*%(pseudoinversaD^2) %*%t(descomposicion$v)
  pseudoinversa = inversadatosTdatos%*%(t(matriz))
  
  pseudoinversa %*% t(datos)%*%etiquetas
}

wreg = weightDecay(cbind(datosgauss3,1), etiquetasgauss3, 0.05/dim(datosgauss3)[1])
wreg

## ------------------------------------------------------------------------
# Función que calcula el error de mínimos cuadrados como aparece en la página
# 139 del libro "Learning from data", evaluando la diferencia cuadrada entre
# las etiquetas reales y las obtenidas con los pesos del weight decay.
errorMinimosCuadrados = function(test, etiquetas_test, pesos){
  errores = (as.vector(test%*%as.vector(pesos)) - etiquetas_test)^2
  sum(errores)
}

# Función que realiza los cálculos pertinentes al experimento propuesto
# y que devuelve una matriz de errores, donde la columna i equivale al
# error de validación del conjunto i, y la fila j contiene los errores
# del experimento j. La última columna contiene la media de los errores.
experimento4 = function(datos, etiquetas){
  d = dim(datos)[2]
  datos = cbind(datos,1)
  errores = matrix(ncol=11)
  errores = errores[-1,]
  for (experimento in 1:1000){
    datos_barajados = sample(1:dim(datos)[1])
    errores_experimento = vector()
    for (indice in 1:10){
      intervalo_inf = (indice-1)*d + (indice-1)*10 + 1
      intervalo_sup = indice*d + indice*10
      train = datos[datos_barajados[-(intervalo_inf:intervalo_sup)],]
      test = datos[datos_barajados[intervalo_inf:intervalo_sup],]
      etiquetas_train = etiquetas[datos_barajados[-(intervalo_inf:intervalo_sup)]]
      etiquetas_test = etiquetas[datos_barajados[intervalo_inf:intervalo_sup]]
      pesos_decay = weightDecay(train, etiquetas_train, 0.05/dim(train)[1])
      
      error = errorMinimosCuadrados(test, etiquetas_test, pesos_decay)
      errores_experimento = c(errores_experimento, error)
    }
    error = sum(errores_experimento)/10
    errores_experimento = c(errores_experimento, error)
    
    errores = rbind(errores, errores_experimento)
  }  
  errores
}

matrizerrores = experimento4(datosgauss3, etiquetasgauss3)

print("Errores")
e1 = mean(matrizerrores[,1])
e1
e2 = mean(matrizerrores[,2])
e2
e3 = mean(matrizerrores[,3])
e3
e4 = mean(matrizerrores[,4])
e4
e5 = mean(matrizerrores[,5])
e5
e6 = mean(matrizerrores[,6])
e6
e7 = mean(matrizerrores[,7])
e7
e8 = mean(matrizerrores[,8])
e8
e9 = mean(matrizerrores[,9])
e9
e10 = mean(matrizerrores[,10])
e10
ecv = mean(matrizerrores[,11])
ecv

print("Varianzas")
v1 = var(matrizerrores[,1])
v1
v2 = var(matrizerrores[,2])
v2
v3 = var(matrizerrores[,3])
v3
v4 = var(matrizerrores[,4])
v4
v5 = var(matrizerrores[,5])
v5
v6 = var(matrizerrores[,6])
v6
v7 = var(matrizerrores[,7])
v7
v8 = var(matrizerrores[,8])
v8
v9 = var(matrizerrores[,9])
v9
v10 = var(matrizerrores[,10])
v10
vcv = mean(c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10))
vcv

## ------------------------------------------------------------------------
# Función que calcula el algoritmo de coordenada descendente de forma
# similar al algoritmo de gradiente descendente pero realizando el gradiente
# explorando primero en la coordenada x y a continuación en la coordenada y
# con la x calculada.
coordenadaDescendente = function(wini, mu, threshold, max_iter, E, duE, dvE){
  
  wold = wini
  iter = 0
  seguir_iterando = T
  wnew = vector(length=2)
  
  while (iter < max_iter & seguir_iterando){
    g=c(duE(wold[1],wold[2]),dvE(wold[1],wold[2]))
    v = -g
    wnew[1] = wold[1] + mu*v[1]
    g=c(duE(wold[1],wold[2]),dvE(wnew[1],wold[2]))
    v = -g
    wnew[2] = wold[2] + mu*v[2]
    if (abs(E(wold[1],wold[2]) - E(wnew[1],wnew[2])) < threshold){
      seguir_iterando = F
    }
    wold = wnew
    iter = iter+1  
  }
  list(minimo=wnew, iter=iter)
}

resultado_coordesc = coordenadaDescendente(c(1,1), 0.1, 10^-4, 50, E, duE, dvE)
resultado_coordesc$minimo
E(resultado_coordesc$minimo[1], resultado_coordesc$minimo[2])

## ------------------------------------------------------------------------
# Función que realiza los cálculos propuestos en el bonus 3, almacenando
# en una lista de salida el error medio y la media de épocas que han tardado
# los experimentos.
bonus3 = function(num_veces){
  errores = vector()
  epocas = vector()
  for (i in 1:num_veces){
    datos = simula_unif(100, 2, c(0,2))
    recta = simula_recta(c(0,2),F)
    etiquetas = asignarEtiquetasSegunRecta(recta, datos)
    
    lista_resultados = regresionLogistica(datos, etiquetas, 0.01, 0.01)
    error = errorRLog(datos, etiquetas, lista_resultados$pesos)
    errores = c(errores, error)
    epocas = c(epocas, lista_resultados$epocas)
  }
  list(epocas=mean(epocas), eout=mean(errores))
}

resultadosbonus3 = bonus3(100)
resultadosbonus3$eout

## ------------------------------------------------------------------------
resultadosbonus3$epocas

