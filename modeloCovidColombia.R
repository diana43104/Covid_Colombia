# Modelo de transmision del Coronavirus Colombia
# Desarrollado por el Observatorio Nacional de Salud del INstituto Nacional de Salud
# Autores. Carlos Castañeda-Orjuela, Diana Díaz-Jimenez, Liliana Hilarión-Gaitán, Jean Carlo Pineda-Lozano
# Consultar el informe en https://issuu.com/inscomunicaciones/docs/modelo_covid-19_colombia_ins_v5


#############################################
# Instructivo de Uso
#############################################


# 1. Cargue las librerias y el directorio donde van a quedar alojadas las salidas.
# 2. Cargue las funciones definidas de la línea 29 a 261
# 3. Genere los resultados de acuerdo con los parámetros que sean de su interés en el apartado de resultados



# Se carga las librerias y se define la ruta(Setworking - stwd)


# Se cargan las librerías necesarias
library(deSolve)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(readxl)
library(reshape2)


# Se Cargan los parámetros de distribución esenciales del modelo

# Número reproductivo básico (Zhang S, Estimation of the reproductive number of Novel Coronavirus (COVID-19) and the probable outbreak size on the Diamond Princess cruise ship: A data-driven analysis. Int J Infect Dis. 2020;1–9.)
mR0 = 2.28
sdR0 = 0.117349095 # Este parámetro corresponde a la desviación estándar del valor del mR0

# Duración de periodo infeccioso (Linton NM. Incubation Period and Other Epidemiological Characteristics of 2019 Novel Coronavirus Infections with Right Truncation : A Statistical Analysis of Publicly Available Case Data. Clin Med (Northfield Il). 2019;9:1–9)
dur = 5.8

# Parámetros de letalidad (Verity R, et al. Estimates of the severity of COVID-19 disease. medRxiv. 2020;2020.03.09.20033357)
CFRm <- 0.011418095
# CFRl <- 0.010039516
# CFRh <- 0.013018324
CFRsd <- 0.00075991
CFRalfa <- CFRm * ((CFRm * (1 - CFRm) / (CFRsd ^ 2)) - 1)
CFRbeta <- CFRalfa * ((1 - CFRm) / CFRm)



## Se crea una función SIR  (sin considerar mortalidad)
sir <- function(time, state, parameters) {

    with(as.list(c(state, parameters)), {

        dS <- -Beta * S * I
        dI <-  Beta * S * I - Gamma * I
        dR <-                 Gamma * I

    return(list(c(dS, dI, dR)))
  })
}


# Se crea la función de trasmisión COVID 2019
# Se crea una función determinística. 

CovidDet <- function(pob = 38837139, nInitPtes = 5, nRcuperados = 0, R0 = 2.28, nDays = 100, dur = 6){

        initPatients <- nInitPtes

# Se definen los parámetros de las ecuaciones diferenciales
        tRec <-  1 / dur
        tInf <-  R0 * tRec

### Se especifican los parámetros

## Proporciones en cada compartimento segun  Susceptible , Infectado , Recuperado
        init <- c(S = 1 - ((initPatients + nRcuperados) /  pob), I = initPatients /  pob, R = nRcuperados /  pob)
## beta:  parametro de infección; gamma:  parametro de recuperación
        parameters <- c(Beta = tInf, Gamma = tRec)
## Marco temporal
        times      <- seq(0, nDays, by = 1)
## Solucionando ode (General Solver for Ordinary Differential Equations)
        out <- ode(y = init, times = times, func = sir, parms = parameters)
## Cambiando a data frame
        out <- as.data.frame(out)
## Eliminando la variable  tiempo
        out$time <- NULL
## Generando resultado en terminos poblacionales
        out <- out * pob
## Se calcula el total de infectados a cada periodo (infectados activos + recuperados)
        out$totalInf <- out$I + out$R
##  Se reportan los casos estimados (acumulados) para cada día
        estCases <- out$totalInf
        return(estCases)
}

# Se crea una función de reporte de resultados en csv
resumir <- function(objeto, nSim){
      resultados <- mutate(objeto, meanCases = rowMeans(objeto[ , 2:(nSim + 1)]), LCI = rowQuantiles(as.matrix(objeto[ , 2:(nSim + 1)]), probs = 0.025), SCI = rowQuantiles(as.matrix(objeto[ , 2:(nSim + 1)]), probs = 0.975))
      resultados <- resultados[c("day", "meanCases", "LCI", "SCI")]
      resultados$meanCases <- round(resultados$meanCases, 0)
      resultados$LCI <- round(resultados$LCI, 0)
      resultados$SCI <- round(resultados$SCI, 0)
      return(resultados)
}

# Se usa una población fija (pob)
# Se incluye un número incicial de infectados (nInitPtes)
# Se define el R0 (mR0) y su desviación estándar (sdR0)
# El número de días de infección se estima como una distribución uniforme entre 5 y 7 días (minTransm y maxTransm)
# Número de días para la simulación nDays
# Todo los modelos corren con 10 mil iteraciones en su componente estocástico

Covid19Model <- function(name = "Colombia", pob = 38837139, nInitPtes = 5, nRcuperados = 0, mR0 = 2.28, sdR0 = 0.117349095, minTransm = 5, maxTransm = 7, nDays = 100, nSim = 10000){
  
    resultados <- data.frame(day = 0:nDays)
    muertes <- resultados
  # Componente probabilístico
    for(i in 1:nSim){
        set.seed(12345 + i)

        R0 <- rnorm(n = 1, mean = mR0, sd = sdR0) #  Sheng Zhang, 2020
        dur <- round(runif(n = 1, min = minTransm, max = maxTransm), digits = 0) # Se define entre un rango entre 5 y 7 días
        let <- rbeta(1, shape1 = CFRalfa, shape2 = CFRbeta)

        outInfectados <- CovidDet(pob = pob, nInitPtes = nInitPtes, nRcuperados = nRcuperados, R0 = R0, dur = dur, nDays = nDays)
        
        # Se incluye el cálculo de la letalidad (se calculan las muertes sobre el total de infectados, posteriormente se deben ajustar las muertes por el % de asintomáticos)
        outMuertes <- outInfectados * let
        # Se compilan los resultados de infectados y muerte
        resultados[i + 1] <- round(outInfectados, digits = 0)
        muertes[i + 1] <- round(outMuertes, digits = 0)
  }
    # Se compilan los resultados de infectados
    resultados <- resumir(resultados, nSim = nSim)
    # Se imprimen los resulatos de infectados en archivo csv
    write.csv(resultados, paste(name, nInitPtes, "r0", mR0, "Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
    # Se imprime la gráfica de la curva epidémica
    pdf(paste(name, mR0, "Trasmi", minTransm,"-",maxTransm,".pdf", sep = ""))
    p <- ggplot(resultados, aes(x = day, y = meanCases))  + geom_ribbon(aes(ymin = LCI, ymax = SCI), fill = "lightblue") + geom_line(aes()) + xlab("Día") + ylab("Casos")
    print(p)
    dev.off()
    # Se compilan los resultados de infectados
    muertes <- resumir(muertes, nSim = nSim)
    # Se imprimen los resultados de muertes en archivo csv
    write.csv(muertes, paste("Muertes", name, nInitPtes, "r0", mR0, "Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
    # Se retornan las salidas
    return(list(casos = resultados, muertes = muertes))
}



# Modelo con efectividad de diferentes intervenciones

# Función de evaluación de efectividad para periodos intermitentes
Covid19ModelAccIntermitente <- function(name = "Colombia", pob = 38837139, nInitPtes = 5, mR0 = 2.28, sdR0 = 0.117349095, Intervention = "Cuarentena", Efect = 0.69, nDaysNOInt = 28, pInt = 14,  pNOInt = 7, minTransm = 5, maxTransm = 7, nDaysTotal = 300, nSim = 10000){

  # Componente probabilístico
  # modelando la primera parte
  resultados <- data.frame(day = 0:nDaysTotal)
  muertes <- resultados
  
  for(i in 1:nSim){
      set.seed(12345 + i)
      # Se define el número de periodos a evaluar  
      periodos <- ceiling(((nDaysTotal - nDaysNOInt)/(pInt + pNOInt)) * 2 + 1)
    
      daysPeriod <- nDaysNOInt
      R0 <- rnorm(n = 1, mean = mR0, sd = sdR0) #  Sheng Zhang, 2020
      dur <- round(runif(n = 1, min = minTransm, max = maxTransm), digits = 0)
      let <- rbeta(1, shape1 = CFRalfa, shape2 = CFRbeta)
      tRec <-  1 / dur
    
      infectados <- nInitPtes / pob
      susceptibles <- 1 - infectados
      recuperados <- 0
      startT <- 0
      finishT <- nDaysNOInt
      
      cummDays <- nDaysNOInt
      for(p in 1:periodos){
          # Se define si es un periodo de intervención o no
        if(p %% 2 !=  0){
          tInf <-  R0 * tRec
        }else{
          tInf <-  R0 * tRec * (1 - Efect)
        }
      
        if(cummDays > nDaysTotal){
          
          daysPeriod <- daysPeriod - (cummDays - nDaysTotal)
          finishT <- nDaysTotal
          
        }
      ## Proporción  en cada compartimento de acuerdo al número inicial de pacientes: Susceptible , Infectado , Recuperado 

        init <- c(S = susceptibles, I = infectados, R = recuperados)
## beta: parámetro infección; gamma: parámetro de recuperado
        parameters <- c(Beta = tInf, Gamma = tRec)

## Periodo de tiempo
        times <- seq(startT, finishT, by = 1)

## Solucionando ode (General Solver for Ordinary Differential Equations)
        out <- ode(y = init, times = times, func = sir, parms = parameters)
## Cambiando a data frame
        out <- as.data.frame(out)
## Eliminando la variable  tiempo
        out$time <- NULL
## Generando resultado en terminos poblacionales
        out <- out * pob
        out$totalInf <- out$I + out$R
    
        tempOut <- round(out$totalInf, digits = 0)
        lastRow <- out[nrow(out), ] / pob

        susceptibles <- lastRow$S
        infectados <- lastRow$I
        recuperados <- lastRow$R
        
        if(p == 1){
          tempRes <- tempOut
        }else{
          # Se elimina la primera fila
          tempOut <- tempOut[-1]
          tempRes <- c(tempRes, tempOut)
        }
        tempMuertes <- round(tempRes * let)
        # Se ajusta la duración del siguiente periodo
        if(p %% 2 !=  0){
          daysPeriod <- pInt
        }else{
          daysPeriod <- pNOInt
        }
        startT <- finishT
        finishT <- finishT + daysPeriod
        cummDays <- cummDays + daysPeriod
      }
      
      resultados[i + 1] <- tempRes
      muertes[i + 1] <- tempMuertes
    }
  
  #Se compilan los resultados de infectados
    resultados <- resumir(resultados, nSim = nSim)
  # Se imprimen los resulatos de infectados en archivo csv
    write.csv(resultados, paste(name, nInitPtes, "r0", mR0, "-", Intervention, Efect, pInt, "x", pNOInt, "días Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
   # Se imprime la gráfica de la curva epidémica
    pdf(paste(name, mR0,"-",Intervention, Efect, pInt, "x", pNOInt, "días Trasmi", minTransm,"-",maxTransm, ".pdf", sep = ""))
    p <- ggplot(resultados, aes(x = day, y = meanCases))  + geom_ribbon(aes(ymin = LCI, ymax = SCI), fill = "lightblue") + geom_line(aes()) + xlab("Día") + ylab("Casos")
    print(p)
    dev.off()
    
    # Se compilan los resultados de infectados
    muertes <- resumir(muertes, nSim = nSim)
    # Se imprimen los resultados de muertes en archivo csv
    write.csv(muertes, paste("Muertes", name, nInitPtes, "r0", mR0, "-", Intervention, Efect, pInt, "x", pNOInt, "días Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
    # Se retornan las salidas
    return(list(resultados, muertes))
}



  
#######################
#Resultados
#######################



# Los parámetros podrán modificarse de acuerdo a la evidencia cientifica que se ajusta a su contexto local.

# Resultados sin intervención 

# name = el nombre con el queva a quedar guardado el archivo
#pob = Debe ser modificado por la población urbana de cada entidad territorial para el caso de colombia se consideró que el total de la población urbana era el 77,1% 
#mR0 = Número reproductivo básico
#nInitPtes = numero inicial de infectados para colombia es de 5 para las entidades territoriales de 1.
#sdR0 = desviación estandar del R0 
#minTransm = Periodo mínimo de transmisión
#maxTransm = Periodo máximo de transmisión
#NDays = número de dias que corre el modelo


# Se genera los resultados sin ninguna intervención, como fue descrito previamente puede ser modificados los parámentros
Covid19Model(name = "Colombia_prueba", pob = 38837139, nInitPtes = 5, mR0 = 2.28, sdR0 = 0.117349095, minTransm = 5, maxTransm = 7,  nDays = 100)


# Resultados con intervención 

# así como en la función anterior, se puede modificar los parametros mencionado y los siguientes

# Intervention = la intervención que se está evaluando
# Efect = la efectividad de la intervención 
# nDaysNOInt = número de días sin intervención 
# pInt = periodo con internveción
# pNOInt = periodo sin intervención
# nDaysTotal = total de días simulados de la epidemia

# prueba 1. Se generan resultados con la intervención "Cierre escuelas + aislamiento + cuarentena + aislamiento social70" con una efectividad del 47%, 
Covid19ModelAccIntermitente(name = "Colombia", pob = 38837139, nInitPtes = 5, Intervention = "Cierre escuelas + aislamiento + cuarentena + aislamiento social70", Efect = 0.47, nDaysNOInt = 18, pInt = 33,  pNOInt = 249, nDaysTotal = 300)

# prueba 2. Se generan resultados con la intervención "Cuarentena 25 marzo 11 mayo abrir R 1,0" con una efectividad del 56%, 
Covid19ModelAccIntermitente(name = "Colombia", pob = 38837139, nInitPtes = 5, Intervention = "Cuarentena 25 marzo 11 mayo abrir R 1,0", Efect = 0.56, nDaysNOInt = 25, pInt = 47,  pNOInt = 228, nDaysTotal = 300)




