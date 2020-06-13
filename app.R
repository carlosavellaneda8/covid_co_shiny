###########################################################################
## Autores: Angélica Moreno, Lina Buitrago y Carlos Avellaneda
## Descripción: Implementación del modelo SIR en R y Shiny
###########################################################################
# TODO:
# - Agregar el resto de curvas del modelo SIR en escenarios
###########################################################################

rm(list=ls())
library(shiny)
library(shinydashboard)
library(plotly)
library(deSolve)
library(dplyr)
library(matrixStats)
library(reshape2)

###########################################################################
# Global ------------------------------------------------------------------
###########################################################################
## Parametros: ####

# Numero reproductivo basico (Zhang S, Estimation of the reproductive number of Novel Coronavirus (COVID-19) and the probable outbreak size on the Diamond Princess cruise ship: A data-driven analysis. Int J Infect Dis. 2020;1–9.)
# mR0 = 2.28
sdR0 = 0.117349095 # Este parametro corresponde a la desviacion estandar del valor del mR0

# Duracion de periodo infeccioso (Linton NM. Incubation Period and Other Epidemiological Characteristics of 2019 Novel Coronavirus Infections with Right Truncation: A Statistical Analysis of Publicly Available Case Data. Clin Med (Northfield Il). 2019;9:1–9)
dur = 5.8

# Parametros de letalidad (Verity R, et al. Estimates of the severity of COVID-19 disease. medRxiv. 2020;2020.03.09.20033357)
CFRm <- 0.011418095
# CFRl <- 0.010039516
# CFRh <- 0.013018324
CFRsd <- 0.00075991

## Funciones: ####

## Se crea una funcion SIR  (sin considerar mortalidad)
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -Beta * S * I
    dI <-  Beta * S * I - Gamma * I
    dR <-                 Gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}


# Se crea la funcion de trasmision COVID 2019
# Se crea una funcion deterministica. 

CovidDet <- function(pob = 38837139, nInitPtes = 5, nRcuperados = 0, R0 = 2.28, nDays = 100, dur = 6){
  
  initPatients <- nInitPtes
  
  # Se definen los parametros de las ecuaciones diferenciales
  tRec <-  1 / dur
  tInf <-  R0 * tRec
  
  ### Se especifican los parametros
  
  ## Proporciones en cada compartimento segun  Susceptible , Infectado , Recuperado
  init <- c(S = 1 - ((initPatients + nRcuperados) /  pob), I = initPatients /  pob, R = nRcuperados /  pob)
  ## beta:  parametro de infeccion; gamma:  parametro de recuperacion
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
  ##  Se reportan los casos estimados (acumulados) para cada dia
  # estCases <- out$totalInf
  return(out)
}

# Se crea una funcion de reporte de resultados en csv
resumir <- function(objeto, nSim){
  resultados <- mutate(objeto, meanCases = rowMeans(objeto[ , 2:(nSim + 1)]), LCI = rowQuantiles(as.matrix(objeto[ , 2:(nSim + 1)]), probs = 0.025), SCI = rowQuantiles(as.matrix(objeto[ , 2:(nSim + 1)]), probs = 0.975))
  resultados <- resultados[c("day", "meanCases", "LCI", "SCI")]
  resultados$meanCases <- round(resultados$meanCases, 0)
  resultados$LCI <- round(resultados$LCI, 0)
  resultados$SCI <- round(resultados$SCI, 0)
  return(resultados)
}

# Se usa una poblacion fija (pob)
# Se incluye un numero incicial de infectados (nInitPtes)
# Se define el R0 (mR0) y su desviacion estandar (sdR0)
# El numero de dias de infeccion se estima como una distribucion uniforme entre 5 y 7 dias (minTransm y maxTransm)
# Numero de dias para la simulacion nDays
# Todo los modelos corren con 10 mil iteraciones en su componente estocastico

Covid19Model <- function(name = "Colombia", pob = 38837139, nInitPtes = 5, nRcuperados = 0, mR0 = 2.28, sdR0 = 0.117349095, minTransm = 5, maxTransm = 7, nDays = 100, nSim = 10000){
  
  susceptibles <- data.frame(day = 0:nDays)
  infectados <- susceptibles
  recuperados <- susceptibles
  muertes <- susceptibles
  
  # Componente probabilistico
  for(i in 1:nSim){
    set.seed(12345 + i)
    
    R0 <- rnorm(n = 1, mean = mR0, sd = sdR0) #  Sheng Zhang, 2020
    dur <- round(runif(n = 1, min = minTransm, max = maxTransm), digits = 0) # Se define entre un rango entre 5 y 7 dias
    let <- rbeta(1, shape1 = CFRalfa, shape2 = CFRbeta)
    
    outInfectados <- CovidDet(pob = pob, nInitPtes = nInitPtes, nRcuperados = nRcuperados, R0 = R0, dur = dur, nDays = nDays)
    
    # Se incluye el calculo de la letalidad (se calculan las muertes sobre el total de infectados, posteriormente se deben ajustar las muertes por el % de asintomaticos)
    outMuertes <- outInfectados$R * let
    # Se compilan los resultados de infectados y muerte
    susceptibles[i + 1] <- round(outInfectados$S, digits = 0)
    infectados[i + 1] <- round(outInfectados$I, digits = 0)
    recuperados[i + 1] <- round(outInfectados$R - outMuertes, digits = 0)
    muertes[i + 1] <- round(outMuertes, digits = 0)
  }
  # Se compilan los resultados de infectados
  susceptibles <- resumir(susceptibles, nSim = nSim)
  infectados <- resumir(infectados, nSim = nSim)
  recuperados <- resumir(recuperados, nSim = nSim)
  muertes <- resumir(muertes, nSim = nSim)
  # Se imprimen los resultados de muertes en archivo csv
  # write.csv(muertes, paste("Muertes", name, nInitPtes, "r0", mR0, "Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
  # Se retornan las salidas
  return(list(susceptibles = susceptibles,
              infectados = infectados,
              recuperados = recuperados,
              muertes = muertes))
}



# Modelo con efectividad de diferentes intervenciones

# Funcion de evaluacion de efectividad para periodos intermitentes
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
  # write.csv(resultados, paste(name, nInitPtes, "r0", mR0, "-", Intervention, Efect, pInt, "x", pNOInt, "días Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
  # Se imprime la gráfica de la curva epidémica
  # pdf(paste(name, mR0,"-",Intervention, Efect, pInt, "x", pNOInt, "días Trasmi", minTransm,"-",maxTransm, ".pdf", sep = ""))
  # p <- ggplot(resultados, aes(x = day, y = meanCases))  + geom_ribbon(aes(ymin = LCI, ymax = SCI), fill = "lightblue") + geom_line(aes()) + xlab("Día") + ylab("Casos")
  # print(p)
  # dev.off()
  
  # Se compilan los resultados de infectados
  muertes <- resumir(muertes, nSim = nSim)
  # Se imprimen los resultados de muertes en archivo csv
  # write.csv(muertes, paste("Muertes", name, nInitPtes, "r0", mR0, "-", Intervention, Efect, pInt, "x", pNOInt, "días Trasmi", minTransm,"-",maxTransm, ".csv", sep = ""), row.names = FALSE)
  # Se retornan las salidas
  return(list(resultados, muertes))
}


###########################################################################
# UI ----------------------------------------------------------------------
###########################################################################
ui <- dashboardPage(
  dashboardHeader(title = "Modelo SIR"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Escenario base", tabName = "base"),
      menuItem("Comparación de escenarios", tabName = "comparison")
    )
  ),
  dashboardBody(
    tabItems(

      # UI Modelo base -------------------------------------------------------------

      tabItem(
        tabName = "base",
        fluidPage(
          titlePanel("Modelo SIR COVID-19"),
          sidebarLayout(
            sidebarPanel(
              selectInput(
                "city",
                "Seleccione la unidad geográfica:",
                c("Colombia", "Bogotá"),
                selected = "Colombia"
              ),
              sliderInput(
                inputId = "r_o",
                label = "Seleccione el parametro R0 del modelo:",
                min = 0.5,
                max = 5,
                value = 2.28,
                step = .001
              ),
              sliderInput(
                inputId = "dur",
                label = "Seleccione el rango de días de duración del periodo infeccioso:",
                min = 0,
                max = 14,
                value = c(5, 7),
                step = .1
              ),
              sliderInput(
                inputId = "let",
                label = "Seleccione el porcentaje de letalidad",
                min = 0.005,
                max = .5,
                value = 0.011418095,
                step = .001
              ),
              numericInput(
                inputId ="ndays",
                label = "Seleccione el numero de dias a evaluar:",
                min = 10,
                max = 500,
                value = 300,
                step = 10
              ),
              selectInput(
                "features",
                "Seleccione las variables a graficar:",
                c("Susceptibles", "Infectados", "Recuperados", "Muertes"),
                selected = "Infectados",
                multiple = TRUE
              ),
              helpText("Nota: Aplicación de Shiny basada en el modelo SIR publicado",
                       "por el Instituto Nacional de Salud. Para mayor información,",
                       "consultar el enlace: https://www.ins.gov.co/Noticias/Paginas/Coronavirus.aspx")
            ),
            mainPanel(
              plotlyOutput("report_cases")
            )
          )
        )
      ),
      
      # UI Escenarios -------------------------------------------------------------
      
      tabItem(
      tabName = "comparison",
      fluidPage(
        titlePanel("Escenarios"),
        sidebarLayout(
          sidebarPanel(
            selectInput(
              "select_scenario",
              "Seleccione un escenario:",
              c("Escenario 1", "Escenario 2")
            ),
            helpText(
              "Nota: El escenario 1 asume una intervención de cuarentena durante",
              "33 días, con una efectividad de 47%. El escenario 2 asume una ",
              "intervención de 47 días con una efectividad de 56%.",
              "Aplicación de Shiny basada en el modelo SIR publicado",
              "por el Instituto Nacional de Salud. Para mayor información,",
              "consultar el enlace: https://www.ins.gov.co/Noticias/Paginas/Coronavirus.aspx"
            )
          ),
          mainPanel(
            # plotlyOutput("comp_S"),
            # plotlyOutput("comp_I"),
            plotlyOutput("comp_casos"),
            plotlyOutput("comp_muertes")
            )
          )
        )
      )
    )
  )
)

###########################################################################
# Server ------------------------------------------------------------------
###########################################################################
server <- function(input, output){

# Modelo base -------------------------------------------------------------

  get_data <- reactive({
    CFRm <- input$let
    CFRalfa <<- CFRm * ((CFRm * (1 - CFRm) / (CFRsd ^ 2)) - 1)
    CFRbeta <<- CFRalfa * ((1 - CFRm) / CFRm)
    if(input$city == "Colombia"){
      output <- Covid19Model(mR0 = input$r_o, nDays = input$ndays, nSim = 500, minTransm = input$dur[1], maxTransm = input$dur[2])
    } else if(input$city == "Bogotá"){
      output <- Covid19Model(pob = 7743955, mR0 = input$r_o, nDays = input$ndays, nSim = 500, minTransm = input$dur[1], maxTransm = input$dur[2])
    }
    output
  })
  
  output$report_cases <- renderPlotly({
    # browser()
    output <- get_data()
    output <- lapply(names(output), function(x){
      names(output[[x]]) <- c(
        'day',
        x,
        paste0(x, '_lb'),
        paste0(x, '_ub')
      )
      output[[x]]$day <- NULL
      output[[x]]
    })
    df <- do.call(cbind, output)
    df$day <- 1:nrow(df)
    
    if(is.null(df)) return(NULL)
    if(is.null(input$features)) return(NULL)
    
    p <- plot_ly(data = df, x = ~day)
    for(feat in input$features){
      feat2 <- tolower(feat)
      p <- p %>%
        add_trace(y = df[, feat2], mode = 'lines', name = feat)
    }
    p <- p %>%
      layout(hovermode = 'compare', title = input$city)
  })
  
  # Escenarios --------------------------------------------------------------
  
  get_scenarios <- reactive({
    if(input$select_scenario == 'Escenario 1'){
      scenario <- Covid19ModelAccIntermitente(
        name = "Colombia", pob = 38837139, nInitPtes = 5,Intervention = "Scenario 1",
        Efect = 0.47, nDaysNOInt = 18, pInt = 33,  pNOInt = 249, nDaysTotal = 300,
        nSim = 500)
    } else {
      scenario <- Covid19ModelAccIntermitente(
        name = "Colombia", pob = 38837139, nInitPtes = 5, Intervention = "Scenario 2",
        Efect = 0.56, nDaysNOInt = 25, pInt = 47,  pNOInt = 228, nDaysTotal = 300,
        nSim = 500)

    }
    
    names(scenario) <- c("casos", "muertes")
    scenario <- lapply(names(scenario), function(x){
      names(scenario[[x]]) <- c(
        'day',
        x,
        paste0(x, '_lb'),
        paste0(x, '_ub')
      )
      scenario[[x]]$day <- NULL
      scenario[[x]]
    })
    scenario <- do.call(cbind, scenario)
    scenario$day <- 1:nrow(scenario)

    base_df <- Covid19Model(mR0 = 2.28, nDays = 300, nSim = 500,
                            minTransm = 5, maxTransm = 7)
    
    base_df <- lapply(names(base_df), function(x){
      names(base_df[[x]]) <- c(
        'day',
        x,
        paste0(x, '_lb'),
        paste0(x, '_ub')
      )
      base_df[[x]]$day <- NULL
      base_df[[x]]
    })
    base_df <- do.call(cbind, base_df)
    base_df$casos <- base_df$infectados + base_df$recuperados + base_df$muertes
    base_df <- base_df[, c('casos','muertes')]
    names(base_df) <- paste0('base_', names(base_df))
    
    final_df <- cbind(scenario, base_df)
    
    final_df
  })
  
  output$comp_casos <- renderPlotly({
    # browser()
    final_df <- get_scenarios()

    plot_ly(final_df, x=~day) %>%
      add_trace(y=~casos, mode = 'lines', name = input$select_scenario) %>%
      add_trace(y = ~base_casos, mode = 'lines', name = 'Base') %>%
      layout(title='Total de Casos', hovermode = 'compare')
  })
  
  output$comp_muertes <- renderPlotly({
    # browser()
    final_df <- get_scenarios()
    
    plot_ly(final_df, x=~day) %>%
      add_trace(y=~muertes, mode = 'lines', name = input$select_scenario) %>%
      add_trace(y = ~base_muertes, mode = 'lines', name = 'Base') %>%
      layout(title='Total de Muertes', hovermode = 'compare')
  })
}

shinyApp(ui, server)
