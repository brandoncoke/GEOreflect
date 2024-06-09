library(tidyverse)
library(shiny)
library(leaflet)

df <- data.frame(
  Number_Total = sample(c(5, 6, 1, 3)),
  Species = sample(c("Ilione trifaria", "Pherbellia argyrotarsis", 
                     "Euthycera seguyi", "Ilione trifaria")),
  X= sample(c(37, 28, 21, 30)),
  Y= sample(c(-5, -16, -10, -15))
)

ui <- (fluidPage(titlePanel("Species Checker"),  
                 sidebarLayout(
                   sidebarPanel(
                     selectizeInput('species', 'Choose species', 
                                    choices = df$Species, multiple = TRUE, 
                                    options = list(placeholder = 'select species'))
                   ),
                   mainPanel(
                     leafletOutput("CountryMap", width = 600, height = 600))
                 )
))

server <- function(input, output, session) {
  
  map_data <- reactive({
    #req(input$species)
    df[df$Species %in% input$species, ]
  })
  
  output$CountryMap <- renderLeaflet({
    
    leaflet() %>% addTiles() %>% 
      setView(lng = 20, lat = 40, zoom = 2) %>%
      addCircles(lng = map_data() %>% pull(Y), lat = map_data() %>% pull(X), weight = 10, 
                 radius = sqrt(map_data() %>% pull(Number_Total))*15000, 
                 popup = map_data() %>% pull(Species))
    
  })
  
}
# Create Shiny app ----
shinyApp(ui, server)

