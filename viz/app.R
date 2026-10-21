# One page for the tree and the step 10-11 summaries.
#
# From the project root, in the fungiflow-viz image:
#   docker build -t fungiflow-viz viz
#   docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -p 3838:3838 -v "$PWD":/data fungiflow-viz
# Then open http://127.0.0.1:3838

library(shiny)

source("tree.R", local = FALSE)
source("plots.R", local = FALSE)

ui <- fluidPage(
  titlePanel("FungiFlow"),
  sidebarLayout(
    sidebarPanel(
      selectInput("view", "Picture",
                  choices = c("Tree" = "tree",
                              "Organism heatmap" = "heatmap",
                              "Organisms per barcode" = "persample",
                              "Same ITS in the top 10" = "its")),
      conditionalPanel(
        "input.view == 'tree'",
        sliderInput("bootstrap", "Show bootstrap at or above",
                    min = 0, max = 100, value = 70, step = 5),
        selectInput("color_by", "Color tips by",
                    choices = c("Phylum" = "phylum", "Class" = "class", "Reads" = "reads")),
        sliderInput("label_size", "Label size",
                    min = 2, max = 6, value = 3.2, step = 0.2)
      ),
      conditionalPanel(
        "input.view == 'heatmap'",
        sliderInput("top_n", "Organisms shown",
                    min = 10, max = 40, value = 20, step = 5)
      ),
      downloadButton("pdf", "Download PDF"),
      p(style = "margin-top: 1.5em; color: #555;",
        textOutput("note"))
    ),
    mainPanel(plotOutput("picture", height = "900px"))
  )
)

server <- function(input, output) {
  picture <- reactive({
    switch(input$view,
      tree = draw_tree(project_root(),
                       bootstrap = input$bootstrap,
                       color_by = input$color_by,
                       label_size = input$label_size),
      heatmap = draw_heatmap(project_root(), top_n = input$top_n),
      persample = draw_per_sample(project_root()),
      its = draw_same_its(project_root())
    )
  })

  output$picture <- renderPlot(picture(), res = 110, height = function() {
    if (identical(input$view, "heatmap")) return(280 + 22 * input$top_n)
    if (identical(input$view, "tree")) return(900)
    620
  })

  output$note <- renderText({
    switch(input$view,
      tree = "Top 30 organisms from step 12. Grey dots are nodes below the bootstrap cutoff.",
      heatmap = "One row per organism, summed across its OTUs. Grey cells are absences.",
      persample = "The number on each bar is the organism count. Orange is found only in that barcode.",
      its = "Share of reads whose ITS matches the representative at 98.5% or more."
    )
  })

  output$pdf <- downloadHandler(
    filename = function() paste0("fungiflow_", input$view, ".pdf"),
    content = function(file) {
      dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
      ggsave(file, picture(), width = 12, height = 9, device = "pdf")
    }
  )
}

shinyApp(ui, server)
