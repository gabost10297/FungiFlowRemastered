# FungiFlow app: upload barcodes, run the pipeline, view the pictures.
#
# Start it from the project root with:
#   bash start.sh
# then open http://127.0.0.1:3838 (stop with: docker compose -f docker/compose.yml down).
# Uploads go to runs/<name>/raw_data/; run_pipeline.sh fills runs/<name>/results/.

library(shiny)

source("tree.R", local = FALSE)
source("plots.R", local = FALSE)
source("runs.R", local = FALSE)

options(shiny.maxRequestSize = 20 * 1024^3)

PROJECT <- project_root()

run_tab <- tabPanel(
  "Run",
  sidebarLayout(
    sidebarPanel(
      h4("New run"),
      textInput("new_run", "Name (letters, digits, _ and -)", placeholder = "air_2025"),
      actionButton("create", "Create"),
      hr(),
      selectInput("run", "Run", choices = NULL),
      fileInput("fastq", "Add barcodes (.fastq.gz)", multiple = TRUE, accept = ".gz"),
      tableOutput("files"),
      actionButton("start", "Start pipeline", class = "btn-primary"),
      actionButton("stop", "Stop", class = "btn-danger"),
      p(style = "margin-top: 1em; color: #555;",
        "Step 06 takes hours. You can close the browser; the run continues.",
        "Starting again continues where it stopped.")
    ),
    mainPanel(
      h4(textOutput("run_title")),
      tableOutput("status"),
      h4("Log"),
      verbatimTextOutput("log")
    )
  )
)

pictures_tab <- tabPanel(
  "Pictures",
  sidebarLayout(
    sidebarPanel(
      selectInput("picture_run", "Run", choices = NULL),
      selectInput("view", "Picture",
                  choices = c("Tree" = "tree",
                              "Circular tree with reads" = "circular",
                              "Organism heatmap" = "heatmap",
                              "PCA of barcodes" = "pca",
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
      conditionalPanel(
        "input.view == 'circular'",
        sliderInput("circular_label_size", "Label size",
                    min = 1.5, max = 4, value = 2.6, step = 0.1)
      ),
      conditionalPanel(
        "input.view == 'pca'",
        sliderInput("loadings", "Organism arrows",
                    min = 1, max = 15, value = 6, step = 1)
      ),
      downloadButton("pdf", "Download PDF"),
      p(style = "margin-top: 1.5em; color: #555;",
        textOutput("note"))
    ),
    mainPanel(plotOutput("picture", height = "900px"))
  )
)

ui <- navbarPage("FungiFlow", run_tab, pictures_tab)

server <- function(input, output, session) {
  refresh <- reactiveVal(0)
  tick <- reactiveTimer(3000)
  created <- reactiveVal(NULL)

  observe({
    refresh()
    runs <- list_runs(PROJECT)
    selected <- isolate(if (is.null(created())) input$run else created())
    created(NULL)
    updateSelectInput(session, "run", choices = runs,
                      selected = if (isTRUE(selected %in% runs)) selected else runs[1])
  })

  observe({
    refresh()
    tick()
    runs <- runs_with_results(PROJECT)
    selected <- isolate(input$picture_run)
    if (identical(runs, isolate(session$userData$picture_runs))) return()
    session$userData$picture_runs <- runs
    updateSelectInput(session, "picture_run", choices = runs,
                      selected = if (isTRUE(selected %in% runs)) selected else runs[1])
  })

  observeEvent(input$create, {
    name <- trimws(input$new_run)
    if (!valid_run_name(name)) {
      showNotification("Use only letters, digits, _ and -.", type = "error")
      return()
    }
    dir.create(file.path(PROJECT, "runs", name, "raw_data"), recursive = TRUE, showWarnings = FALSE)
    created(name)
    refresh(refresh() + 1)
    updateTextInput(session, "new_run", value = "")
  })

  observeEvent(input$fastq, {
    req(input$run)
    if (is_running(input$run)) {
      showNotification("This run is in progress. Stop it before adding barcodes.", type = "error")
      return()
    }
    result <- add_files(PROJECT, input$run, input$fastq)
    if (length(result$added)) {
      showNotification(paste("Added", length(result$added), "file(s)."), type = "message")
    }
    if (length(result$skipped)) {
      showNotification(paste("Skipped (not .fastq.gz):", paste(result$skipped, collapse = ", ")),
                       type = "warning", duration = 10)
    }
    refresh(refresh() + 1)
  })

  observeEvent(input$start, {
    req(input$run)
    if (!nrow(run_files(PROJECT, input$run))) {
      showNotification("Add barcodes first.", type = "error")
      return()
    }
    if (is_running(input$run)) {
      showNotification("Already running.", type = "warning")
      return()
    }
    tryCatch({
      start_run(PROJECT, input$run)
      showNotification("Pipeline started.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = NULL))
    refresh(refresh() + 1)
  })

  observeEvent(input$stop, {
    req(input$run)
    stop_run(PROJECT, input$run)
    showNotification("Stopped.", type = "warning")
    refresh(refresh() + 1)
  })

  output$files <- renderTable({
    refresh()
    req(input$run)
    run_files(PROJECT, input$run)
  })

  running <- reactive({
    tick()
    refresh()
    req(input$run)
    is_running(input$run)
  })

  output$run_title <- renderText({
    req(input$run)
    paste0(input$run, if (running()) ": running" else "")
  })

  output$status <- renderTable({
    tick()
    req(input$run)
    status <- run_status(PROJECT, input$run)
    if (!running()) status$state[status$state == "running"] <- "interrupted"
    status
  })

  output$log <- renderText({
    tick()
    req(input$run)
    run_log(PROJECT, input$run)
  })

  picture_dir <- reactive({
    req(input$picture_run)
    run_dir(PROJECT, input$picture_run)
  })

  picture <- reactive({
    dir <- picture_dir()
    tryCatch(
      switch(input$view,
        tree = draw_tree(dir,
                         bootstrap = input$bootstrap,
                         color_by = input$color_by,
                         label_size = input$label_size),
        circular = draw_circular(dir, label_size = input$circular_label_size),
        heatmap = draw_heatmap(dir, top_n = input$top_n),
        pca = draw_pca(dir, loadings = input$loadings),
        persample = draw_per_sample(dir),
        its = draw_same_its(dir)
      ),
      error = function(e) validate(need(FALSE, conditionMessage(e)))
    )
  })

  output$picture <- renderPlot(picture(), res = 110, height = function() {
    if (identical(input$view, "heatmap")) return(280 + 22 * input$top_n)
    if (input$view %in% c("tree", "circular")) return(900)
    if (identical(input$view, "pca")) return(750)
    620
  })

  output$note <- renderText({
    switch(input$view,
      tree = "Top 30 organisms from step 12. Grey dots are nodes below the bootstrap cutoff.",
      circular = "The step 12 tree with one ring per barcode. Grey cells are absences.",
      heatmap = "One row per organism, summed across its OTUs. Grey cells are absences.",
      pca = "Barcodes close together have similar organism shares. Point size is the organism count.",
      persample = "The number on each bar is the organism count. Orange is found only in that barcode.",
      its = "Share of reads whose ITS matches the representative at 98.5% or more."
    )
  })

  output$pdf <- downloadHandler(
    filename = function() {
      run <- gsub("[^A-Za-z0-9_-]", "", gsub(" ", "_", input$picture_run))
      paste0("fungiflow_", run, "_", input$view, ".pdf")
    },
    content = function(file) {
      ggsave(file, picture(), width = 12, height = 9, device = "pdf")
    }
  )
}

shinyApp(ui, server)
