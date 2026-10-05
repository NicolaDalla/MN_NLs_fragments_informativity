library(shiny)
library(igraph)
library(tidyverse)
library(Spectra)
library(MsBackendMgf)
library(qgam)
library(MsCoreUtils)
library(ggplot2)
library(plotly)
library(visNetwork)
library(viridis)

source("MN_info.R")

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  titlePanel("Interactive MS/MS Molecular Network & Feature Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Data Input"),
      fileInput("msms_file", "Spectra Object (.rds or .mgf)", accept = c(".rds", ".mgf")),
      fileInput("matrix_file", "Similarity Matrix (.rds or .csv)", accept = c(".rds", ".csv")),
      
      hr(),
      h4("2. Network Parameters"),
      numericInput("thr", "Similarity Cutoff (thr):", value = 0.7, min = 0, max = 1, step = 0.05),
      numericInput("tolerance", "m/z Tolerance:", value = 0.005, step = 0.001),
      numericInput("min_presence", "Min Feature Prevalence:", value = 10, min = 1),
      numericInput("perm", "Permutations (N):", value = 50, min = 10, max = 500),
      
      hr(),
      h4("3. QGAM Settings (Optional)"),
      numericInput("Q", "Quantile (Q):", value = 0.5, min = 0.1, max = 0.9, step = 0.1),
      numericInput("K", "Basis Dimension (K):", value = 5, min = 3, max = 20),
      
      br(),
      actionButton("run_btn", "Run Analysis", class = "btn-primary btn-block")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "Interactive Feature Plot",
          br(),
          helpText("Click on any point to view that feature in the 'Selected Feature Network' tab."),
          plotlyOutput("qgam_plot", height = "550px")
        ),
        
        tabPanel(
          "Selected Feature Network",
          br(),
          fluidRow(
            column(8, selectInput("feature_select", "Selected Feature:",
                                  choices = c("Run analysis first" = ""), width = "100%")),
            column(4, br(), uiOutput("feature_info_badge"))
          ),
          hr(),
          visNetworkOutput("network_plot", height = "600px")
        ),
        
        tabPanel(
          "Custom Search Network",
          br(),
          wellPanel(
            h5(strong("Enter custom values to highlight nodes")),
            fluidRow(
              column(4, textInput("custom_fi", "Fragment ions m/z (e.g. 123.04, 147.08):", value = "")),
              column(4, textInput("custom_nl", "Neutral losses (e.g. 176.03, 162.05):", value = "")),
              column(4, numericInput("custom_tol", "Match tolerance (m/z):", value = 0.005, step = 0.001))
            ),
            helpText("Fragment ion = red | Neutral loss = blue | Both = purple | No match = grey")
          ),
          hr(),
          visNetworkOutput("custom_network_plot", height = "600px")
        ),
        
        tabPanel(
          "Merged Results Table",
          br(),
          downloadButton("download_csv", "Export Table (.csv)", class = "btn-success"),
          br(), br(),
          tableOutput("results_table")
        )
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================
server <- function(input, output, session) {
  
  processed_spectra <- reactiveVal(NULL)
  processed_matrix  <- reactiveVal(NULL)
  alldf_data        <- reactiveVal(NULL)
  
  # ----------------------------------------------------------------------------
  # RUN ANALYSIS (always executes on button press, independent of the open tab)
  # ----------------------------------------------------------------------------
  observeEvent(input$run_btn, {
    if (is.null(input$msms_file) || is.null(input$matrix_file)) {
      showNotification("Please upload both the spectra file and the similarity matrix.",
                       type = "warning")
      return()
    }
    
    tryCatch({
      withProgress(message = "Processing network features...", value = 0, {
        
        incProgress(0.1, detail = "Loading spectra...")
        ext_spec <- tolower(tools::file_ext(input$msms_file$name))
        msms_spectra <- if (ext_spec == "rds") {
          readRDS(input$msms_file$datapath)
        } else {
          Spectra(input$msms_file$datapath, source = MsBackendMgf::MsBackendMgf())
        }
        
        id_var <- if ("feature_id" %in% spectraVariables(msms_spectra)) msms_spectra$feature_id else msms_spectra$id
        if (is.null(id_var)) id_var <- seq_along(msms_spectra)
        msms_spectra$id         <- as.character(id_var)
        msms_spectra$feature_id <- msms_spectra$id
        
        # Neutral losses (kept because MN_info may read msms$nl)
        neutral_loss <- function(x, precursorMz, ...) {
          x[, "mz"] <- precursorMz - x[, "mz"]
          x[order(x[, "mz"]), , drop = FALSE]
        }
        nl_proc <- addProcessing(msms_spectra, neutral_loss, spectraVariables = "precursorMz")
        msms_spectra$nl <- mz(nl_proc)
        
        incProgress(0.2, detail = "Loading matrix...")
        ext_mat <- tolower(tools::file_ext(input$matrix_file$name))
        sim_matrix <- if (ext_mat == "rds") {
          readRDS(input$matrix_file$datapath)
        } else {
          as.matrix(read.csv(input$matrix_file$datapath, header = TRUE, row.names = 1))
        }
        
        if (nrow(sim_matrix) != length(msms_spectra) || ncol(sim_matrix) != length(msms_spectra)) {
          stop("Matrix size (", nrow(sim_matrix), " x ", ncol(sim_matrix),
               ") does not match the number of spectra (", length(msms_spectra), ").")
        }
        colnames(sim_matrix) <- rownames(sim_matrix) <- msms_spectra$id
        
        processed_spectra(msms_spectra)
        processed_matrix(sim_matrix)
        
        incProgress(0.4, detail = "Calculating fragment ions...")
        exfi <- MN_info(method = "FIs", matrix = sim_matrix, msms = msms_spectra,
                        perm = input$perm, thr = input$thr, min_presence = input$min_presence,
                        Q = input$Q, K = input$K, tolerance = input$tolerance)
        
        incProgress(0.7, detail = "Calculating neutral losses...")
        exnl <- MN_info(method = "NLs", matrix = sim_matrix, msms = msms_spectra,
                        perm = input$perm, thr = input$thr, min_presence = input$min_presence,
                        Q = input$Q, K = input$K, tolerance = input$tolerance)
        
        nls <- exnl %>% dplyr::select(value, number, z_score, excess) %>%
          mutate(origin = "NL", log_number = log10(number))
        mzs <- exfi %>% dplyr::select(value, number, z_score, excess) %>%
          mutate(origin = "mz", log_number = log10(number))
        
        alldf <- rbind(nls, mzs) %>%
          arrange(desc(z_score)) %>%
          mutate(key = paste0(origin, "_", value))
        
        feature_choices <- setNames(
          alldf$key,
          paste0(alldf$origin, ": ", round(alldf$value, 4), " (Z-score: ", round(alldf$z_score, 2), ")")
        )
        updateSelectInput(session, "feature_select",
                          choices = feature_choices, selected = feature_choices[1])
        
        alldf_data(alldf)
        incProgress(1, detail = "Done!")
      })
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)), type = "error", duration = NULL)
    })
  })
  
  # ----------------------------------------------------------------------------
  # SCATTER PLOT + CLICK -> Tab 2
  # ----------------------------------------------------------------------------
  output$qgam_plot <- renderPlotly({
    df <- alldf_data()
    validate(need(!is.null(df), "Upload the files and press 'Run Analysis'."))
    
    p <- ggplot(df, aes(
      x = log_number, y = z_score, color = z_score, key = key,
      text = paste0(
        "<b>Origin:</b> ", origin, "<br>",
        "<b>Value (m/z or NL):</b> ", round(value, 4), "<br>",
        "<b>Node count:</b> ", number, "<br>",
        "<b>Z-modularity:</b> ", round(z_score, 3)
      )
    )) +
      geom_point(aes(shape = origin), size = 3.5, alpha = 0.85) +
      scale_color_viridis_c(
        option = "plasma", 
        direction = 1,        # High values map to bright yellow
        name = "Z-score"
      )+
      scale_shape_manual(values = c("mz" = 16, "NL" = 17)) +
      labs(title = "Diagnostic Z-Modularity Distribution",
           x = "Log10(Occurrence Frequency)", y = "Z-Modularity Score") +
      theme_minimal(base_size = 13)
    
    ggplotly(p, tooltip = "text", source = "qgam_click") %>%
      event_register("plotly_click")
  })
  
  click_ev <- reactive({
    suppressWarnings(event_data("plotly_click", source = "qgam_click"))
  })
  
  observeEvent(click_ev(), {
    ev <- click_ev()
    df <- alldf_data()
    req(ev, df)
    
    sel_key <- NULL
    if (!is.null(ev$key) && length(ev$key) > 0 && nzchar(ev$key[[1]])) {
      sel_key <- as.character(ev$key[[1]])
    } else {
      m <- df[abs(df$log_number - ev$x) < 1e-4 & abs(df$z_score - ev$y) < 1e-4, ]
      if (nrow(m) > 0) sel_key <- m$key[1]
    }
    
    if (!is.null(sel_key) && sel_key %in% df$key) {
      updateSelectInput(session, "feature_select", selected = sel_key)
      updateTabsetPanel(session, "main_tabs", selected = "Selected Feature Network")
    }
  })
  
  # ----------------------------------------------------------------------------
  # SHARED GRAPH DATA (graph, layout, fragment and NL lists)
  # ----------------------------------------------------------------------------
  graph_data <- reactive({
    msms <- processed_spectra()
    mat  <- processed_matrix()
    validate(need(!is.null(msms) && !is.null(mat), "Press 'Run Analysis' first."))
    
    mat[is.na(mat)]      <- 0
    mat[mat < input$thr] <- 0
    diag(mat)            <- 0
    
    g <- graph_from_adjacency_matrix(mat, mode = "max", weighted = TRUE, diag = FALSE)
    g <- delete_vertices(g, V(g)[degree(g) == 0])
    validate(need(vcount(g) > 0, "No edges above the similarity cutoff."))
    
    ids <- V(g)$name
    idx <- match(ids, msms$id)
    
    frag_all <- as.list(mz(msms))
    pmz      <- precursorMz(msms)
    nl_all   <- Map(function(f, p) p - f, frag_all, pmz)
    
    set.seed(123)
    lay <- layout_with_fr(g)
    
    nodes <- data.frame(
      id          = ids,
      precursorMz = pmz[idx],
      rt          = if ("rtime" %in% spectraVariables(msms)) rtime(msms)[idx] / 60 else NA_real_,
      x           = lay[, 1] * 100,
      y           = lay[, 2] * 100,
      stringsAsFactors = FALSE
    )
    
    list(nodes = nodes,
         edges = igraph::as_data_frame(g, what = "edges"),
         frag  = frag_all[idx],
         nl    = nl_all[idx])
  })
  
  has_match <- function(val_list, targets, tol) {
    if (length(targets) == 0) return(rep(FALSE, length(val_list)))
    vapply(val_list, function(v) {
      v <- as.numeric(v)
      if (length(v) == 0) return(FALSE)
      any(vapply(targets, function(t) any(abs(v - t) <= tol, na.rm = TRUE), logical(1)))
    }, logical(1))
  }
  
  draw_network <- function(nodes, edges, legend = NULL) {
    vn <- visNetwork(nodes, edges) %>%
      visNodes(shape = "dot", borderWidth = 1) %>%
      visEdges(color = list(color = "#d3d3d3", highlight = "#2b7ce9"), width = 1.5) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = TRUE) %>%
      visInteraction(dragNodes = TRUE, zoomView = TRUE) %>%
      visPhysics(enabled = FALSE)
    if (!is.null(legend)) vn <- vn %>% visLegend(addNodes = legend, useGroups = FALSE)
    vn
  }
  
  # ----------------------------------------------------------------------------
  # TAB 2
  # ----------------------------------------------------------------------------
  output$network_plot <- renderVisNetwork({
    validate(need(nzchar(input$feature_select), "Press 'Run Analysis' and select a feature."))
    gd <- graph_data()
    
    sel_origin <- sub("_.*$", "", input$feature_select)
    sel_value  <- as.numeric(sub("^(mz|NL)_", "", input$feature_select))
    
    vals <- if (sel_origin == "NL") gd$nl else gd$frag
    hit  <- has_match(vals, sel_value, input$tolerance)
    
    nodes <- gd$nodes
    nodes$color.background <- ifelse(hit, "firebrick", "lightgrey")
    nodes$color.border     <- ifelse(hit, "darkred", "grey50")
    nodes$size             <- ifelse(hit, 22, 10)
    nodes$title <- paste0(
      "<b>Node ID:</b> ", nodes$id, "<br>",
      "<b>Precursor m/z:</b> ", round(nodes$precursorMz, 4), "<br>",
      "<b>Retention time:</b> ", round(nodes$rt, 2), " min<br>",
      "<b>Contains selected feature:</b> ", ifelse(hit, "YES", "No")
    )
    
    legend <- data.frame(label = c("Contains feature", "Other"),
                         color = c("firebrick", "lightgrey"), shape = "dot")
    draw_network(nodes, gd$edges, legend)
  })
  
  # ----------------------------------------------------------------------------
  # TAB 3
  # ----------------------------------------------------------------------------
  output$custom_network_plot <- renderVisNetwork({
    gd  <- graph_data()
    tol <- if (is.null(input$custom_tol) || is.na(input$custom_tol)) 0.005 else input$custom_tol
    
    parse_num_inputs <- function(txt) {
      if (is.null(txt) || !nzchar(trimws(txt))) return(numeric(0))
      vals <- suppressWarnings(as.numeric(unlist(strsplit(trimws(txt), "[,;[:space:]]+"))))
      vals[!is.na(vals)]
    }
    
    has_fi <- has_match(gd$frag, parse_num_inputs(input$custom_fi), tol)
    has_nl <- has_match(gd$nl,   parse_num_inputs(input$custom_nl), tol)
    
    nodes <- gd$nodes
    nodes$color.background <- dplyr::case_when(
      has_fi & has_nl ~ "#8A2BE2",
      has_fi          ~ "#CD5C5C",
      has_nl          ~ "#1E4FD8",
      TRUE            ~ "#D3D3D3"
    )
    nodes$color.border <- ifelse(has_fi | has_nl, "#333333", "#A9A9A9")
    nodes$size         <- ifelse(has_fi | has_nl, 22, 10)
    nodes$title <- paste0(
      "<b>Node ID:</b> ", nodes$id, "<br>",
      "<b>Precursor m/z:</b> ", round(nodes$precursorMz, 4), "<br>",
      "<b>Retention time:</b> ", round(nodes$rt, 2), " min<br>",
      "<b>Matched FI:</b> ", ifelse(has_fi, "YES", "No"), "<br>",
      "<b>Matched NL:</b> ", ifelse(has_nl, "YES", "No")
    )
    
    legend <- data.frame(label = c("Fragment ion", "Neutral loss", "Both", "No match"),
                         color = c("#CD5C5C", "#1E4FD8", "#8A2BE2", "#D3D3D3"), shape = "dot")
    draw_network(nodes, gd$edges, legend)
  })
  
  # ----------------------------------------------------------------------------
  # Badge, table, download
  # ----------------------------------------------------------------------------
  output$feature_info_badge <- renderUI({
    req(nzchar(input$feature_select))
    origin <- sub("_.*$", "", input$feature_select)
    value  <- as.numeric(sub("^(mz|NL)_", "", input$feature_select))
    tags$div(
      style = "padding: 8px; background-color: #f8f9fa; border-radius: 5px; border-left: 4px solid #d9534f;",
      tags$b("Highlighted feature: "), origin, " @ ", round(value, 4), " m/z"
    )
  })
  
  output$results_table <- renderTable({
    req(alldf_data())
    head(alldf_data() %>% dplyr::select(-key) %>% arrange(desc(z_score)), 35)
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("diagnostic_features_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(alldf_data() %>% dplyr::select(-key), file, row.names = FALSE)
  )
}

shinyApp(ui = ui, server = server)
