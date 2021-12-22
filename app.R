#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(plotly)
library(ggcorrplot)
library(shiny)
library(shinythemes)
library(Biostrings)
library(tibble)
library(dplyr)
library(colorRamps)
library(DT)
library(seqRFLP)
library(msa)
library(seqinr)
library(DT)
library(shinythemes)
library(xlsx)
library(openxlsx)

csv_db <- read.csv('data/bioinformatics_database.csv')
csv_db <- head(csv_db, n=20)
granges <- unique(csv_db$company)

ui <- fluidPage(
    theme = shinytheme("cerulean"),
    # Application title
    titlePanel("Anàlisis filogenètic Biofar"),

    sidebarLayout(
        selectInput(
            inputId='granja',
            label='escull granja',
            choices = granges,
            selected = granges,
            multiple=TRUE
        ),

        # Show a plot of the generated distribution
        mainPanel(h1("Matriu de distàncies"),
                  dataTableOutput("distances"),
                  h1('Matriu de similituds'),
                  dataTableOutput('similarity'),
                  br(),
                  downloadButton("downloadData", "Descarregar taules format Excel"),
                  br(),
                  h1("Arbre filogenètic"),
                  plotOutput("phylogenetic_tree"))



    )
)

server <- function(input, output) {

    generate_table_distance <- function(granges){
        csv_db <- read.csv('data/bioinformatics_database.csv')
        csv_db <- head(csv_db, n=20)
        csv_db$seq_len <- unlist(lapply(csv_db$sequence, function(x) nchar(x)))
        csv_db <- csv_db[csv_db$company %in% granges,]
        names <- csv_db$identification
        sequences <- csv_db$sequencenames <- csv_db$identification
        sequences <- csv_db$sequence
        df_db <- data.frame(names,sequences)
        fasta_db <- dataframe2fas(df_db, file='data/df_db.fasta')
        mySeqs <- readAAStringSet('data/df_db.fasta')
        ### Perform a multiple sequence alignment
        myAln <- msa(mySeqs)
        # this uses all the default settings and the CLUSTALW algorithm
        ### Turn your alignment into a tree
        # convert the alignment for the seqinr package
        myAln2 <- msaConvert(myAln, type="seqinr::alignment")
        # this object is a list object with 4 elements

        # generate a distance matrix using seqinr package
        d <- dist.alignment(myAln2, "identity")


        return(d)

    }

     output$similarity <- renderDataTable({
        d <- generate_table_distance(input$granja)
        # have a look at the output
        df_d <- as.data.frame(as.matrix(d))
        df_d <- 1-df_d
        df_d <- round(df_d, 3)


        df_d <- rownames_to_column(df_d, "identification")
        select_csv_db <- csv_db[c('company','explotation','identification','date', 'orf5','seq_len')]
        joined_df_d <- dplyr::left_join(df_d, select_csv_db)
        joined_df_d <- joined_df_d %>% select(seq_len, everything())
        joined_df_d <- joined_df_d %>% select(orf5, everything())
        joined_df_d <- joined_df_d %>% select(date, everything())
        joined_df_d <- joined_df_d %>% select(explotation, everything())
        joined_df_d <- joined_df_d %>% select(company, everything())
        brks <- quantile(df_d[-1], probs = seq(.05, .95, .05), na.rm = TRUE)
        clrs <- green2red(length(brks)+1)
        clrs <- rev(clrs)
        df_color <- datatable(joined_df_d) %>% formatStyle(names(joined_df_d)[c(-1,-2,-3,-4,-5)], backgroundColor = styleInterval(brks, clrs))
        df_color
    })


    output$distances <- renderDataTable({
        d <- generate_table_distance(input$granja)
        # have a look at the output
        df_d <- as.data.frame(as.matrix(d))
        df_d <- round(df_d, 3)

        df_d <- rownames_to_column(df_d, "identification")
        select_csv_db <- csv_db[c('company','explotation','identification','date', 'orf5','seq_len')]
        joined_df_d <- dplyr::left_join(df_d, select_csv_db)
        joined_df_d <- joined_df_d %>% select(seq_len, everything())
        joined_df_d <- joined_df_d %>% select(orf5, everything())
        joined_df_d <- joined_df_d %>% select(date, everything())
        joined_df_d <- joined_df_d %>% select(explotation, everything())
        joined_df_d <- joined_df_d %>% select(company, everything())
        brks <- quantile(df_d[-1], probs = seq(.05, .95, .05), na.rm = TRUE)
        clrs <- green2red(length(brks)+1)
        df_color <- datatable(joined_df_d) %>% formatStyle(names(joined_df_d)[c(-1,-2,-3,-4,-5)], backgroundColor = styleInterval(brks, clrs))
        df_color})


    df_d_s <- reactive({
        d <- generate_table_distance(input$granja)
        # similarity download
        df_s <- as.data.frame(as.matrix(d))
        df_s <- round(df_s, 3)

        df_s <- rownames_to_column(df_s, "identification")
        select_csv_db <- csv_db[c('company','explotation','identification','date', 'orf5','seq_len')]
        joined_df_s <- dplyr::left_join(df_s, select_csv_db)
        joined_df_s <- joined_df_s %>% select(seq_len, everything())
        joined_df_s <- joined_df_s %>% select(orf5, everything())
        joined_df_s <- joined_df_s %>% select(date, everything())
        joined_df_s <- joined_df_s %>% select(explotation, everything())
        joined_df_s <- joined_df_s %>% select(company, everything())

        # distance download
        df_d <- as.data.frame(as.matrix(d))
        df_d <- 1-df_d
        df_d <- round(df_d, 3)

        df_d <- rownames_to_column(df_d, "identification")
        select_csv_db <- csv_db[1:5]
        joined_df_d <- dplyr::left_join(df_d, select_csv_db)
        joined_df_d <- joined_df_d %>% select(orf5, everything())
        joined_df_d <- joined_df_d %>% select(date, everything())
        joined_df_d <- joined_df_d %>% select(explotation, everything())
        joined_df_d <- joined_df_d %>% select(company, everything())

        df_d_s <- list(joined_df_d, joined_df_s)

        return(df_d_s)})

    output$phylogenetic_tree <- renderPlot({
        d <- generate_table_distance(input$granja)
        h_cluster <- hclust(d, method = "average", members = NULL)
        plot(h_cluster, ann=FALSE)
    })

    output$downloadData <- downloadHandler(
        filename = function(){paste0(gsub(' ', '_', Sys.time()), 'shiny_phylogenetic_analysis.xlsx')},
        content = function(file) {
            wb <- createWorkbook()
            addWorksheet(wb, "distàncies")
            df<- df_d_s()[[2]]


            conditionalFormatting(wb, "distàncies",
                                  cols = 1:ncol(df), rows = 1:nrow(df)+1,
                                  style = green2red(3),
                                  type = "colourScale"
            )

            addStyle(wb, sheet = 'distàncies',
                     style = createStyle(valign='center', halign = "center", border='TopBottomLeftRight', borderStyle='thin'),
                     rows = 0:nrow(df)+1, cols = 1:ncol(df), gridExpand = TRUE)

            addStyle(wb, sheet = 'distàncies',
                     style = createStyle(valign='center', halign = "center", border='TopBottomLeftRight', borderStyle='thin', textDecoration='bold'),
                     rows = 0:1, cols = 1:ncol(df), gridExpand = TRUE)

            addStyle(wb, sheet = 'distàncies',
                     style = createStyle(textRotation=45, border='TopBottomLeftRight', borderStyle='thin'),
                     rows = 0:1, cols = 6:ncol(df), gridExpand = TRUE)

            setColWidths(wb, sheet='distàncies', cols = 1:ncol(df), widths = "auto")
            setRowHeights(wb, 'distàncies', rows = 1, heights = 40)
            writeData(wb, "distàncies", df, colNames = TRUE) ## write data.frame

            addWorksheet(wb, "similituds")
            df<- df_d_s()[[1]]


            conditionalFormatting(wb, "similituds",
                                  cols = 1:ncol(df), rows = 1:nrow(df)+1,
                                  style = rev(green2red(3)),
                                  type = "colourScale"
            )

            addStyle(wb, sheet = 'similituds',
                     style = createStyle(valign='center', halign = "center", border='TopBottomLeftRight', borderStyle='thin'),
                     rows = 0:nrow(df)+1, cols = 1:ncol(df), gridExpand = TRUE)

            addStyle(wb, sheet = 'similituds',
                     style = createStyle(valign='center', halign = "center", border='TopBottomLeftRight', borderStyle='thin', textDecoration='bold'),
                     rows = 0:1, cols = 1:ncol(df), gridExpand = TRUE)

            addStyle(wb, sheet = 'similituds',
                     style = createStyle(textRotation=45, border='TopBottomLeftRight', borderStyle='thin'),
                     rows = 0:1, cols = 6:ncol(df), gridExpand = TRUE)

            setColWidths(wb, sheet='similituds', cols = 1:ncol(df), widths = "auto")
            setRowHeights(wb, 'similituds', rows = 1, heights = 40)
            writeData(wb, "similituds", df, colNames = TRUE) ## write data.frame

            saveWorkbook(wb, file=file, TRUE)
        },
    )

}



# Run the application
shinyApp(ui = ui, server = server)
