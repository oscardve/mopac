library(shiny)
library(visNetwork)
library(shinyBS)
library(igraph)

# library(shinythemes)
# library(shinyFiles)
# library(ggplot2)
# library(DT)
# library(STRINGdb)
# library(sqldf)
# library(rhandsontable)

shinyUI(navbarPage(title="MoPAC",id="MoPAC",theme=shinythemes::shinytheme("cerulean"),

   ### INPUT ###

   tabPanel("FASTQ alignment",
        column(4,
               h3("Objectives"),
               tags$ul(
                 tags$li("To map the reads in all FASTQ files contained in a specified directory to a custom sgRNA library."),
                 tags$li("To visualize the sequencing depth and percentage of reads mapped in each file."),
                 tags$li("To annotate the conditions and replicates of the experiment.")
               ),
               h3("Instructions"),
               tags$ol(
                 tags$li("Select a working directory in which all output is to be stored."),
                 tags$li("Select the directory containing the FASTQ files to be mapped."),
                 tags$li("Load an sgRNA library file containing the following two columns: 'sgRNA' and 'Gene'."),
                 "If the library includes gene categories, please add a third column called 'Category' containing the following identifiers: 'Essential', 'Nonessential', 'Nontargeting' or 'Other'.",
                 tags$li("Select the starting and ending position of the spacer in the FASTQ file sequences, and specify whether the spacer is to be reversed and/or complemented."),
                 tags$li("Click on RUN."),
                 tags$li("Once finished, go to the 'File annotation' tab and follow the instructions shown.")# fill out the condition and replicate information in the table provided. Finally, click on the red button above the table to save it.")
               ),
          wellPanel(
          tags$style(".popover{max-width: 100%;}"),
          fluidRow(column(8,uiOutput("dir_out_text1")),
                   column(3,shinyFiles::shinyDirButton("dir_out",label="Browse...",title="Please select a working directory:",class="btn-secondary"))),
          textOutput("dir_out_text2"),HTML("<br>"),
          fluidRow(column(8,uiOutput("dir_fastq_text1")),
                    column(3,shinyFiles::shinyDirButton("dir_fastq",label="Browse...",title="Please select the directory containing FASTQ files:",class="btn-secondary"))),
           textOutput("dir_fastq_text2"),HTML("<br>"),
          fileInput("file.library",label=p("Please load an sgRNA library file:",downloadLink("q3",label="?",class="btn btn-primary btn-xs"))),
          sliderInput("spacer.range","Spacer location:",min=1,max=50,value=c(1,19),step=1),
           fluidRow(
             column(6,checkboxInput("spacer.rev",label="Reverse.",value=T)),
             column(6,checkboxInput("spacer.comp",label="Complement.",value=T)),
             textAreaInput("essentials",label=p("Add core-essential genes:")),
             textAreaInput("nonessentials",label=p("Add non-essential genes:")),
             textAreaInput("nontargeting",label=p("Add 'Gene' label of non-targeting sgRNAs:"))
          ),
          fluidRow(align="center",actionButton("button_run",label="RUN",class="btn-primary"))
        )),
      column(8,tabsetPanel(
        tabPanel("Charts",
               plotOutput("mapped2"),
               plotOutput("mapped1")
        ),
        tabPanel("File annotation",
           HTML("<br>"),uiOutput("text_anno1"),HTML("<br>"),
           fluidRow(align="center",actionButton("save_anno","Click here when finished.",class="btn-danger")),
           br(),
           rhandsontable::rHandsontableOutput("anno_table")),
        tabPanel("Table",
                 wellPanel(fluidRow(
                   column(6,radioButtons("fastq.level",label="Level:",choices=list("FASTQ files","Conditions"),selected="FASTQ files",inline=T)),
                   DT::dataTableOutput(outputId="data_input"))
                 ))
      )),
      tags$head(tags$style(".shiny-output-error{color: blue;}")),
      shinyBS::bsPopover("q1",content=HTML(paste("sgRNA  &nbsp;&nbsp;&nbsp;&nbsp; Gene (optional) &nbsp;&nbsp;&nbsp;&nbsp; Category (optional) &nbsp;&nbsp;&nbsp;&nbsp; sample1 &nbsp;&nbsp;&nbsp;&nbsp; sample2 &nbsp;&nbsp;&nbsp;&nbsp;",
                                                 "sgrna1 &nbsp;&nbsp;&nbsp;&nbsp; gene1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; category1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .............. &nbsp;&nbsp;&nbsp;&nbsp; ..............",
                                                 "sgrna2 &nbsp;&nbsp;&nbsp;&nbsp; gene2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; category2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .............. &nbsp;&nbsp;&nbsp;&nbsp; ..............",
                                                 "sgrna3 &nbsp;&nbsp;&nbsp;&nbsp; gene3 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; category3 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; .............. &nbsp;&nbsp;&nbsp;&nbsp; ..............",
                                                 sep="<br/>")),title="Format of table of counts (click to download example).",placement="right",options=list(container="body")),
      shinyBS::bsPopover("q3",content=HTML(paste("sgRNA  &nbsp;&nbsp;&nbsp;      Gene &nbsp;&nbsp;&nbsp;&nbsp;      Category (optional)",
                                                 "sgrna1 &nbsp;&nbsp;&nbsp;      gene1 &nbsp;&nbsp;&nbsp;      category1",
                                                 "sgrna2 &nbsp;&nbsp;&nbsp;      gene2 &nbsp;&nbsp;&nbsp;      category2",
                                                 "sgrna3 &nbsp;&nbsp;&nbsp;      gene3 &nbsp;&nbsp;&nbsp;      category3",
                                                 sep="<br/>")),title="sgRNA Library Format (click to download example)",placement="right",options=list(container="body")),
      shinyBS::bsPopover("q4",content=HTML(paste("sgRNA  &nbsp;&nbsp;&nbsp;      Gene &nbsp;&nbsp;&nbsp;&nbsp;      Category (optional)",
                                                 "sgrna1 &nbsp;&nbsp;&nbsp;      gene1 &nbsp;&nbsp;&nbsp;      category1",
                                                 "sgrna2 &nbsp;&nbsp;&nbsp;      gene2 &nbsp;&nbsp;&nbsp;      category2",
                                                 "sgrna3 &nbsp;&nbsp;&nbsp;      gene3 &nbsp;&nbsp;&nbsp;      category3",
                                                 sep="<br/>")),title="sgRNA Library Format (click to download example)",placement="right",options=list(container="body"))
   ),
   tabPanel("Quality control",
          column(4,
                 h3("Objectives"),
                 tags$ul(
                   tags$li("To visualize the quality of the experiment."),
                   tags$li("To annotate the control samples for each condition."),
                   tags$li("To compute the sgRNA log-fold-change through a default or re-optimized pseudocount.")
                 ),
                 h3("Instructions"),
                 tags$ol(
                   tags$li("If you skipped the FASTQ Alignment module:"),
                   tags$ul(
                     tags$li("Select a working directory in which all output is to be stored."),
                     tags$li("Upload a table of read counts, which must include a column named 'sgRNA'."),
                     tags$li("If the table of read counts is missing the 'Gene' column, please also upload an sgRNA library file containing the following two columns: 'sgRNA' and 'Gene'."),
                     tags$li("If the library includes gene categories, please add a column called 'Category' to either the table of read counts or the library file, containing the following identifiers: 'Essential', 'Nonessential', 'Nontargeting' or 'Other'."),
                     tags$li("Click on RUN.")
                   ),
                   tags$li("Otherwise, proceed directly to the 'Condition annotation' tab and follow the instructions shown.")# specify which is the control sample (i.e. 'Day 0') for each condition in the table provided. Finally, click on the red button above the table to save it.")
                 ),
            wellPanel(
            tags$style(".popover{max-width: 100%;}"),
            fluidRow(column(8,uiOutput("dir_qc_text1")),
                     column(3,shinyFiles::shinyDirButton("dir_qc",label="Browse...",title="Select working directory",class="btn-secondary"))),
            textOutput("dir_qc_text2"),HTML("<br>"),
            conditionalPanel(
              condition="output.fastq == 0",
              fileInput("file.counts",label=p("Read counts file:",downloadLink("q1",label="?",class="btn btn-primary btn-xs"))),
              fileInput("file.library2",label=p("sgRNA library file:",downloadLink("q4",label="?",class="btn btn-primary btn-xs"))),
              textAreaInput("essentials1",label=p("Add core-essential genes:")),
              textAreaInput("nonessentials1",label=p("Add non-essential genes:")),
              textAreaInput("nontargeting1",label=p("Add 'Gene' label of non-targeting sgRNAs:"))
            ),
            conditionalPanel(
              condition="output.fastq != 0",
              uiOutput("counts_qc"),HTML("<br>"),
              uiOutput("library_qc"),HTML("<br>")
            ),
            conditionalPanel(
              condition="output.fastq == 0",
              fluidRow(align="center",actionButton("button_qc",label="RUN",class="btn-primary"))
            )
          )),
        column(8,tabsetPanel(
          tabPanel("Charts",
                   plotOutput("classification"),
                   plotOutput("sizes"),
                   plotOutput("depths"),
                   plotOutput("zeros"),
                   plotOutput("gini")
          ),
          tabPanel("Condition annotation",
                   HTML("<br>"),uiOutput("text_anno2"),HTML("<br>"),
                   tags$ul(
                     tags$li("If you wish to re-optimize the pseudocount, uncheck the corresponding box (not recommended if the amount of control genes is small)."),
                     tags$li("If you wish to generate a PDF report of the quality control (and pandoc is installed), check the corresponding box.")
                   ),
                   wellPanel(fluidRow(
                    column(4,checkboxInput("pseudo_default",label="Use default pseudo-count = 0.15.",value=TRUE),
                    conditionalPanel(condition="!input.pseudo_default",
                                    checkboxInput("pseudo_optimize",label="Optimize pseudo-count.",value=FALSE),
                                    conditionalPanel(condition="!input.pseudo_optimize",
                                                     sliderInput("pseudo.count",label="Choose pseudo-count:",value=0.15,min=0.01,max=0.3,step=0.01)))),
                    column(6,checkboxInput("qc.report",label="Generate PDF report (requires pandoc installation).",value=FALSE))
                   )),
                   fluidRow(align="center",actionButton("save_anno1","Click here when finished.",class="btn-danger")),
                   br(),
                   rhandsontable::rHandsontableOutput("anno1_table")),
          tabPanel("Distributions",
                   wellPanel(
                     selectInput("lfc.condition",choices="",label="Condition:"),
                     fluidRow(
                       column(4,radioButtons("lfc.style",label="Visualization:",choices=c("Box plot","Density plot"),selected="Box plot",inline=T)),
                       column(4,radioButtons("folds.level",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=T)),
                       column(4,radioButtons("log.lfc",label="Stage:",choices=c("Log-read-count","Log-fold-change"),selected="Log-fold-change",inline=T))
                     )),
                   plotOutput("folds",width="100%",height="auto"),
                   plotOutput("ssmd")
          ),
          tabPanel("Correlation/PCA",
                   wellPanel(fluidRow(
                     column(4,radioButtons("corr_pca",label="Analysis:",choices=list("Correlation","Principal component analysis"),selected="Correlation",inline=F)),
                     column(4,radioButtons("corr_pca_method",label="Method:",choices=c("pearson","spearman"),selected="pearson",inline=F)),
                     column(4,radioButtons("corr.level",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=F))
                   )),
                   conditionalPanel(condition="input.corr_pca=='Correlation'",
                                    plotOutput("correlation",width="100%",height="auto")),
                   conditionalPanel(condition="input.corr_pca=='Principal component analysis'",
                                    plotOutput("pca1",width="100%",height="auto"),
                                    plotOutput("pve"))
          ),
          tabPanel("Reproducibility",
                   wellPanel(selectInput("qc.condition",choices="",label="Condition:"),
                             fluidRow(
                               column(6,radioButtons("reproducibility.level",label="Level:",choices=c("sgRNA","Gene"),selected="Gene",inline=T)),
                               column(6,radioButtons("reproducibility.stage",label="Stage:",choices=c("Log-read-count","Log-fold-change"),selected="Log-fold-change",inline=T))
                             )),
                   plotOutput("reproducibility",width="100%",height="auto")
          ),
          tabPanel("Table",
                   wellPanel(fluidRow(
                     column(4,radioButtons("qc.level",label="Level:",choices=list("sgRNA","Gene"),selected="sgRNA",inline=T)),
                     column(4,radioButtons("qc.stage",label="Stage:",choices=c("Log-read-count","Log-fold-change"),selected="Log-fold-change",inline=T)),
                     column(4,checkboxInput("qc.averaged",label="Replicate-averaged.",value=T))
                   )),
                   DT::dataTableOutput(outputId="data_qc")
          )
        ))),
   tabPanel("Gene essentiality measure",
        column(4,
               h3("Objectives"),
               tags$ul(
                 tags$li("To perform a 2-tail a-RRA analysis of gene essentiality. Within MoPAC's framework, this is necessary if there are no control genes available for the normalization of two conditions when computing the differential essentiality scores, or if the control genes should be filtered when computing the essentiality scores."),
                 tags$li("To compute a condition-specific measure of gene essentiality score through either default or re-optimized sgRNA weights.")
               ),
               h3("Instructions"),
               tags$ol(
                 tags$li("Choose whether the sgRNA weights should be re-optimized (not recommended if the amount of control genes is small)."),
                 tags$li("If a subset of the samples are to be used, uncheck the corresponding box."),
                 tags$li("If you wish to modify the internal RRA settings, uncheck the corresponding box."),
                 tags$li("Click on RUN.")
               ),
         wellPanel(
            tags$style(".popover{max-width: 100%;}"),
            uiOutput("dir_ea_text1"),
            textOutput("dir_ea_text2"),HTML("<br>"),
            radioButtons("empirical.weight",label="Weight calculation:",choices=c("Averaged","Pre-trained","Self-trained"),selected="Averaged",inline=T),
            checkboxInput("use_reps",label="Use all samples for the analysis.",value=T),
            conditionalPanel(condition="!input.use_reps",
                             selectizeInput("use_reps1",label="Choose samples to proceed with the analysis:",choices=NULL,multiple=TRUE,selected=1,width="100%")),
            checkboxInput("rra_default",label="Use default 2-RRA settings.",value=TRUE),
            conditionalPanel(condition="!input.rra_default",
                             checkboxInput("read.filtered",label="Filter control genes with 2-tail RRA.",value=F),
                             sliderInput("pvalue",label="P value threshold:",value=0.05,min=0.01,max=0.1,step=0.01),
                             sliderInput("fraction",label="Minimum fraction of conditions satisfied:",value=0.8,min=0.1,max=1,step=0.1)),
            fluidRow(align="center",actionButton("button_ea",label="RUN",class="btn-primary")),
            uiOutput("inforra")
          )),
        column(8,tabsetPanel(
          tabPanel("Scores",
                   wellPanel(fluidRow(
                     column(8,selectInput("ea.condition",choices="",label="Condition:")),
                     column(4,radioButtons("ea.style",label="Visualization:",choices=c("Density plot","Box plot"),selected="Box plot",inline=F))
                   )),
                   plotOutput("scores",width="100%",height="auto")
          ),
         tabPanel("Filtering",
           wellPanel(
             selectInput("rra.condition",choices="",label="Visualize essentiality p values of condition:"),
             selectizeInput("rragene",label="Highlight genes:",choices=NULL,multiple=TRUE,selected=1,width="100%")
           ),

        fluidRow(
          column(6,div(style="position:relative",plotOutput("rra",width="100%",height="auto",hover=hoverOpts("rra_hover",delay=10,delayType="debounce"),click="clickrra",brush=brushOpts(id="plotrra_brush",resetOnNew=FALSE)),uiOutput("plotrrainfo"))),
          column(6,div(style="position:relative",plotOutput("rrazoom",width="100%",height="auto",hover=hoverOpts("rrazoom_hover",delay=10,delayType="debounce"),click="clickrrazoom"),uiOutput("plotrrazoominfo")))
         )),
          tabPanel("Table",
                   wellPanel(fluidRow(
                     column(6,checkboxInput("ea.averaged",label="Replicate-averaged.",value=T))
                   )),
                   DT::dataTableOutput(outputId="data_ea")
          ))
        )
),

tabPanel("Differential essentiality analysis",
              column(4,
                     h3("Objectives"),
                     tags$ul(
                       tags$li("To normalize the distribution of one pair of conditions."),
                       tags$li("To generate differential gene essentiality scores for the pair of conditions selected.")
                     ),
                     h3("Instructions"),
                     tags$ol(
                       tags$li("Select the two conditions to analyze."),
                       tags$li("Choose whether to use the RRA-generated list of nondifferential genes for the normalization, or a custom list of controls located in the sgRNA library file."),
                       tags$li("Choose whether you wish to generate an html report and have 'plotly' install."),
                       tags$li("Click on RUN.")
                     ),
                wellPanel(
                tags$style(".popover{max-width: 100%;}"),
                # fluidRow(column(8,uiOutput("dir_rra_text1")),
                #          column(3,shinyFiles::shinyDirButton("dir_rra",label="Browse...",title="Select working directory",class="btn-secondary"))),
                uiOutput("dir_dea_text1"),
                textOutput("dir_dea_text2"),HTML("<br>"),
                  selectInput("condition1",choices="",label="First condition:"),
                  selectInput("condition2",choices="",label="Second condition:"),
                  selectizeInput("dagene",label="Highlight genes:",choices=NULL,multiple=TRUE,selected=1,width="100%"),
                  radioButtons("normalization",label="Normalization based on:",choices=list("Unenriched","Controls"),selected="Unenriched",inline=T),
                  checkboxInput("da.report",label="Generate HTML report (requires pandoc and plotly installation).",value=FALSE),
                  fluidRow(align="center",actionButton("button_da",label="RUN",class="btn-primary")),
                  uiOutput("infoscore")
              )),
         column(8,tabsetPanel(
           tabPanel("Volcano",
                    wellPanel(radioButtons("p.volcano",label="P-value source:",choices=c("sgRNAs","Replicates"),selected="sgRNAs",inline=T)),
            fluidRow(
             column(6,div(style="position:relative",plotOutput("plot3",width="100%",height="auto",hover=hoverOpts("plot3_hover",delay=10,delayType="debounce"),click="click3",brush=brushOpts(id="plot3_brush",resetOnNew=FALSE)),uiOutput("plot3info"))),
             column(6,div(style="position:relative",plotOutput("plot4",width="100%",height="auto",hover=hoverOpts("plot4_hover",delay=10,delayType="debounce"),click="click4"),uiOutput("plot4info")))
           )),
           tabPanel("Scores",fluidRow(
             column(6,div(style="position:relative",plotOutput("plot1",width="100%",height="auto",hover=hoverOpts("plot1_hover",delay=10,delayType="debounce"),click="click1",brush=brushOpts(id="plot1_brush",resetOnNew=FALSE)),uiOutput("plot1info"))),
             column(6,div(style="position:relative",plotOutput("plot2",width="100%",height="auto",hover=hoverOpts("plot2_hover",delay=10,delayType="debounce"),click="click2"),uiOutput("plot2info")))),
             plotOutput("weights"),
             plotOutput("scaling"),
             plotOutput("shifting")
           ),
           tabPanel("sgRNA-LFC",fluidRow(
             column(6,div(style="position:relative",plotOutput("plot5",width="100%",height="auto",hover=hoverOpts("plot5_hover",delay=10,delayType="debounce"),click="click5",brush=brushOpts(id="plot5_brush",resetOnNew=FALSE)),uiOutput("plot5info"))),
             column(6,div(style="position:relative",plotOutput("plot6",width="100%",height="auto",hover=hoverOpts("plot6_hover",delay=10,delayType="debounce"),click="click6"),uiOutput("plot6info")))
           )),

           tabPanel("Table",
                    wellPanel(selectInput("data_da",label="Dataset:",choices=c("Normalized pair","Pre-normalized replicates","Pre-normalized conditions"),selected="Normalized pair")),
                    DT::dataTableOutput(outputId="table_da")
           )
         ))
         ),

tabPanel("Network analysis",
         column(4,
                h3("Objectives"),
                tags$ul(
                  tags$li("To search for the top-ranking genes in the STRINGdb database."),
                  tags$li("To perform network enrichment analysis for clusters in the top-ranking genes.")
                ),
                h3("Instructions"),
                tags$ol(
                  tags$li("Select one of the pair of conditions analyzed in the differential essentiality analysis module, and which should be differentially essential."),
                  tags$li("Choose a threshold for significance of gene differential essentiality and a gene clustering algorithm."),
                  tags$li("Click on RUN."),
                  tags$li("Select a network in which to search for enrichment and one of the clusters in the network plotted on the right."),
                  tags$li("Choose whether you wish to introduce a significance threshold for gene interactions in the STRING database."),
                  tags$li("Click on 'Find pathways'.")
                ),
                wellPanel(
                  tags$style(".popover{max-width: 100%;}"),
                  uiOutput("dir_net_text1"),
                  textOutput("dir_net_text2"),HTML("<br>"),
                  radioButtons("list_source",label="Source of significant gene list:",choices=list("Differential score","Custom list"),selected="Differential score",inline=T),
                  conditionalPanel(condition="input.list_source=='Differential score'",
                    selectInput("comparisons",choices="",label="Select condition comparison:"),
                    selectInput("comparisons_essential",choices="",label="Select differentially essential condition:"),
                    # selectInput("net_algorithm",choices=c("fastgreedy","walktrap","spinglass","edge.betweenness"),label="Gene clustering algorithm:",selected="fastgreedy"), #uncomment to add clustering
                    numericInput("comparisons_sd",label="Minimum std. devs. from median",value=3,step=0.1,min=0),
                    uiOutput("comparisons_number")
                  ),
                  # HTML("<br>"),
                  conditionalPanel(condition="input.list_source=='Custom list'",textAreaInput("gene.list",label="Custom gene list:",height="300px")),
                  sliderInput("comparisons_tresh",label="Minimum STRING interaction score",value=700,step=1,min=0,max=1000),
                  fluidRow(align="center",actionButton("button_net",label="RUN",class="btn-primary")),
                  uiOutput("infonode")
               ),
               uiOutput("infonet"),HTML("<br>"),
               wellPanel(
                 selectInput("net_database",label="Select a network database for enrichment analysis:",choices=c("Kyoto Encyclopedia of Genes and Genomes","Process","Component","Function","Pfam","InterPro","Tissue","Disease"),selected="Kyoto Encyclopedia of Genes and Genomes"),
                 # selectInput("net_cluster",label="Network cluster to analyze:",choices="All",selected="All"), #uncomment to add clustering
                 radioButtons("net_method",label="Network enrichment analysis method:",choices=list("STRINGdb","NEAT (slower)"),selected="STRINGdb",inline=T),
                 checkboxInput("net_default",label="Use full annotation list.",value=T),
                 conditionalPanel(condition="!input.net_default",
                                  sliderInput("net_minscore","Score threshold to filter annotations:",value=0,step=1,min=0,max=5),
                                  checkboxInput("net_iea",label="Use electronic-inferred annotations",value=T)),
                 fluidRow(align="center",actionButton("button_paths",label="Find pathways",class="btn-primary")),
                 HTML("<br><b>Sources of network analysis data and algorithms:</b><br><br>Szklarczyk, D. et al. The STRING database in 2017: Quality-controlled protein-protein association networks, made broadly accessible. Nucleic Acids Res. 45, D362–D368 (2017).<br><br>Signorelli, M. et al. NEAT: an efficient network enrichment analysis test. BMC Bioinformatics 17, 352 (2016).")
               )
         ),
         column(8,
            visNetwork::visNetworkOutput("network1",width="100%",height="1000px"),
            DT::dataTableOutput(outputId="table_net")
         )
)
)
)

# tags$style(type="text/css",".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; }")
