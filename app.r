library(shiny)
library(shinyWidgets)
library(shinyjs)
library(kableExtra)
library(knitr)
library(Rcpp)
library(fontawesome)
library(Rcpp)
library(BH)
library("bslib")

sourceCpp("Functions/pt.cpp", cacheDir = "tmp_cache")
source("Functions/onesample.r")
source("Functions/onesample_e.r")
source("Functions/twosample.r")
source("Functions/twosample_e.r")
source("Functions/Correlation.r")
source("Functions/Correlation_e.r")
source("Functions/ANOVA.r")
source("Functions/ANOVAe.r")
source("Functions/binomial.r")
source("Functions/binomial_e.r")

ui <- navbarPage(id = "id",
                 "\\(\\text{BayesPower}_{AlphaVersion}\\)",
  navbarMenu(
    "\\(\\text{ Standardized Mean Difference}\\)",
    
tabPanel("\\(\\text{One-sample/paired t-test}\\)",withMathJax(),
        sidebarLayout(sidebarPanel(
          prettyRadioButtons(
                 "Modet1",
                 "\\(\\text{Select Mode}\\)",
                 choices = list("\\(\\text{Sample size determination}\\)" = 1, 
                                "\\(\\text{Fixed N}\\)" = 2, 
                                "\\(\\text{BF calculator}\\)" = 3),
                 selected = 1,
                 inline = T
               )
,
          fluidRow(
            column(width = 6, prettyRadioButtons(
              inputId = "h0t1",
              label = em("\\(\\mathcal{H}_0:\\)"),  # Inline LaTeX for H0
              choices = list(
                "\\(\\delta = 0\\)" = 1,  # Inline LaTeX for delta = 0
                "\\(\\delta \\in \\{-\\epsilon, \\epsilon\\}\\)" = 2   # Inline LaTeX for delta = 1
              ),
              inline = T,
              selected = 1  # Default selected option
            ))
            ,
            
            column(width = 6,
                   conditionalPanel(
                     condition = "input.h0t1 == 1",  # Show when delta = 0 is selected
                     prettyRadioButtons(
                       inputId = "h1t1",
                       label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                       choices = list(
                         "\\(\\delta ≠ 0\\)" = 1,  # Inline LaTeX for delta ≠ 0
                         "\\(\\delta > 0\\)" = 2,   # Inline LaTeX for delta > 0
                         "\\(\\delta < 0\\)" = 3    # Inline LaTeX for delta < 0
                       ),
                       inline = T,
                       selected = 1  # Default selected option
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.h0t1 == 2",  # Show when delta subset {-ε, ε} is selected
                     prettyRadioButtons(
                       inputId = "h1t1e",
                       label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                       choices = list(
                         "\\(\\delta \\not\\in \\{-\\epsilon, \\epsilon\\}\\)" = 1,  # Inline LaTeX for delta ≠ {-ε, ε}
                         "\\(\\delta > \\epsilon\\)" = 2,   # Inline LaTeX for delta > ε
                         "\\(\\delta < \\epsilon\\)" = 3    # Inline LaTeX for delta < ε
                       ),
                       inline = T,
                       selected = 1  # Default selected option
                     )
                   )
            )
          ),
          fluidRow(
            column(width=6,
                   conditionalPanel( condition =  "input.h1t1e == 2 && input.h0t1 == 2",
                                     em( "\\( -\\epsilon = 0 \\)")),
                   conditionalPanel( condition =  "(input.h1t1e == 1 || input.h1t1e == 3) && input.h0t1 == 2",
                                     sliderInput(
                                       inputId = "lbt1e",
                                       label =  "\\( -\\epsilon \\)",
                                       min = -.5,
                                       ticks = FALSE,
                                       max = -0.01,
                                       value = -0.2,  # Default value
                                       step = 0.01    # Step size
                                     )
                                     )),
            column(width=6,
                   conditionalPanel( condition =  "input.h1t1e == 3 && input.h0t1 == 2",
                                     em( "\\( \\epsilon = 0 \\)")),
                   conditionalPanel( condition =  "(input.h1t1e == 1 || input.h1t1e == 2) && input.h0t1 == 2",
                                     sliderInput(
                                       inputId = "ubt1e",
                                       label = "\\(\\epsilon \\)",
                                       min = .01,
                                       max = .5,
                                       value = 0.2,  # Default value
                                       step = 0.01,   # Step size
                                       ticks = FALSE)
                   ))),
          prettyRadioButtons(
            inputId = "modelt1",
            label = em("\\(\\text{ Analysis Prior Distribution}\\)"),  # Inline LaTeX for H1
            choices = list(
              "\\( \\text{Scaled t} \\)" = 1,  # Inline LaTeX for delta ≠ 0
              "\\( \\text{Normal} \\)" = 2,   # Inline LaTeX for delta > 0
              "\\( \\text{Moment} \\)" = 3    # Inline LaTeX for delta < 0
            ),
            inline = T,
            selected = 1  # Default selected option
          ),
          fluidRow(
            column(width=4,
                   conditionalPanel( condition = "input.h0t1 == 1",
                                     sliderInput(
                                       inputId = "lt1",
                                       label =  "\\( \\text{Location} \\)",
                                       min = -2,
                                       ticks = FALSE,
                                       max = 2,
                                       value = 0,  # Default value
                                       step = 0.01    # Step size
                                     )
                   ),
                   conditionalPanel("input.h0t1 == 2",em("\\( \\text{Location = 0} \\)"))),
            column(width=4,
                  
                                     sliderInput(
                                       inputId = "st1",
                                       label =  "\\( \\text{Scale} \\)",
                                       min = .01,
                                       ticks = FALSE,
                                       max = 2,
                                       value = .707,  # Default value
                                       step = 0.001    # Step size
                                     )
                   ),
            column(width=4,
                   conditionalPanel( condition = "input.modelt1 == 1",
                   sliderInput(
                     inputId = "dft1",
                     label =  "\\( \\text{df} \\)",
                     min = 1,
                     ticks = FALSE,
                     max = 100,
                     value = 1,  # Default value
                     step = 1    # Step size
                   )))),conditionalPanel(condition = "input.Modet1 == 1|input.Modet1 == 2",
          prettyRadioButtons(
            "prior",
            "\\( \\text{Design prior is the same as analysis prior: } \\)",
            choices = list("\\( \\text{Yes} \\)" = 1, 
                           "\\( \\text{No} \\)"  = 2),
            selected = 1,
            inline = T
          ),conditionalPanel( condition = "input.prior == 2",
                              prettyRadioButtons(
                                inputId = "modelt1d",
                                label = em("\\( \\text{Design prior distribution} \\)"),  # Inline LaTeX for H1
                                choices = list(
                                  "\\( \\text{Scaled t} \\)" = 1,  # Inline LaTeX for delta ≠ 0
                                  "\\( \\text{Normal} \\)" = 2,   # Inline LaTeX for delta > 0
                                  "\\( \\text{Moment} \\)" = 3,    # Inline LaTeX for delta < 0
                                  "\\( \\text{Point} \\)" = 4
                                ),
                                inline = T,
                                selected = 1  # Default selected option
                              ),
                              fluidRow(
                                column(width=4,
                                       conditionalPanel( condition = "input.h0t1 == 1|input.modelt1d==4",
                                                         sliderInput(
                                                           inputId = "lt1d",
                                                           label =  "\\( \\text{Location} \\)",
                                                           min = -2,
                                                           ticks = FALSE,
                                                           max = 2,
                                                           value = 0,  # Default value
                                                           step = 0.01    # Step size
                                                         )),
                                       conditionalPanel(condition = "input.h0t1 == 2 && input.modelt1d != 4",em("\\( \\text{Location = 0} \\)"))
                                       ),
                                column(width=4,
                                       conditionalPanel(condition = "input.modelt1d == 1 || input.modelt1d == 2 || input.modelt1d == 3",
                                       sliderInput(
                                         inputId = "st1d",
                                         label =  "\\( \\text{Scale} \\)",
                                         min = .01,
                                         ticks = FALSE,
                                         max = 2,
                                         value = .707,  # Default value
                                         step = 0.001    # Step size
                                       ))
                                ),
                                column(width=4,
                                       conditionalPanel( condition = "input.modelt1d == 1",
                                                         sliderInput(
                                                           inputId = "dft1d",
                                                           label =  "\\( \\text{df} \\)",
                                                           min = 1,
                                                           ticks = FALSE,
                                                           max = 100,
                                                           value = 1,  # Default value
                                                           step = 1    # Step size
                                                         ))))
                              
                              
                              
                              
                              )),conditionalPanel("input.modelt1d == 4 &(input.h1t1 == 2 |input.h1t1 == 3|input.h1t1e == 2|input.h1t1 == 3)",
                                                  em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
                              ),conditionalPanel(condition = "input.Modet1 == 1",
                             em("\\(\\text{Controlloing for the probability of}\\)"),
                             
                             fluidRow(
                               column(width = 6,sliderInput("powert1", "\\( \\text{True Positive Evidence:} \\)",
                                                            min =.5, max = .99, value = .8,step = .01,ticks = FALSE)),
                               column(width = 6,
                                      conditionalPanel("input.h0t1 == 1|input.h0t1 == 2",
                                      sliderInput("alphat1", "\\( \\text{False Positive Evidence:} \\)",
                                                            min =.001, max = .05, value = .05,step = .001,ticks = FALSE)
                                      
                                      ))),
                             
                             
                             
                             
                             
                             
                             
          ), conditionalPanel(condition = "input.Modet1 == 1|input.Modet1 == 2",sliderInput("bt1", "\\( \\text{Bound of compelling evidence:} \\)", min=1,max = 20,value = 3,ticks = FALSE)),
conditionalPanel(condition = "input.Modet1 == 2",
                 numericInput(inputId ="nt1",
  label="\\( \\text{Sample Size:} \\)",
  value = 50
  
)),


 conditionalPanel(condition = "input.Modet1 == 1|input.Modet1 == 2",
                           actionButton("runt1", label = "\\( \\text{Run} \\)"),
                  conditionalPanel(condition ="input.Modet1 == 1",
                                  em(span("\\(\\text{Note: Error when the required N > 100,000} \\)", style = "color: red;")))),



conditionalPanel(condition = "input.Modet1 == 3",
                 fluidRow(
                   column(width = 6,numericInput(
                     inputId = "t1df",
                     label="\\( \\text{Degree of freedom:} \\)",
                     value = 50
                     
                   )),
                   column(width = 6,numericInput(
                     inputId = "t1tval",
                     label="\\( \\text{t-value:} \\)",
                     value = 2
                     
                   ))
                   
                 ),
                 actionButton("cal1", label = "\\( \\text{Calculate} \\)"),
                 htmlOutput("BFt1"))
          
               
             ),
          mainPanel(
            navset_card_underline(
              nav_panel(id = "mt1",em("\\(\\text{Result}\\)"),
                        fluidRow(
                          column(6, plotOutput("priort1")),
                          column(6,   htmlOutput("resultt1"))
                        ),
                        plotOutput("bfrt1")),
             nav_panel(id = "pt1",em("\\(\\text{Power Curve}\\)"),
                       plotOutput("PCt1"))
                        
              ))
    
             )
   ),tabPanel("\\(\\text{Independent samples t-test}\\)",withMathJax(),sidebarLayout(sidebarPanel(
     prettyRadioButtons(
       "Modet2",
       "\\(\\text{Select Mode}\\)",
       choices = list("\\(\\text{Sample size determination}\\)" = 1, 
                      "\\(\\text{Fixed N}\\)" = 2, 
                      "\\(\\text{BF calculator}\\)" = 3),
       selected = 1,
       inline = T
     )
     ,
     fluidRow(
       column(width = 6, prettyRadioButtons(
         inputId = "h0t2",
         label = em("\\(\\mathcal{H}_0:\\)"),  # Inline LaTeX for H0
         choices = list(
           "\\(\\delta = 0\\)" = 1,  # Inline LaTeX for delta = 0
           "\\(\\delta \\in \\{-\\epsilon, \\epsilon\\}\\)" = 2   # Inline LaTeX for delta = 1
         ),
         inline = T,
         selected = 1  # Default selected option
       ))
       ,
       
       column(width = 6,
              conditionalPanel(
                condition = "input.h0t2 == 1",  # Show when delta = 0 is selected
                prettyRadioButtons(
                  inputId = "h1t2",
                  label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                  choices = list(
                    "\\(\\delta ≠ 0\\)" = 1,  # Inline LaTeX for delta ≠ 0
                    "\\(\\delta > 0\\)" = 2,   # Inline LaTeX for delta > 0
                    "\\(\\delta < 0\\)" = 3    # Inline LaTeX for delta < 0
                  ),
                  inline = T,
                  selected = 1  # Default selected option
                )
              ),
              
              conditionalPanel(
                condition = "input.h0t2 == 2",  # Show when delta subset {-ε, ε} is selected
                prettyRadioButtons(
                  inputId = "h1t2e",
                  label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                  choices = list(
                    "\\(\\delta \\not\\in \\{-\\epsilon, \\epsilon\\}\\)" = 1,  # Inline LaTeX for delta ≠ {-ε, ε}
                    "\\(\\delta > \\epsilon\\)" = 2,   # Inline LaTeX for delta > ε
                    "\\(\\delta < \\epsilon\\)" = 3    # Inline LaTeX for delta < ε
                  ),
                  inline = T,
                  selected = 1  # Default selected option
                )
              )
       )
     ),
     fluidRow(
       column(width=6,
              conditionalPanel( condition =  "input.h1t2e == 2 && input.h0t2 == 2",
                                em( "\\( -\\epsilon = 0 \\)")),
              conditionalPanel( condition =  "(input.h1t2e == 1 || input.h1t2e == 3) && input.h0t2 == 2",
                                sliderInput(
                                  inputId = "lbt2e",
                                  label =  "\\( -\\epsilon \\)",
                                  min = -.5,
                                  ticks = FALSE,
                                  max = -0.01,
                                  value = -0.2,  # Default value
                                  step = 0.01    # Step size
                                )
              )),
       column(width=6,
              conditionalPanel( condition =  "input.h1t2e == 3 && input.h0t2 == 2",
                                em( "\\( \\epsilon = 0 \\)")),
              conditionalPanel( condition =  "(input.h1t2e == 1 || input.h1t2e == 2) && input.h0t2 == 2",
                                sliderInput(
                                  inputId = "ubt2e",
                                  label = "\\(\\epsilon \\)",
                                  min = .01,
                                  max = .5,
                                  value = 0.2,  # Default value
                                  step = 0.01,   # Step size
                                  ticks = FALSE)
              ))),
     prettyRadioButtons(
       inputId = "modelt2",
       label = em("\\(\\text{ Analysis Prior Distribution}\\)"),  # Inline LaTeX for H1
       choices = list(
         "\\( \\text{Scaled t} \\)" = 1,  # Inline LaTeX for delta ≠ 0
         "\\( \\text{Normal} \\)" = 2,   # Inline LaTeX for delta > 0
         "\\( \\text{Moment} \\)" = 3    # Inline LaTeX for delta < 0
       ),
       inline = T,
       selected = 1  # Default selected option
     ),
     fluidRow(
       column(width=4,
              conditionalPanel( condition = "input.h0t2 == 1",
                                sliderInput(
                                  inputId = "lt2",
                                  label =  "\\( \\text{Location} \\)",
                                  min = -2,
                                  ticks = FALSE,
                                  max = 2,
                                  value = 0,  # Default value
                                  step = 0.01    # Step size
                                )
              ),
              conditionalPanel("input.h0t2 == 2",em("\\( \\text{Location = 0} \\)"))),
       column(width=4,
              
              sliderInput(
                inputId = "st2",
                label =  "\\( \\text{Scale} \\)",
                min = .01,
                ticks = FALSE,
                max = 2,
                value = 1,  # Default value
                step = 0.01    # Step size
              )
       ),
       column(width=4,
              conditionalPanel( condition = "input.modelt2 == 1",
                                sliderInput(
                                  inputId = "dft2",
                                  label =  "\\( \\text{df} \\)",
                                  min = 1,
                                  ticks = FALSE,
                                  max = 100,
                                  value = 1,  # Default value
                                  step = 1    # Step size
                                )))),conditionalPanel(condition = "input.Modet2 == 1|input.Modet2 == 2",
                                                      prettyRadioButtons(
                                                        "priort2",
                                                        "\\( \\text{Design prior is the same as analysis prior: } \\)",
                                                        choices = list("\\( \\text{Yes} \\)" = 1, 
                                                                       "\\( \\text{No} \\)"  = 2),
                                                        selected = 1,
                                                        inline = T
                                                      ),conditionalPanel( condition = "input.priort2 == 2",
                                                                          prettyRadioButtons(
                                                                            inputId = "modelt2d",
                                                                            label = em("\\( \\text{Design prior distribution} \\)"),  # Inline LaTeX for H1
                                                                            choices = list(
                                                                              "\\( \\text{Scaled t} \\)" = 1,  # Inline LaTeX for delta ≠ 0
                                                                              "\\( \\text{Normal} \\)" = 2,   # Inline LaTeX for delta > 0
                                                                              "\\( \\text{Moment} \\)" = 3,    # Inline LaTeX for delta < 0
                                                                              "\\( \\text{Point} \\)" = 4
                                                                            ),
                                                                            inline = T,
                                                                            selected = 1  # Default selected option
                                                                          ),
                                                                          fluidRow(
                                                                            column(width=4,
                                                                                   conditionalPanel( condition = "input.h0t2 == 1|input.modelt2d==4",
                                                                                                     sliderInput(
                                                                                                       inputId = "lt2d",
                                                                                                       label =  "\\( \\text{Location} \\)",
                                                                                                       min = -2,
                                                                                                       ticks = FALSE,
                                                                                                       max = 2,
                                                                                                       value = 0,  # Default value
                                                                                                       step = 0.01    # Step size
                                                                                                     )),
                                                                                   conditionalPanel(condition = "input.h0t2 == 2 && input.modelt2d != 4",em("\\( \\text{Location = 0} \\)"))
                                                                            ),
                                                                            column(width=4,
                                                                                   conditionalPanel(condition = "input.modelt2d == 1 || input.modelt2d == 2 || input.modelt2d == 3",
                                                                                                    sliderInput(
                                                                                                      inputId = "st2d",
                                                                                                      label =  "\\( \\text{Scale} \\)",
                                                                                                      min = .01,
                                                                                                      ticks = FALSE,
                                                                                                      max = 2,
                                                                                                      value = 1,  # Default value
                                                                                                      step = 0.01    # Step size
                                                                                                    ))
                                                                            ),
                                                                            column(width=4,
                                                                                   conditionalPanel( condition = "input.modelt2d == 1",
                                                                                                     sliderInput(
                                                                                                       inputId = "dft2d",
                                                                                                       label =  "\\( \\text{df} \\)",
                                                                                                       min = 1,
                                                                                                       ticks = FALSE,
                                                                                                       max = 100,
                                                                                                       value = 1,  # Default value
                                                                                                       step = 1    # Step size
                                                                                                     ))))
                                                                          
                                                                          
                                                                          
                                                                          
                                                      )),conditionalPanel("input.modelt2d == 4 &(input.h1t2 == 2 |input.h1t2 == 3|input.h1t2e == 2|input.h1t2 == 3)",
                                                                          em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
                                                      ),conditionalPanel(condition = "input.Modet2 == 1",
                                                                                                                                                                              em("\\(\\text{Controlloing for the probability of}\\)"),
                                                                                                                                                                              
                                                                                                                                                                              fluidRow(
                                                                                                                                                                                column(width = 6,sliderInput("powert2", "\\( \\text{True Positive Evidence:} \\)",
                                                                                                                                                                                                             min =.5, max = .99, value = .8,step = .01,ticks = FALSE)),
                                                                                                                                                                                column(width = 6,
                                                                                                                                                                                       conditionalPanel("input.h0t2 == 1|input.h0t2 == 2",
                                                                                                                                                                                                        sliderInput("alphat2", "\\( \\text{False Positive Evidence:} \\)",
                                                                                                                                                                                                                    min =.001, max = .05, value = .05,step = .001,ticks = FALSE)
                                                                                                                                                                                                        
                                                                                                                                                                                       ))),
                                                                                                                                                                              
                                                                                                                                                                              
                                                                                                                                                                              
                                                                                                                                                                              
                                                                                                                                                                              
                                                                                                                                                                              
                                                                                                                                                                              
                                                                          ), 
     conditionalPanel(condition = "input.Modet2 == 1|input.Modet2 == 2",sliderInput("bt2", "\\( \\text{Bound of compelling evidence:} \\)", min=1,max = 20,value = 3,ticks = FALSE)),
     conditionalPanel(condition = "input.Modet2 == 1|input.Modet2 == 3",sliderInput("rt2", "\\( \\text{The ratio of sample sizes in two groups } N_2/N_1 \\)", min=1,max = 10,value = 1,ticks = FALSE)),
     conditionalPanel(condition = "input.Modet2 == 2",
                      fluidRow(
                        column(width=4,
                      numericInput(inputId ="n1t2",
                                   label="\\( \\text{Sample Size } N_1: \\)",
                                   value = 50)),
                      column(width=4,numericInput(inputId ="n2t2",
                                   label="\\( \\text{Sample Size } N_2: \\)",
                                   value = 50)))
                      
                      
                      
                      ),
     
     
     conditionalPanel(condition = "input.Modet2 == 1|input.Modet2 == 2",
                      actionButton("runt2", label = "\\( \\text{Run} \\)"),
                      conditionalPanel(condition ="input.Modet2 == 1",
                                       em(span("\\(\\text{Note: Error when the required N > 100,000} \\)", style = "color: red;")))
                      
                      
                      
                      ),
     
     
     
     conditionalPanel(condition = "input.Modet2 == 3",
                      fluidRow(
                        column(width = 6,numericInput(
                          inputId = "t2df",
                          label="\\( \\text{Degree of freedom:} \\)",
                          value = 50
                          
                        )),
                        column(width = 6,numericInput(
                          inputId = "t2tval",
                          label="\\( \\text{t-value:} \\)",
                          value = 2
                          
                        ))
                        
                      ),
                      actionButton("cal1", label = "\\( \\text{Calculate} \\)"),
                      htmlOutput("BFt2"))
     
     
   ),
   mainPanel(
     navset_card_underline(
       nav_panel(id = "mt2",em("\\(\\text{Result}\\)"),
                 fluidRow(
                   column(6, plotOutput("priort2")),
                   column(6,   htmlOutput("resultt2"))
                 ),
                 plotOutput("bfrt2")),
       nav_panel(id = "pt2",em("\\(\\text{Power Curve}\\)"),
                 plotOutput("PCt2"))
       
     ))
   
   )
              
              
              
              
              
              
              )




),
  tabPanel("\\(\\text{Correlation}\\)",withMathJax(),
           sidebarLayout(sidebarPanel(
             prettyRadioButtons(
               "Moder",
               "\\(\\text{Select Mode}\\)",
               choices = list("\\(\\text{Sample size determination}\\)" = 1, 
                              "\\(\\text{Fixed N}\\)" = 2, 
                              "\\(\\text{BF calculator}\\)" = 3),
               selected = 1,
               inline = T
             )
             ,
             fluidRow(
               column(width = 6, prettyRadioButtons(
                 inputId = "h0r",
                 label = em("\\(\\mathcal{H}_0:\\)"),  # Inline LaTeX for H0
                 choices = list(
                   "\\(\\rho = \\rho_0\\)"  = 1,  
                   "\\(\\rho \\in (\\rho_0 - \\epsilon, \\rho_0 + \\epsilon)\\)" = 2  
                 ),
                 inline = T,
                 selected = 1  # Default selected option
               ))
               ,
               
               column(width = 6,
                      conditionalPanel(
                        condition = "input.h0r == 1",  # Show when delta = 0 is selected
                        prettyRadioButtons(
                          inputId = "h1r",
                          label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                          choices = list(
                            "\\(\\rho ≠ \\rho_0\\)" = 1,  # Inline LaTeX for delta ≠ 0
                            "\\(\\rho > \\rho_0\\)" = 2,   # Inline LaTeX for delta > 0
                            "\\(\\rho < \\rho_0\\)" = 3    # Inline LaTeX for delta < 0
                          ),
                          inline = T,
                          selected = 1  # Default selected option
                        )
                      ),
                      
                      conditionalPanel(
                        condition = "input.h0r == 2",  # Show when delta subset {-ε, ε} is selected
                        prettyRadioButtons(
                          inputId = "h1re",
                          label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                          choices = list(
                            "\\(\\rho \\not\\in \\{\\rho_0 -\\epsilon, \\rho_0 +\\epsilon\\}\\)" = 1,  # Inline LaTeX for delta ≠ {-ε, ε}
                            "\\(\\rho > \\rho_0 + \\epsilon\\)" = 2,   # Inline LaTeX for delta > ε
                            "\\(\\rho < \\rho_0 - \\epsilon\\)" = 3    # Inline LaTeX for delta < ε
                          ),
                          inline = T,
                          selected = 1  # Default selected option
                        )
                      )
               )
             ),
             fluidRow(
               column(width=4,sliderInput(
                 inputId = "h0pho",
                 label =  "\\(\\rho_0\\)",
                 min = -1,
                 ticks = FALSE,
                 max = 1,
                 value = 0,  # Default value
                 step = 0.01    # Step size
               )),
                 column(width=4,
                      conditionalPanel( condition =  "(input.h1re == 1 || input.h1re == 3) && input.h0r == 2",
                                        sliderInput(
                                          inputId = "lbre",
                                          label =  "\\( -\\epsilon \\)",
                                          min = -.5,
                                          ticks = FALSE,
                                          max = -0.01,
                                          value = -0.2,  # Default value
                                          step = 0.01    # Step size
                                        ),htmlOutput("r_lower")
                      )),
               column(width=4,
                      conditionalPanel( condition =  "(input.h1re == 1 || input.h1re == 2) && input.h0r == 2",
                                        sliderInput(
                                          inputId = "ubre",
                                          label = "\\(\\epsilon \\)",
                                          min = .01,
                                          max = .5,
                                          value = 0.2,  # Default value
                                          step = 0.01,   # Step size
                                          ticks = FALSE),htmlOutput("r_upper")
                      ))),
             prettyRadioButtons(
               inputId = "modelr",
               label = em("\\(\\text{ Analysis Prior Distribution}\\)"),  # Inline LaTeX for H1
               choices = list(
                 "\\( \\text{Default Stretched Beta} \\)" = 1,  
                 "\\( \\text{Stretched Beta} \\)" = 2,  
                 "\\( \\text{Moment} \\)" = 3  
               ),
               inline = T,
               selected = 1  # Default selected option
             ),
             fluidRow(
               column(width=4,
                      conditionalPanel( condition = "input.modelr == 1",
                                        sliderInput(
                                          inputId = "kr",
                                          label =  "\\(k \\)",
                                          min = .01,
                                          ticks = FALSE,
                                          max = 10,
                                          value = 1,  # Default value
                                          step = 0.01    # Step size
                                        )
                      ),conditionalPanel( condition = "input.modelr == 3",
                                          sliderInput(
                                            inputId = "sr",
                                            label =  "\\(Scale \\)",
                                            min = 0.01,
                                            ticks = FALSE,
                                            max = 1,
                                            value = .01,  # Default value
                                            step = 0.01    # Step size
                                          )
                      )
                      
                      
                      ),
               column(width=4,
                      conditionalPanel( condition = "input.modelr == 2",
                                                sliderInput(
                                                  inputId = "ralpha",
                                                  label =  "\\(\\alpha \\)",
                                                  min = 0.01,
                                                  ticks = FALSE,
                                                  max = 10,
                                                  value = 1,  # Default value
                                                  step = 0.01    # Step size
                                                )
               )
               ),
               column(width=4,
                      conditionalPanel( condition = "input.modelr == 2",
                                        sliderInput(
                                          inputId = "rbeta",
                                          label =  "\\(\\beta \\)",
                                          min = 0.01,
                                          ticks = FALSE,
                                          max = 10,
                                          value = 1,  # Default value
                                          step = 0.01    # Step size
                                        )
                      )
                      
                      )
               
               ),conditionalPanel(condition = "input.Moder == 1|input.Moder == 2",
                                                              prettyRadioButtons(
                                                                "priorr",
                                                                "\\( \\text{Design prior is the same as analysis prior: } \\)",
                                                                choices = list("\\( \\text{Yes} \\)" = 1, 
                                                                               "\\( \\text{No} \\)"  = 2),
                                                                selected = 1,
                                                                inline = T
                                                              ),conditionalPanel( condition = "input.priorr == 2",
                                                                                  prettyRadioButtons(
                                                                                    inputId = "modelrd",
                                                                                    label = em("\\( \\text{Design prior distribution} \\)"),  # Inline LaTeX for H1
                                                                                    choices = list(
                                                                                      "\\( \\text{Default Stretched Beta} \\)" = 1,  
                                                                                      "\\( \\text{Stretched Beta} \\)" = 2,  
                                                                                      "\\( \\text{Moment} \\)" = 3 ,
                                                                                      "\\( \\text{Point} \\)" = 4
                                                                                    ),
                                                                                    inline = T,
                                                                                    selected = 1  # Default selected option
                                                                                  ),conditionalPanel(condition = "input.modelrd == 4",
                                                                                                     sliderInput(
                                                                                                       inputId = "h0phod",
                                                                                                       label =  "\\(\\rho_1\\)",
                                                                                                       min = -1,
                                                                                                       ticks = FALSE,
                                                                                                       max = 1,
                                                                                                       value = .3,  # Default value
                                                                                                       step = 0.01    # Step size
                                                                                                     ))
                                                                                  ,
                                                                                  fluidRow(
                                                                                    column(width=4,
                                                                                           conditionalPanel( condition = "input.modelrd == 1",
                                                                                                             sliderInput(
                                                                                                               inputId = "rkd",
                                                                                                               label =   "\\( k \\)",
                                                                                                               min = 0.01,
                                                                                                               ticks = FALSE,
                                                                                                               max = 10,
                                                                                                               value = 1,  # Default value
                                                                                                               step = 0.01    # Step size
                                                                                                             )),
                                                                                           
                                                                                           conditionalPanel( condition = "input.modelrd == 3",
                                                                                                             sliderInput(
                                                                                                               inputId = "rsd",
                                                                                                               label =   "\\( Scale \\)",
                                                                                                               min = 0.01,
                                                                                                               ticks = FALSE,
                                                                                                               max = 1,
                                                                                                               value = .1,  # Default value
                                                                                                               step = 0.01    # Step size
                                                                                                             ))
                                                                                           
                                                                                           
                                                                                    ),
                                                                                    column(width=4,
                                                                                           conditionalPanel(condition = "input.modelrd == 2",
                                                                                                            sliderInput(
                                                                                                              inputId = "ralphad",
                                                                                                              label =  "\\(\\alpha \\)",
                                                                                                              min = 0.01,
                                                                                                              ticks = FALSE,
                                                                                                              max = 10,
                                                                                                              value = .1,  # Default value
                                                                                                              step = 0.01    # Step size
                                                                                                            ))
                                                                                    ),
                                                                                    column(width=4,
                                                                                           conditionalPanel( condition = "input.modelrd == 2",
                                                                                                             sliderInput(
                                                                                                               inputId = "rbetad",
                                                                                                               label =  "\\(\\beta \\)",
                                                                                                               min = .01,
                                                                                                               ticks = FALSE,
                                                                                                               max = 10,
                                                                                                               value = .1,  # Default value
                                                                                                               step = 0.01    # Step size
                                                                                                             ))))
                                                                                  
                                                                                  
                                                                                  
                                                                                  
                                                              )),conditionalPanel("input.modelrd == 4 &(input.h1r == 2 |input.h1r == 3|input.h1re == 2|input.h1re == 3)",
                                                                                  em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
                                                              ),conditionalPanel(condition = "input.Moder == 1",
                                                                                                                                                                                      em("\\(\\text{Controlloing for the probability of}\\)"),
                                                                                                                                                                                      
                                                                                                                                                                                      fluidRow(
                                                                                                                                                                                        column(width = 6,sliderInput("powerr", "\\( \\text{True Positive Evidence:} \\)",
                                                                                                                                                                                                                     min =.5, max = .99, value = .8,step = .01,ticks = FALSE)),
                                                                                                                                                                                        column(width = 6,
                                                                                                                                                                                               conditionalPanel("input.h0r == 1|input.h0r == 2",
                                                                                                                                                                                                                sliderInput("alphapr", "\\( \\text{False Positive Evidence:} \\)",
                                                                                                                                                                                                                            min =.001, max = .05, value = .05,step = .001,ticks = FALSE)
                                                                                                                                                                                                                
                                                                                                                                                                                               ))),
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                                                                                                                      
                                                                                  ), conditionalPanel(condition = "input.Moder == 1|input.Moder == 2",sliderInput("br", "\\( \\text{Bound of compelling evidence:} \\)", min=1,max = 20,value = 3,ticks = FALSE)),
             conditionalPanel(condition = "input.Moder == 2",
                              numericInput(inputId ="nr",
                                           label="\\( \\text{Sample Size:} \\)",
                                           value = 50
                                           
                              )),
             
             
             conditionalPanel(condition = "input.Moder == 1|input.Moder == 2",
                             actionButton("runr", label = "\\( \\text{Run} \\)"),
                              conditionalPanel(condition ="input.Moder == 1",
                                               em(span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;")))
                              ),
             
             
             
             conditionalPanel(condition = "input.Moder == 3",
                              fluidRow(
                                column(width = 6,numericInput(
                                  inputId = "rdf",
                                  label="\\( \\text{Degree of freedom:} \\)",
                                  value = 50
                                  
                                )),
                                column(width = 6,sliderInput(
                                  inputId = "rval",
                                  label =  "\\(\\text{Pearson's:} \\)",
                                  min = -1,
                                  ticks = FALSE,
                                  max = 1,
                                  value = 0,  # Default value
                                  step = 0.01    # Step size
                                ))
                                
                              ),
                              actionButton("calr", label = "\\( \\text{Calculate} \\)"),
                              htmlOutput("BFrv"))
             
             
           ),
           mainPanel(
             navset_card_underline(
               nav_panel(id = "mr",em("\\(\\text{Result}\\)"),
                         fluidRow(
                           column(6, plotOutput("prior_r")),
                           column(6,   htmlOutput("resultr"))
                         ),
                         plotOutput("bfrr")),
               nav_panel(id = "pr",em("\\(\\text{Power Curve}\\)"),
                         plotOutput("PCr"))
               
             ))
           
           )
  ),
tabPanel(em("\\(\\text{Regression}\\)"),withMathJax(),
         sidebarLayout(sidebarPanel(
           fluidRow( 
             column(width = 4,prettyRadioButtons(inputId = "ANOREG",
                                                 label = em("\\(\\text{Type of Analysis}\\)"),  # Inline LaTeX for H1
                                                 choices = list(
                                                   "\\(\\text{ANOVA}\\)" = 1,  # Inline LaTeX for delta ≠ 0
                                                   "\\(\\text{Regression}\\)" = 2 # Inline LaTeX for delta < 0
                                                 ),
                                                 inline = T,
                                                 selected = 1  # Default selected option
             )),
             column(width = 6,prettyRadioButtons(
               "Modef",
               "\\(\\text{Select Mode}\\)",
               choices = list("\\(\\text{Sample size determination}\\)" = 1, 
                              "\\(\\text{Fixed N}\\)" = 2, 
                              "\\(\\text{BF calculator}\\)" = 3),
               selected = 1,
               inline = T
             ))),conditionalPanel(condition = "input.ANOREG == 2&(input.Modef ==1|input.Modef ==2)",
                         fluidRow(
                           column(width=6,
                                         sliderInput(
                                           inputId = "pf",
                                           label =  "\\(p\\text{ predictor - reduced model :}\\)",
                                           min = 1,
                                           ticks = FALSE,
                                           max = 100,
                                           value = 1,  # Default value
                                           step = 1    # Step size
                                         )),
                           column(width=6,
                                  sliderInput(
                                    inputId = "kf",
                                    label =  "\\(k\\text{ predictor - full model :}\\)",
                                    min = 1,
                                    ticks = FALSE,
                                    max = 100,
                                    value = 1,  # Default value
                                    step = 1    # Step size
                                  ))
                                  
                                  
                                  
                                  
                                  )         
                                  
                                  
                                  
                                  
                                  ),conditionalPanel(condition = "input.ANOREG == 1&(input.Modef ==1|input.Modef ==2)",
                                  fluidRow( 
                                    column(width = 6,prettyRadioButtons(inputId = "redf",
                                                                        label = em("\\(\\text{Reduced model}\\)"),  # Inline LaTeX for H1
                                                                        choices = list(
                                                                          "\\(\\text{Intercept}\\)" = 1,  # Inline LaTeX for delta ≠ 0
                                                                          "\\(\\text{One-factor }\\)" = 2,
                                                                          "\\(\\text{Two-factor - main effect}\\)" = 3
                                                                        ),
                                                                        inline = T,
                                                                        selected = 1  # Default selected option
                                    )),
                                    column(width = 6,
                                           conditionalPanel(condition = "input.redf == 1",
                                                            
                                                            prettyRadioButtons( "full1", "\\(\\text{Full Model}\\)",
                                                                                choices = list("\\(\\text{One-factor}\\)" = 2,
                                                                                               "\\(\\text{Two-factor - main effect}\\)" = 3,
                                                                                               "\\(\\text{Two-factor - interaction}\\)" = 4), selected = 2,inline = T)),
                                           conditionalPanel(condition = "input.redf == 2",
                                                            
                                                            prettyRadioButtons( "full2", "\\(\\text{Full Model}\\)",
                                                                                choices = list(
                                                                                  "\\(\\text{Two-factor - main effect}\\)" = 3), selected = 3,inline = T)),
                                           conditionalPanel(condition = "input.redf == 3",
                                                            
                                                            prettyRadioButtons( "full3", "\\(\\text{Full Model}\\)",
                                                                                choices = list(
                                                                                  "\\(\\text{Two-factor - interaction}\\)" = 4), selected = 4,inline = T))
                                           
                                           
                                           
                                           
                                    ))),
           
           conditionalPanel(condition = "input.ANOREG == 1&(input.Modef ==1|input.Modef ==2)",
           fluidRow(
             column(width=6,
                    sliderInput(
                      inputId = "f1",
                      label =  "\\(\\text{Factor 1 level:}\\)",
                      min = 2,
                      ticks = FALSE,
                      max = 10,
                      value = 2,  # Default value
                      step = 1    # Step size
                    )),
             
             
             column(width=6,
                    sliderInput(
                      inputId = "f2",
                      label =  "\\(\\text{Factor 2 level:}\\)",
                      min = 2,
                      ticks = FALSE,
                      max = 10,
                      value = 2,  # Default value
                      step = 1    # Step size
                    ))))
             
             
             
           ,
           
           fluidRow(
             column(width = 6, prettyRadioButtons(
               inputId = "h0f",
               label = em("\\(\\mathcal{H}_0:\\)"),  # Inline LaTeX for H0
               choices = list(
                 "\\(\\lambda^2 = 0 \\)" = 1,  
                 "\\(\\lambda^2 \\in \\{0, \\epsilon\\}\\)" = 2   
               ),
               inline = T,
               selected = 1  # Default selected option
             ))
             ,
             
             column(width = 6,
                    conditionalPanel(
                      condition = "input.h0f == 1",  
                      prettyRadioButtons(
                        inputId = "h1f",
                        label = em("\\(\\mathcal{H}_1:\\)"),  
                        choices = list( 
                          "\\(\\lambda^2 > 0 \\)" = 1   # Inline LaTeX for delta < 0
                        ),
                        inline = T,
                        selected = 1  # Default selected option
                      )
                    ),
                    
                    conditionalPanel(
                      condition = "input.h0f == 2",  # Show when delta subset {-ε, ε} is selected
                      prettyRadioButtons(
                        inputId = "h1fe",
                        label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                        choices = list( # Inline LaTeX for delta > ε
                          "\\(\\lambda^2 > \\epsilon\\)" = 1    # Inline LaTeX for delta < ε
                        ),
                        inline = T,
                        selected = 1  # Default selected option
                      ),sliderInput(
  inputId = "epsilinff",
  label =  "\\(\\epsilon :\\)",
  min = 0.01,
  ticks = FALSE,
  max = .25,
  value = .1,  # Default value
  step = 0.01    # Step size
)
                    )
             )
           )
           
           ,prettyRadioButtons(
             inputId = "modelf",
             label = em("\\(\\text{ Analysis Prior Distribution}\\)"),  # Inline LaTeX for H1
             choices = list(
               "\\( \\text{Effect size prior} \\)" = 1,  # Inline LaTeX for delta ≠ 0
               "\\( \\text{Moment prior (must df ≥ 3)} \\)" = 2 
             ),
             inline = T,
             selected = 1  # Default selected option
           ),
           
 fluidRow(
             column(width=4,
                    
                    conditionalPanel(condition = "input.modelf == 1",sliderInput(
                                        inputId = "rf",
                                        label =  "\\(\\text{r scale:}\\)",
                                        min = 0,
                                        ticks = FALSE,
                                        max = 3,
                                        value = 1,  # Default value
                                        step = 0.01    # Step size
                                      ))
                    
                    
                    )
                    ,
             column(width=4,
                    
                    sliderInput(
                      inputId = "fsdf",
                      label =  "\\(\\mathcal{f}^2 :\\)",
                      min = .01,
                      ticks = FALSE,
                      max = .5,
                      value = .1,  # Default value
                      step = 0.01    # Step size
                    )
             ),
             column(width=4,sliderInput(
                                        inputId = "dff",
                                        label =  "\\(\\text{df :}\\)",
                                        min = 1,
                                        ticks = FALSE,
                                        max = 100,
                                        value = 1,  # Default value
                                        step = 1    # Step size
                                      ))),
 conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",
                  prettyRadioButtons(
                    "priorf",
                    "\\( \\text{Design prior is the same as analysis prior: } \\)",
                    choices = list("\\( \\text{Yes} \\)" = 1, 
                                   "\\( \\text{No} \\)"  = 2),
                    selected = 1,
                    inline = T
                  ),
                  
                  
conditionalPanel(condition =    "input.priorf == 2",
                  
                  prettyRadioButtons(
                    inputId = "modelfd",
                    label = em("\\(\\text{ Analysis Prior Distribution}\\)"),  # Inline LaTeX for H1
                    choices = list(
                      "\\( \\text{Effect size prior} \\)" = 1,  # Inline LaTeX for delta ≠ 0
                      "\\( \\text{Moment prior (must df ≥ 3)} \\)" = 2 ,
                      "\\( \\text{Point} \\)" = 3
                    ),inline = T,selected = 1  ),
                  
fluidRow(
             column(width=4,
                    
conditionalPanel(condition = "input.modelfd == 1",
                 sliderInput(
                                        inputId = "rfd",
                                        label =  "\\(\\text{r scale:}\\)",
                                        min = 0,
                                        ticks = FALSE,
                                        max = 3,
                                        value = 1,  # Default value
                                        step = 0.01    # Step size
                                      )),
conditionalPanel(condition = "input.modelfd == 3",sliderInput(
  inputId = "lfd",
  label =  "\\(\\lambda^2:\\)",
  min = 0,
  ticks = FALSE,
  max = .5,
  value = .1,  # Default value
  step = 0.01    # Step size
))




)
                    ,
             column(width=4,
                    conditionalPanel(condition = "input.modelfd == 1|input.modelfd == 2",
                    sliderInput(
                      inputId = "fsdfd",
                      label =  "\\(\\mathcal{f}^2 :\\)",
                      min = .01,
                      ticks = FALSE,
                      max = .5,
                      value = .1,  # Default value
                      step = 0.01    # Step size
                    ))
             ),
             column(width=4,
                    conditionalPanel(condition = "input.modelfd == 1|input.modelfd == 2",
                    sliderInput(inputId = "dffd",
                                        label =  "\\(\\text{df :}\\)",
                                        min = 1,
                                        ticks = FALSE,
                                        max = 100,
                                        value = 1,  # Default value
                                        step = 1    # Step size
                                      ))
                    )
)                  
                  
                  
                  
                  
                  
                  
                  
                  
                  ),
conditionalPanel(condition = "input.Modef == 1",
fluidRow(
  column(width = 6,sliderInput("powerf", "\\( \\text{True Positive Evidence:} \\)",
                               min =.5, max = .99, value = .8,step = .01,ticks = FALSE)),
  column(width = 6,
         sliderInput("alphaf", "\\( \\text{False Positive Evidence:} \\)",
                                      min =.001, max = .05, value = .05,step = .001,ticks = FALSE)
                          
         ))),
conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",sliderInput("bff", "\\( \\text{Bound of compelling evidence:} \\)", min=1,max = 20,value = 3,ticks = FALSE)),

conditionalPanel(condition = "input.Modef == 2",numericInput(inputId ="nf",
             label="\\( \\text{Sample Size } N: \\)",
             value = 50)),
conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",
                 actionButton("runf", label = "\\( \\text{Run} \\)"),
                 conditionalPanel(condition ="input.Modef == 1",
                                  em(span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;")))
                 
                 
                 ),



),conditionalPanel(condition = "input.Modef == 3",
                   fluidRow(
                     column(width = 4,numericInput(
                       inputId = "df1f",
                       label="\\( \\mathcal{df}_1: \\)",
                       value = 1
                       
                     )),
                     column(width = 4,numericInput(
                       inputId = "df2f",
                       label="\\( \\mathcal{df}_2: \\)",
                       value = 30
                       
                     )),
                     column(width = 4,numericInput(
                       inputId = "fval",
                       label="\\( f\\text{-value:} \\)",
                       value = 1
                       
                     ))
                     
                   ),
                   actionButton("calf", label = "\\( \\text{Calculate} \\)"),
                   htmlOutput("BFcalf")),

em("\\( \\text{Recommended hyperparameters:} \\)"),
htmlOutput("prior_suggest")

           
           
           
           
           
           
           
           
           
           
           
           
         ), mainPanel(
           navset_card_underline(
             nav_panel(id = "mf",em("\\(\\text{Result}\\)"),
                       fluidRow(
                         column(6, plotOutput("priorff")),
                         column(6,   htmlOutput("resultf"))
                       ),
                       plotOutput("bff")),
             nav_panel(id = "pf2",em("\\(\\text{Power Curve}\\)"),
                       plotOutput("PCf"))
             
           ))
         
         
         )),


navbarMenu(
  "\\(\\text{Proportion}\\)",
  tabPanel("\\(\\text{One proportion - binomial}\\)",withMathJax(),
           sidebarLayout(
    sidebarPanel(prettyRadioButtons(
      "Modebin",
      "\\(\\text{Select Mode}\\)",
      choices = list("\\(\\text{Sample size determination}\\)" = 1, 
                     "\\(\\text{Fixed N}\\)" = 2, 
                     "\\(\\text{BF calculator}\\)" = 3),
      selected = 1,
      inline = T
    ),
    fluidRow(
      column(width = 6, prettyRadioButtons(
        inputId = "h0bin",
        label = em("\\(\\mathcal{H}_0:\\)"),  # Inline LaTeX for H0
        choices = list(
          "\\(p = p_0\\)"  = 1,  # Inline LaTeX for delta = 0
          "\\(p \\in (p_0 - \\epsilon, p_0 + \\epsilon)\\)" = 2   # Inline LaTeX for delta = 1
        ),
        inline = T,
        selected = 1  # Default selected option
      ))
      ,
      
      column(width = 6,
             conditionalPanel(
               condition = "input.h0bin == 1",  # Show when delta = 0 is selected
               prettyRadioButtons(
                 inputId = "h1bin",
                 label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                 choices = list(
                   "\\(p ≠ p_0\\)" = 1,  # Inline LaTeX for delta ≠ 0
                   "\\(p > p_0\\)" = 2,   # Inline LaTeX for delta > 0
                   "\\(p < p_0\\)" = 3    # Inline LaTeX for delta < 0
                 ),
                 inline = T,
                 selected = 1  # Default selected option
               )
             ),
             
             conditionalPanel(
               condition = "input.h0bin == 2",  # Show when delta subset {-ε, ε} is selected
               prettyRadioButtons(
                 inputId = "h1bine",
                 label = em("\\(\\mathcal{H}_1:\\)"),  # Inline LaTeX for H1
                 choices = list(
                   "\\(p \\not\\in \\{p_0 -\\epsilon, p_0 +\\epsilon\\}\\)" = 1,  # Inline LaTeX for delta ≠ {-ε, ε}
                   "\\(p > p_0 + \\epsilon\\)" = 2,   # Inline LaTeX for delta > ε
                   "\\(p < p_0 - \\epsilon\\)" = 3    # Inline LaTeX for delta < ε
                 ),
                 inline = T,
                 selected = 1  # Default selected option
               )
             )
      )
    ),
    fluidRow(
      column(width=4,sliderInput(
        inputId = "h0prop",
        label =  "\\(p_0\\)",
        min = .01,
        ticks = FALSE,
        max = .99,
        value = .5,  # Default value
        step = 0.01    # Step size
      )),
      column(width=4,
             conditionalPanel( condition =  "(input.h1bine == 1 || input.h1bine == 3) && input.h0bin == 2",
                               sliderInput(
                                 inputId = "lbbine",
                                 label =  "\\( -\\epsilon \\)",
                                 min = -.5,
                                 ticks = FALSE,
                                 max = -0.01,
                                 value = -0.2,  # Default value
                                 step = 0.01    # Step size
                               ),htmlOutput("bin_lower")
             )),
      column(width=4,
             conditionalPanel( condition =  "(input.h1bine == 1 || input.h1bine == 2) && input.h0bin == 2",
                               sliderInput(
                                 inputId = "ubbine",
                                 label = "\\(\\epsilon \\)",
                                 min = .01,
                                 max = .5,
                                 value = 0.2,  # Default value
                                 step = 0.01,   # Step size
                                 ticks = FALSE),htmlOutput("bin_upper")
             ))),
    
    prettyRadioButtons(
      inputId = "modelbin",
      label = em("\\(\\text{ Analysis Prior Distribution}\\)"),  # Inline LaTeX for H1
      choices = list(
        "\\( \\text{Beta} \\)" = 1,  
        "\\( \\text{Moment} \\)" = 2
      ),
      inline = T,
      selected = 1  # Default selected option
    ),
    fluidRow(
      column(width=4,
             conditionalPanel( condition = "input.modelbin == 1",
                               sliderInput(
                                 inputId = "alphabin",
                                 label =  "\\(\\alpha \\)",
                                 min = 0.01,
                                 ticks = FALSE,
                                 max = 100,
                                 value = 1,  # Default value
                                 step = 1    # Step size
                               )
             )
             ,conditionalPanel( condition = "input.modelbin == 2",
                                 sliderInput(
                                   inputId = "sbin",
                                   label =  "\\(Scale \\)",
                                   min = 0.01,
                                   ticks = FALSE,
                                   max = 3,
                                   value = 1,  # Default value
                                   step = .01    # Step size
                                 )
             )
             
             
      ),
      column(width=4,
             conditionalPanel( condition = "input.modelbin == 1",
                               sliderInput(
                                 inputId = "betabin",
                                 label =  "\\(\\beta \\)",
                                 min = 0.01,
                                 ticks = FALSE,
                                 max = 100,
                                 value = 1,  # Default value
                                 step = 0.01    # Step size
                               )
             ))),
    
    conditionalPanel(condition = "input.Modebin == 1|input.Modebin == 2",
                     prettyRadioButtons(
                       "priorbin",
                       "\\( \\text{Design prior is the same as analysis prior: } \\)",
                       choices = list("\\( \\text{Yes} \\)" = 1, 
                                      "\\( \\text{No} \\)"  = 2),
                       selected = 1,
                       inline = T
                     ),
      conditionalPanel( condition = "input.priorbin == 2",
      prettyRadioButtons(inputId = "modelbind",
                         label = em("\\( \\text{Design prior distribution} \\)"),  # Inline LaTeX for H1
                         choices = list("\\( \\text{Beta} \\)" = 1,  
                                        "\\( \\text{Moment} \\)" = 2,
                                        "\\( \\text{Point} \\)" = 3
                                           ),
                                           inline = T,
                                           selected = 1  # Default selected option
                                         ),
      conditionalPanel(condition = "input.modelbind == 3",
                       sliderInput( inputId = "h0bind",
                                    label =  "\\(p_1\\)",
                                    min = .01,
                                    ticks = FALSE,
                                    max = .99,
                                    value = .5,  # Default value
                                    step = 0.01)),
fluidRow(
  column(width=4,
         conditionalPanel( condition = "input.modelbind == 1",
                           sliderInput(
                             inputId = "alphabind",
                             label =  "\\(\\alpha \\)",
                             min = 0.01,
                             ticks = FALSE,
                             max = 100,
                             value = 1,  # Default value
                             step = 1    # Step size
                           )),
                                                  
conditionalPanel( condition = "input.modelbind == 2",
                  sliderInput(
                    inputId = "sbind",
                    label =  "\\(Scale \\)",
                    min = 0.01,
                    ticks = FALSE,
                    max = 3,
                    value = 1,  # Default value
                    step = .01    # Step size
                  )
                  )),
column(width=4,
       conditionalPanel(condition = "input.modelbind == 1",
                        sliderInput(
                          inputId = "betabind",
                          label =  "\\(\\beta \\)",
                          min = 0.01,
                          ticks = FALSE,
                          max = 100,
                          value = 1,  # Default value
                          step = 0.01    # Step size
                        ))))
)),

conditionalPanel("input.modelbind == 3 &(input.h1bin == 2 |input.h1bin == 3|input.h1bine == 2|input.h1bine == 3)",
                 em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
),
conditionalPanel(condition = "input.Modebin == 1",
                 em("\\(\\text{Controlloing for the probability of}\\)"),
fluidRow(
  column(width = 6,
         sliderInput("powerbin", "\\( \\text{True Positive Evidence:} \\)",
       min =.5, max = .99, value = .8,step = .01,ticks = FALSE)),
       column(width = 6,
         sliderInput("FP_bin", "\\( \\text{False Positive Evidence:} \\)",
              min =.001, max = .05, value = .05,step = .001,ticks = FALSE))),
                                                                                                                                                                                      
                                                                   ),
conditionalPanel(condition = "input.Modebin == 1|input.Modebin == 2",
                 sliderInput("bbin", "\\( \\text{Bound of compelling evidence:} \\)", min=1,max = 20,value = 3,ticks = FALSE)),

conditionalPanel(condition = "input.Modebin == 2|input.Modebin == 3",
fluidRow(
  column(width = 6,
         numericInput(inputId ="nbin",
                      label="\\( \\text{Sample Size:} \\)",
                      value = 50
                      
         )

         ),
  
  column(width = 6, 
         conditionalPanel(condition = "input.Modebin == 3",
         numericInput(inputId ="xbin",
                                 label="\\( \\text{Number of Success:} \\)",
                                 value = 25
                                 
  )))
  
) ),
conditionalPanel(condition = "input.Modebin == 1|input.Modebin == 2",
                 actionButton("runbin", label = "\\( \\text{Run} \\)"),
                 conditionalPanel(condition ="input.Modebin == 1",
                                  em(span("\\(\\text{Note: Error when the required N > 50,000} \\)", style = "color: red;")))
                 
                 
                 ),
conditionalPanel(condition = "input.Modebin == 3",
                 actionButton("calbin", label = "\\( \\text{Calculate} \\)"),
                 htmlOutput("BFbin"))

    
    ),
                         
                         
                         
    mainPanel(
      navset_card_underline(
        nav_panel(id = "mbin",em("\\(\\text{Result}\\)"),
                  fluidRow(
                    column(6, plotOutput("prior_bin")),
                    column(6,   htmlOutput("resultbin"))
                  ),
                  plotOutput("bfbin")),
        nav_panel(id = "pbin",em("\\(\\text{Power Curve}\\)"),
                  plotOutput("PCbin"))
        
      )
      
      
      
      
      
    ))
           )
  
  )
)

# Server logic
server <- function(input, output, session) {
  
source("Server/Server_t1.r",local = T)
  
source("Server/Server_t2.r",local = T)
  
source("Server/Server_r.r",local = T)
  
source("Server/Server_f.r",local = T)
  
source("Server/Server_bin.r",local = T)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

