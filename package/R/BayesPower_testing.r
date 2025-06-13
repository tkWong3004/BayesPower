

ui <- navbarPage(id = "id",
                 "\\(\\text{BayesPower}_{1.0}\\)",
  navbarMenu(
    "\\(\\text{Standardized Mean Difference}\\)",

    tabPanel(
      "\\(\\text{One-sample/paired t-test}\\)",

      # Custom font for MathJax
      tags$style(HTML("
    body {
  font-family: sans-serif;
}
  ")),

      withMathJax(),

      sidebarLayout(
        sidebarPanel(
          # Mode selection
          prettyRadioButtons(
            "Modet1",
            "\\(\\text{Select Mode}\\)",
            choices = list(
              "\\(\\text{Sample size determination}\\)" = 1,
              "\\(\\text{Fixed N}\\)" = 2,
              "\\(\\text{BF calculator}\\)" = 3
            ),
            selected = 1,
            inline = TRUE
          ),

          # Hypotheses selection
          fluidRow(
            column(6,
                   prettyRadioButtons(
                     inputId = "h0t1",
                     label = em("\\(\\mathcal{H}_0:\\)"),
                     choices = list(
                       "\\(\\delta = 0\\)" = 1,
                       "\\(\\delta \\in \\{-\\epsilon, \\epsilon\\}\\)" = 2
                     ),
                     inline = TRUE,
                     selected = 1
                   )
            ),
            column(6,
                   conditionalPanel("input.h0t1 == 1",
                                    prettyRadioButtons(
                                      inputId = "h1t1",
                                      label = em("\\(\\mathcal{H}_1:\\)"),
                                      choices = list(
                                        "\\(\\delta ≠ 0\\)" = 1,
                                        "\\(\\delta > 0\\)" = 2,
                                        "\\(\\delta < 0\\)" = 3
                                      ),
                                      inline = TRUE,
                                      selected = 1
                                    )
                   ),
                   conditionalPanel("input.h0t1 == 2",
                                    prettyRadioButtons(
                                      inputId = "h1t1e",
                                      label = em("\\(\\mathcal{H}_1:\\)"),
                                      choices = list(
                                        "\\(\\delta \\not\\in \\{-\\epsilon, \\epsilon\\}\\)" = 1,
                                        "\\(\\delta > \\epsilon\\)" = 2,
                                        "\\(\\delta < \\epsilon\\)" = 3
                                      ),
                                      inline = TRUE,
                                      selected = 1
                                    )
                   )
            )
          ),

          # ε inputs
          fluidRow(
            column(6,
                   conditionalPanel("input.h1t1e == 2 && input.h0t1 == 2",
                                    em("\\( -\\epsilon = 0 \\)")
                   ),
                   conditionalPanel("(input.h1t1e == 1 || input.h1t1e == 3) && input.h0t1 == 2",
                                    sliderInput("lbt1e", "\\( -\\epsilon \\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE)
                   )
            ),
            column(6,
                   conditionalPanel("input.h1t1e == 3 && input.h0t1 == 2",
                                    em("\\( \\epsilon = 0 \\)")
                   ),
                   conditionalPanel("(input.h1t1e == 1 || input.h1t1e == 2) && input.h0t1 == 2",
                                    sliderInput("ubt1e", "\\( \\epsilon \\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE)
                   )
            )
          ),

          # Analysis prior
          prettyRadioButtons(
            inputId = "modelt1",
            label = em("\\(\\text{Analysis Prior Distribution}\\)"),
            choices = list(
              "\\(\\text{Scaled t}\\)" = 1,
              "\\(\\text{Normal}\\)" = 2,
              "\\(\\text{Moment}\\)" = 3
            ),
            inline = TRUE,
            selected = 1
          ),

          # Location/scale/df
          fluidRow(
            column(4,
                   conditionalPanel("input.h0t1 == 1",
                                    sliderInput("lt1", "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
                   ),
                   conditionalPanel("input.h0t1 == 2", em("\\(\\text{Location = 0}\\)"))
            ),
            column(4,
                   sliderInput("st1", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 0.707, step = 0.001, ticks = FALSE)
            ),
            column(4,
                   conditionalPanel("input.modelt1 == 1",
                                    sliderInput("dft1", "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                   )
            )
          ),

          # Design prior if different
          conditionalPanel("input.Modet1 == 1 || input.Modet1 == 2",
                           prettyRadioButtons(
                             "prior",
                             "\\(\\text{Design prior is the same as analysis prior:}\\)",
                             choices = list("\\(\\text{Yes}\\)" = 1, "\\(\\text{No}\\)" = 2),
                             selected = 1,
                             inline = TRUE
                           ),
                           conditionalPanel("input.prior == 2",
                                            prettyRadioButtons(
                                              inputId = "modelt1d",
                                              label = em("\\(\\text{Design Prior Distribution}\\)"),
                                              choices = list(
                                                "\\(\\text{Scaled t}\\)" = 1,
                                                "\\(\\text{Normal}\\)" = 2,
                                                "\\(\\text{Moment}\\)" = 3,
                                                "\\(\\text{Point}\\)" = 4
                                              ),
                                              selected = 1,
                                              inline = TRUE
                                            ),
                                            fluidRow(
                                              column(4,
                                                     conditionalPanel("input.h0t1 == 1 || input.modelt1d == 4",
                                                                      sliderInput("lt1d", "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
                                                     ),
                                                     conditionalPanel("input.h0t1 == 2 && input.modelt1d != 4", em("\\(\\text{Location = 0}\\)"))
                                              ),
                                              column(4,
                                                     conditionalPanel("input.modelt1d != 4",
                                                                      sliderInput("st1d", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 0.707, step = 0.001, ticks = FALSE)
                                                     )
                                              ),
                                              column(4,
                                                     conditionalPanel("input.modelt1d == 1",
                                                                      sliderInput("dft1d", "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                                                     )
                                              )
                                            )
                           )
          ),

          # Note for point prior
          conditionalPanel("input.modelt1d == 4 && (input.h1t1 == 2 || input.h1t1 == 3 || input.h1t1e == 2 || input.h1t1e == 3)",
                           em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
          ),

          # Power / error control
          conditionalPanel("input.Modet1 == 1",
                           em("\\(\\text{Controlling for the probability of}\\)"),
                           fluidRow(
                             column(6,
                                    sliderInput("powert1", "\\(\\text{True Positive Evidence:}\\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)
                             ),
                             column(6,
                                    sliderInput("alphat1", "\\(\\text{False Positive Evidence:}\\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE)
                             )
                           )
          ),

          # Bound
          conditionalPanel("input.Modet1 == 1 || input.Modet1 == 2",
                           sliderInput("bt1", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)
          ),

          # Sample size input
          conditionalPanel("input.Modet1 == 2",
                           numericInput("nt1", "\\(\\text{Sample Size:}\\)", value = 50)
          ),

          # Run button + error message
          conditionalPanel("input.Modet1 == 1 || input.Modet1 == 2",
                           actionButton("runt1", label = "\\(\\text{Run}\\)"),
                           conditionalPanel("input.Modet1 == 1",
                                            em(span("\\(\\text{Note: Error when the required N > 10,000}\\)", style = "color: red;"))
                           )
          ),

          # BF calculator mode
          conditionalPanel("input.Modet1 == 3",
                           fluidRow(
                             column(6, numericInput("t1df", "\\(\\text{Degree of freedom:}\\)", value = 50)),
                             column(6, numericInput("t1tval", "\\(\\text{t-value:}\\)", value = 2))
                           ),
                           actionButton("cal1", label = "\\(\\text{Calculate}\\)"),
                           htmlOutput("BFt1")

          ),

          conditionalPanel(
            condition = "input.Modet1 == 1",
            checkboxGroupInput(
              "o_plot_t1",
              label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
              choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
              selected = NULL
            ),
            downloadButton("export_t1", "Download result as PDF")
          )
        ),

        # Main panel with result tabs
        mainPanel( fluidRow(
          column(6, plotOutput("priort1")),
          column(6, htmlOutput("resultt1"))),
          uiOutput("Optional_Plots_t1"))
      )
    )
    ,tabPanel(
      "\\(\\text{Independent samples t-test}\\)",
      withMathJax(),
      sidebarLayout(
        sidebarPanel(
          # Mode selection
          prettyRadioButtons(
            inputId = "Modet2",
            label = "\\(\\text{Select Mode}\\)",
            choices = list(
              "\\(\\text{Sample size determination}\\)" = 1,
              "\\(\\text{Fixed N}\\)" = 2,
              "\\(\\text{BF calculator}\\)" = 3
            ),
            selected = 1,
            inline = TRUE
          ),

          # Hypotheses
          fluidRow(
            column(
              width = 6,
              prettyRadioButtons(
                inputId = "h0t2",
                label = em("\\(\\mathcal{H}_0:\\)"),
                choices = list(
                  "\\(\\delta = 0\\)" = 1,
                  "\\(\\delta \\in \\{-\\epsilon, \\epsilon\\}\\)" = 2
                ),
                selected = 1,
                inline = TRUE
              )
            ),
            column(
              width = 6,
              conditionalPanel(
                condition = "input.h0t2 == 1",
                prettyRadioButtons(
                  inputId = "h1t2",
                  label = em("\\(\\mathcal{H}_1:\\)"),
                  choices = list(
                    "\\(\\delta ≠ 0\\)" = 1,
                    "\\(\\delta > 0\\)" = 2,
                    "\\(\\delta < 0\\)" = 3
                  ),
                  selected = 1,
                  inline = TRUE
                )
              ),
              conditionalPanel(
                condition = "input.h0t2 == 2",
                prettyRadioButtons(
                  inputId = "h1t2e",
                  label = em("\\(\\mathcal{H}_1:\\)"),
                  choices = list(
                    "\\(\\delta \\not\\in \\{-\\epsilon, \\epsilon\\}\\)" = 1,
                    "\\(\\delta > \\epsilon\\)" = 2,
                    "\\(\\delta < \\epsilon\\)" = 3
                  ),
                  selected = 1,
                  inline = TRUE
                )
              )
            )
          ),

          # Epsilon sliders
          fluidRow(
            column(
              width = 6,
              conditionalPanel("input.h1t2e == 2 && input.h0t2 == 2", em("\\(-\\epsilon = 0\\)")),
              conditionalPanel(
                condition = "(input.h1t2e == 1 || input.h1t2e == 3) && input.h0t2 == 2",
                sliderInput("lbt2e", label = "\\(-\\epsilon\\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE)
              )
            ),
            column(
              width = 6,
              conditionalPanel("input.h1t2e == 3 && input.h0t2 == 2", em("\\(\\epsilon = 0\\)")),
              conditionalPanel(
                condition = "(input.h1t2e == 1 || input.h1t2e == 2) && input.h0t2 == 2",
                sliderInput("ubt2e", label = "\\(\\epsilon\\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE)
              )
            )
          ),

          # Analysis prior
          prettyRadioButtons(
            inputId = "modelt2",
            label = em("\\(\\text{Analysis Prior Distribution}\\)"),
            choices = list(
              "\\(\\text{Scaled t}\\)" = 1,
              "\\(\\text{Normal}\\)" = 2,
              "\\(\\text{Moment}\\)" = 3
            ),
            selected = 1,
            inline = TRUE
          ),

          fluidRow(
            column(
              width = 4,
              conditionalPanel(
                "input.h0t2 == 1",
                sliderInput("lt2", label = "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
              ),
              conditionalPanel("input.h0t2 == 2", em("\\(\\text{Location = 0}\\)"))
            ),
            column(
              width = 4,
              sliderInput("st2", label = "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE)
            ),
            column(
              width = 4,
              conditionalPanel(
                "input.modelt2 == 1",
                sliderInput("dft2", label = "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
              )
            )
          ),

          # Design prior toggle
          conditionalPanel(
            "input.Modet2 == 1 || input.Modet2 == 2",
            prettyRadioButtons(
              "priort2",
              "\\(\\text{Design prior is the same as analysis prior:}\\)",
              choices = list("\\(\\text{Yes}\\)" = 1, "\\(\\text{No}\\)" = 2),
              selected = 1,
              inline = TRUE
            ),
            conditionalPanel(
              "input.priort2 == 2",
              prettyRadioButtons(
                inputId = "modelt2d",
                label = em("\\(\\text{Design prior distribution}\\)"),
                choices = list(
                  "\\(\\text{Scaled t}\\)" = 1,
                  "\\(\\text{Normal}\\)" = 2,
                  "\\(\\text{Moment}\\)" = 3,
                  "\\(\\text{Point}\\)" = 4
                ),
                selected = 1,
                inline = TRUE
              ),
              fluidRow(
                column(
                  width = 4,
                  conditionalPanel(
                    "input.h0t2 == 1 || input.modelt2d == 4",
                    sliderInput("lt2d", label = "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
                  ),
                  conditionalPanel("input.h0t2 == 2 && input.modelt2d != 4", em("\\(\\text{Location = 0}\\)"))
                ),
                column(
                  width = 4,
                  conditionalPanel(
                    "input.modelt2d == 1 || input.modelt2d == 2 || input.modelt2d == 3",
                    sliderInput("st2d", label = "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE)
                  )
                ),
                column(
                  width = 4,
                  conditionalPanel(
                    "input.modelt2d == 1",
                    sliderInput("dft2d", label = "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                  )
                )
              )
            )
          ),

          # Hint for point prior direction
          conditionalPanel(
            "input.modelt2d == 4 && (input.h1t2 == 2 || input.h1t2 == 3 || input.h1t2e == 2 || input.h1t2e == 3)",
            em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
          ),

          # Controls for sample size determination
          conditionalPanel("input.Modet2 == 1",
                           em("\\(\\text{Controlling for the probability of}\\)"),
                           fluidRow(
                             column(6, sliderInput("powert2", "\\(\\text{True Positive Evidence:}\\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)),
                             column(6, sliderInput("alphat2", "\\(\\text{False Positive Evidence:}\\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE))
                           )
          ),

          # Shared parameters
          conditionalPanel("input.Modet2 == 1 || input.Modet2 == 2", sliderInput("bt2", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)),
          conditionalPanel("input.Modet2 == 1 || input.Modet2 == 3", sliderInput("rt2", "\\(N_2/N_1\\)", min = 1, max = 10, value = 1, ticks = FALSE)),

          # Fixed sample sizes
          conditionalPanel("input.Modet2 == 2",
                           fluidRow(
                             column(4, numericInput("n1t2", "\\(N_1:\\)", value = 50)),
                             column(4, numericInput("n2t2", "\\(N_2:\\)", value = 50))
                           )
          ),

          # Action buttons
          conditionalPanel("input.Modet2 == 1 || input.Modet2 == 2",
                           actionButton("runt2", label = "\\(\\text{Run}\\)"),
                           conditionalPanel("input.Modet2 == 1", em(span("\\(\\text{Note: Error when required } N > 10,000\\)", style = "color: red;")))
          ),
          conditionalPanel("input.Modet2 == 3",
                           fluidRow(
                             column(6, numericInput("t2df", "\\(\\text{Degrees of freedom:}\\)", value = 50)),
                             column(6, numericInput("t2tval", "\\(\\text{t-value:}\\)", value = 2))
                           ),
                           actionButton("cal1", label = "\\(\\text{Calculate}\\)"),
                           htmlOutput("BFt2")
          ),
          conditionalPanel(
            condition = "input.Modet2 == 1",
            checkboxGroupInput(
              "o_plot_t2",
              label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
              choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
              selected = NULL
            ),
            downloadButton("export_t2", "Download result as PDF")
          )





        ),

        mainPanel(fluidRow(
          column(6, plotOutput("priort2")),
          column(6, htmlOutput("resultt2"))),
          uiOutput("Optional_Plots_t2")





          )
      )
    )
  )
,
tabPanel("\\(\\text{Correlation}\\)", withMathJax(),
         sidebarLayout(
           sidebarPanel(
             prettyRadioButtons(
               inputId = "Moder",
               label = "\\(\\text{Select Mode}\\)",
               choices = list(
                 "\\(\\text{Sample size determination}\\)" = 1,
                 "\\(\\text{Fixed N}\\)" = 2,
                 "\\(\\text{BF calculator}\\)" = 3
               ),
               selected = 1,
               inline = TRUE
             ),

             # H0 and H1 Specification
             fluidRow(
               column(6,
                      prettyRadioButtons(
                        inputId = "h0r",
                        label = em("\\(\\mathcal{H}_0:\\)"),
                        choices = list(
                          "\\(\\rho = \\rho_0\\)" = 1,
                          "\\(\\rho \\in (\\rho_0 - \\epsilon, \\rho_0 + \\epsilon)\\)" = 2
                        ),
                        selected = 1,
                        inline = TRUE
                      )
               ),
               column(6,
                      conditionalPanel(
                        condition = "input.h0r == 1",
                        prettyRadioButtons(
                          inputId = "h1r",
                          label = em("\\(\\mathcal{H}_1:\\)"),
                          choices = list(
                            "\\(\\rho ≠ \\rho_0\\)" = 1,
                            "\\(\\rho > \\rho_0\\)" = 2,
                            "\\(\\rho < \\rho_0\\)" = 3
                          ),
                          selected = 1,
                          inline = TRUE
                        )
                      ),
                      conditionalPanel(
                        condition = "input.h0r == 2",
                        prettyRadioButtons(
                          inputId = "h1re",
                          label = em("\\(\\mathcal{H}_1:\\)"),
                          choices = list(
                            "\\(\\rho \\not\\in \\{\\rho_0 -\\epsilon, \\rho_0 +\\epsilon\\}\\)" = 1,
                            "\\(\\rho > \\rho_0 + \\epsilon\\)" = 2,
                            "\\(\\rho < \\rho_0 - \\epsilon\\)" = 3
                          ),
                          selected = 1,
                          inline = TRUE
                        )
                      )
               )
             ),

             # rho_0 and epsilon inputs
             fluidRow(
               column(4,
                      sliderInput("h0pho", "\\(\\rho_0\\)", min = -.99, max = .99, value = 0, step = 0.01, ticks = FALSE)
               ),
               column(4,
                      conditionalPanel(
                        condition = "(input.h1re == 1 || input.h1re == 3) && input.h0r == 2",
                        sliderInput("lbre", "\\( -\\epsilon \\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE),
                        htmlOutput("r_lower")
                      )
               ),
               column(4,
                      conditionalPanel(
                        condition = "(input.h1re == 1 || input.h1re == 2) && input.h0r == 2",
                        sliderInput("ubre", "\\(\\epsilon \\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE),
                        htmlOutput("r_upper")
                      )
               )
             ),

             # Analysis prior
             prettyRadioButtons(
               inputId = "modelr",
               label = em("\\(\\text{ Analysis Prior Distribution}\\)"),
               choices = list(
                 "\\( \\text{Default Stretched Beta} \\)" = 1,
                 "\\( \\text{Stretched Beta} \\)" = 2,
                 "\\( \\text{Moment} \\)" = 3
               ),
               selected = 1,
               inline = TRUE
             ),

             fluidRow(
               column(4,
                      conditionalPanel("input.modelr == 1",
                                       sliderInput("kr", "\\(k \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                      ),
                      conditionalPanel("input.modelr == 3",
                                       sliderInput("sr", "\\(Scale \\)", min = 0.01, max = 1, value = 0.01, step = 0.01, ticks = FALSE)
                      )
               ),
               column(4,
                      conditionalPanel("input.modelr == 2",
                                       sliderInput("ralpha", "\\(\\alpha \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                      )
               ),
               column(4,
                      conditionalPanel("input.modelr == 2",
                                       sliderInput("rbeta", "\\(\\beta \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                      )
               )
             ),

             # Design prior
             conditionalPanel("input.Moder == 1 | input.Moder == 2",
                              prettyRadioButtons(
                                inputId = "priorr",
                                label = "\\( \\text{Design prior is the same as analysis prior: } \\)",
                                choices = list(
                                  "\\( \\text{Yes} \\)" = 1,
                                  "\\( \\text{No} \\)" = 2
                                ),
                                selected = 1,
                                inline = TRUE
                              ),
                              conditionalPanel("input.priorr == 2",
                                               prettyRadioButtons(
                                                 inputId = "modelrd",
                                                 label = em("\\( \\text{Design prior distribution} \\)"),
                                                 choices = list(
                                                   "\\( \\text{Default Stretched Beta} \\)" = 1,
                                                   "\\( \\text{Stretched Beta} \\)" = 2,
                                                   "\\( \\text{Moment} \\)" = 3,
                                                   "\\( \\text{Point} \\)" = 4
                                                 ),
                                                 selected = 1,
                                                 inline = TRUE
                                               ),
                                               conditionalPanel("input.modelrd == 4",
                                                                sliderInput("h0phod", "\\(\\rho_1\\)", min = -1, max = 1, value = 0.3, step = 0.01, ticks = FALSE)
                                               ),
                                               fluidRow(
                                                 column(4,
                                                        conditionalPanel("input.modelrd == 1",
                                                                         sliderInput("rkd", "\\( k \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                                                        ),
                                                        conditionalPanel("input.modelrd == 3",
                                                                         sliderInput("rsd", "\\( Scale \\)", min = 0.01, max = 1, value = 0.1, step = 0.01, ticks = FALSE)
                                                        )
                                                 ),
                                                 column(4,
                                                        conditionalPanel("input.modelrd == 2",
                                                                         sliderInput("ralphad", "\\(\\alpha \\)", min = 0.01, max = 10, value = 0.1, step = 0.01, ticks = FALSE)
                                                        )
                                                 ),
                                                 column(4,
                                                        conditionalPanel("input.modelrd == 2",
                                                                         sliderInput("rbetad", "\\(\\beta \\)", min = 0.01, max = 10, value = 0.1, step = 0.01, ticks = FALSE)
                                                        )
                                                 )
                                               )
                              ),
                              conditionalPanel("input.modelrd == 4 & (input.h1r == 2 | input.h1r == 3 | input.h1re == 2 | input.h1re == 3)",
                                               em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
                              )
             ),

             # Power planning
             conditionalPanel("input.Moder == 1",
                              em("\\(\\text{Controlling for the probability of}\\)"),
                              fluidRow(
                                column(6,
                                       sliderInput("powerr", "\\( \\text{True Positive Evidence:} \\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)
                                ),
                                column(6,
                                       conditionalPanel("input.h0r == 1 | input.h0r == 2",
                                                        sliderInput("alphapr", "\\( \\text{False Positive Evidence:} \\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE)
                                       )
                                )
                              )
             ),

             # Common for Mode 1 and 2
             conditionalPanel("input.Moder == 1 | input.Moder == 2",
                              sliderInput("br", "\\( \\text{Bound of compelling evidence:} \\)", min = 1, max = 20, value = 3, ticks = FALSE)
             ),

             conditionalPanel("input.Moder == 2",
                              numericInput("nr", "\\( \\text{Sample Size:} \\)", value = 50)
             ),

             conditionalPanel("input.Moder == 1 | input.Moder == 2",
                              actionButton("runr", label = "\\( \\text{Run} \\)"),
                              conditionalPanel("input.Moder == 1",
                                               em(span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;"))
                              )
             ),

             # BF calculator
             conditionalPanel("input.Moder == 3",
                              fluidRow(
                                column(6,
                                       numericInput("rdf", "\\( \\text{Degree of freedom:} \\)", value = 50)
                                ),
                                column(6,
                                       sliderInput("rval", "\\(\\text{Pearson's:} \\)", min = -1, max = 1, value = 0, step = 0.01, ticks = FALSE)
                                )
                              ),
                              actionButton("calr", label = "\\( \\text{Calculate} \\)"),
                              htmlOutput("BFrv")
             ),
             conditionalPanel(
               condition = "input.Moder == 1",
               checkboxGroupInput(
                 "o_plot_r",
                 label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
                 choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
                 selected = NULL
               ),
               downloadButton("export_r", "Download result as PDF")
             )
           ),

           # Main panel with output tabs
           mainPanel(
             fluidRow(
               column(6, plotOutput("prior_r")),
               column(6, htmlOutput("resultr"))),
             uiOutput("Optional_Plots_r")

           )
         )
)
,
tabPanel(em("\\(\\text{Regression}\\)"), withMathJax(),
         sidebarLayout(
           sidebarPanel(
             fluidRow(
               column(width = 4,
                      prettyRadioButtons(
                        inputId = "ANOREG",
                        label = em("\\(\\text{Type of Analysis}\\)"),
                        choices = list(
                          "\\(\\text{ANOVA}\\)" = 1,
                          "\\(\\text{Regression}\\)" = 2
                        ),
                        inline = TRUE,
                        selected = 1
                      )
               ),
               column(width = 6,
                      prettyRadioButtons(
                        inputId = "Modef",
                        label = "\\(\\text{Select Mode}\\)",
                        choices = list(
                          "\\(\\text{Sample size determination}\\)" = 1,
                          "\\(\\text{Fixed N}\\)" = 2,
                          "\\(\\text{BF calculator}\\)" = 3
                        ),
                        selected = 1,
                        inline = TRUE
                      )
               )
             ),

             # Regression model inputs
             conditionalPanel(
               condition = "input.ANOREG == 2&(input.Modef ==1|input.Modef ==2)",
               fluidRow(
                 column(width = 6,
                        sliderInput("pf", "\\(p\\text{ predictor - reduced model :}\\)",
                                    min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                 ),
                 column(width = 6,
                        sliderInput("kf", "\\(k\\text{ predictor - full model :}\\)",
                                    min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                 )
               )
             ),

             # ANOVA model structure
             conditionalPanel(
               condition = "input.ANOREG == 1&(input.Modef ==1|input.Modef ==2)",
               fluidRow(
                 column(width = 6,
                        prettyRadioButtons(
                          inputId = "redf",
                          label = em("\\(\\text{Reduced model}\\)"),
                          choices = list(
                            "\\(\\text{Intercept}\\)" = 1,
                            "\\(\\text{One-factor }\\)" = 2,
                            "\\(\\text{Two-factor - main effect}\\)" = 3
                          ),
                          inline = TRUE,
                          selected = 1
                        )
                 ),
                 column(width = 6,
                        conditionalPanel(
                          condition = "input.redf == 1",
                          prettyRadioButtons("full1", "\\(\\text{Full Model}\\)",
                                             choices = list(
                                               "\\(\\text{One-factor}\\)" = 2,
                                               "\\(\\text{Two-factor - main effect}\\)" = 3,
                                               "\\(\\text{Two-factor - interaction}\\)" = 4
                                             ),
                                             selected = 2,
                                             inline = TRUE
                          )
                        ),
                        conditionalPanel(
                          condition = "input.redf == 2",
                          prettyRadioButtons("full2", "\\(\\text{Full Model}\\)",
                                             choices = list(
                                               "\\(\\text{Two-factor - main effect}\\)" = 3
                                             ),
                                             selected = 3,
                                             inline = TRUE
                          )
                        ),
                        conditionalPanel(
                          condition = "input.redf == 3",
                          prettyRadioButtons("full3", "\\(\\text{Full Model}\\)",
                                             choices = list(
                                               "\\(\\text{Two-factor - interaction}\\)" = 4
                                             ),
                                             selected = 4,
                                             inline = TRUE
                          )
                        )
                 )
               )
             ),

             # Factor levels
             conditionalPanel(
               condition = "input.ANOREG == 1&(input.Modef ==1|input.Modef ==2)",
               fluidRow(
                 column(width = 6,
                        sliderInput("f1", "\\(\\text{Factor 1 level:}\\)", min = 2, max = 10, value = 2, step = 1, ticks = FALSE)
                 ),
                 column(width = 6,
                        sliderInput("f2", "\\(\\text{Factor 2 level:}\\)", min = 2, max = 10, value = 2, step = 1, ticks = FALSE)
                 )
               )
             ),

             # Hypotheses
             fluidRow(
               column(width = 6,
                      prettyRadioButtons(
                        inputId = "h0f",
                        label = em("$$\\mathcal{H}_0:$$"),
                        choices = list(
                          "\\(\\lambda^2 = 0 \\)" = 1,
                          "\\(\\lambda^2 \\in \\{0, \\epsilon\\}\\)" = 2
                        ),
                        inline = TRUE,
                        selected = 1
                      )
               ),
               column(width = 6,
                      conditionalPanel(
                        condition = "input.h0f == 1",
                        prettyRadioButtons(
                          inputId = "h1f",
                          label = em("\\(\\mathcal{H}_1:\\)"),
                          choices = list("\\(\\lambda^2 > 0 \\)" = 1),
                          inline = TRUE,
                          selected = 1
                        )
                      ),
                      conditionalPanel(
                        condition = "input.h0f == 2",
                        prettyRadioButtons(
                          inputId = "h1fe",
                          label = em("\\(\\mathcal{H}_1:\\)"),
                          choices = list("\\(\\lambda^2 > \\epsilon\\)" = 1),
                          inline = TRUE,
                          selected = 1
                        ),
                        sliderInput("epsilinff", "\\(\\epsilon :\\)", min = 0.01, max = 0.25, value = 0.1, step = 0.01, ticks = FALSE)
                      )
               )
             ),

             # Prior
             prettyRadioButtons(
               inputId = "modelf",
               label = em("\\(\\text{Analysis Prior Distribution}\\)"),
               choices = list(
                 "\\( \\text{Effect size prior} \\)" = 1,
                 "\\( \\text{Moment prior (must df ≥ 3)} \\)" = 2
               ),
               inline = TRUE,
               selected = 1
             ),

             fluidRow(
               column(width = 4,
                      conditionalPanel(condition = "input.modelf == 1",
                                       sliderInput("rf", "\\(\\text{r scale:}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE)
                      )
               ),
               column(width = 4,
                      sliderInput("fsdf", "\\(\\mathcal{f}^2 :\\)", min = 0.01, max = 0.5, value = 0.1, step = 0.01, ticks = FALSE)
               ),
               column(width = 4,
                      sliderInput("dff", "\\(\\text{df :}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
               )
             ),

             # Design prior
             conditionalPanel(
               condition = "input.Modef == 1|input.Modef == 2",
               prettyRadioButtons(
                 "priorf",
                 "\\( \\text{Design prior is the same as analysis prior: } \\)",
                 choices = list(
                   "\\( \\text{Yes} \\)" = 1,
                   "\\( \\text{No} \\)"  = 2
                 ),
                 selected = 1,
                 inline = TRUE
               )
             ),

             # Design prior ≠ analysis prior
             conditionalPanel(
               condition = "input.priorf == 2",
               prettyRadioButtons(
                 inputId = "modelfd",
                 label = em("\\(\\text{Design Prior Distribution}\\)"),
                 choices = list(
                   "\\( \\text{Effect size prior} \\)" = 1,
                   "\\( \\text{Moment prior (must df ≥ 3)} \\)" = 2,
                   "\\( \\text{Point} \\)" = 3
                 ),
                 inline = TRUE,
                 selected = 1
               ),
               fluidRow(
                 column(width = 4,
                        conditionalPanel(condition = "input.modelfd == 1",
                                         sliderInput("rfd", "\\(\\text{r scale:}\\)", min = 0, max = 3, value = 1, step = 0.01, ticks = FALSE)
                        ),
                        conditionalPanel(condition = "input.modelfd == 3",
                                         sliderInput("lfd", "\\(\\lambda^2:\\)", min = 0, max = 0.5, value = 0.1, step = 0.01, ticks = FALSE)
                        )
                 ),
                 column(width = 4,
                        conditionalPanel(condition = "input.modelfd == 1|input.modelfd == 2",
                                         sliderInput("fsdfd", "\\(\\mathcal{f}^2 :\\)", min = 0.01, max = 0.5, value = 0.1, step = 0.01, ticks = FALSE)
                        )
                 ),
                 column(width = 4,
                        conditionalPanel(condition = "input.modelfd == 1|input.modelfd == 2",
                                         sliderInput("dffd", "\\(\\text{df :}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                        )
                 )
               )
             ),

             # Sample size determination controls
             conditionalPanel(condition = "input.Modef == 1",
                              fluidRow(
                                column(width = 6,
                                       sliderInput("powerf", "\\( \\text{True Positive Evidence:} \\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)
                                ),
                                column(width = 6,
                                       sliderInput("alphaf", "\\( \\text{False Positive Evidence:} \\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE)
                                )
                              )
             ),

             conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",
                              sliderInput("bff", "\\( \\text{Bound of compelling evidence:} \\)", min = 1, max = 20, value = 3, ticks = FALSE)
             ),
             conditionalPanel(condition = "input.Modef == 2",
                              numericInput("nf", "\\( \\text{Sample Size } N: \\)", value = 50)
             ),
             conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",
                              actionButton("runf", label = "\\( \\text{Run} \\)"),
                              conditionalPanel(condition = "input.Modef == 1",
                                               em(span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;"))
                              )
             ),

             # BF calculator mode
             conditionalPanel(condition = "input.Modef == 3",
                              fluidRow(
                                column(width = 4, numericInput("df1f", label = "\\( \\mathcal{df}_1: \\)", value = 1)),
                                column(width = 4, numericInput("df2f", label = "\\( \\mathcal{df}_2: \\)", value = 30)),
                                column(width = 4, numericInput("fval", label = "\\( f\\text{-value:} \\)", value = 1))
                              ),
                              actionButton("calf", label = "\\( \\text{Calculate} \\)"),
                              htmlOutput("BFcalf")
             ),

             em("\\( \\text{Recommended hyperparameters:} \\)"),
             htmlOutput("prior_suggest"),
             conditionalPanel(
               condition = "input.Modef == 1",
               checkboxGroupInput(
                 "o_plot_f",
                 label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
                 choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
                 selected = NULL
               ),
               downloadButton("export_f", "Download result as PDF")
             )
           ),

           # Main Panel
           mainPanel(fluidRow(
             column(6, plotOutput("priorff")),
             column(6, htmlOutput("resultf"))

           ),uiOutput("Optional_Plots_f")
           )
         )
)
,


navbarMenu(
  "\\(\\text{Proportion}\\)",

  tabPanel("\\(\\text{One proportion - binomial}\\)", withMathJax(),
           sidebarLayout(
             sidebarPanel(
               prettyRadioButtons(
                 "Modebin", "\\(\\text{Select Mode}\\)",
                 choices = list(
                   "\\(\\text{Sample size determination}\\)" = 1,
                   "\\(\\text{Fixed N}\\)" = 2,
                   "\\(\\text{BF calculator}\\)" = 3
                 ),
                 selected = 1, inline = TRUE
               ),

               fluidRow(
                 column(6,
                        prettyRadioButtons(
                          inputId = "h0bin", label = em("\\(\\mathcal{H}_0:\\)"),
                          choices = list(
                            "\\(p = p_0\\)" = 1,
                            "\\(p \\in (p_0 - \\epsilon, p_0 + \\epsilon)\\)" = 2
                          ),
                          selected = 1, inline = TRUE
                        )),

                 column(6,
                        conditionalPanel(
                          condition = "input.h0bin == 1",
                          prettyRadioButtons(
                            inputId = "h1bin", label = em("\\(\\mathcal{H}_1:\\)"),
                            choices = list(
                              "\\(p ≠ p_0\\)" = 1,
                              "\\(p > p_0\\)" = 2,
                              "\\(p < p_0\\)" = 3
                            ),
                            selected = 1, inline = TRUE
                          )
                        ),

                        conditionalPanel(
                          condition = "input.h0bin == 2",
                          prettyRadioButtons(
                            inputId = "h1bine", label = em("\\(\\mathcal{H}_1:\\)"),
                            choices = list(
                              "\\(p \\not\\in \\{p_0 -\\epsilon, p_0 +\\epsilon\\}\\)" = 1,
                              "\\(p > p_0 + \\epsilon\\)" = 2,
                              "\\(p < p_0 - \\epsilon\\)" = 3
                            ),
                            selected = 1, inline = TRUE
                          )
                        )
                 )
               ),

               fluidRow(
                 column(4, sliderInput("h0prop", "\\(p_0\\)", min = .01, max = .99, value = .5, step = .01, ticks = FALSE)),

                 column(4,
                        conditionalPanel("input.h0bin == 2 && (input.h1bine == 1 || input.h1bine == 3)",
                                         sliderInput("lbbine", "\\(-\\epsilon\\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE),
                                         htmlOutput("bin_lower"))
                 ),

                 column(4,
                        conditionalPanel("input.h0bin == 2 && (input.h1bine == 1 || input.h1bine == 2)",
                                         sliderInput("ubbine", "\\(\\epsilon\\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE),
                                         htmlOutput("bin_upper"))
                 )
               ),

               prettyRadioButtons(
                 inputId = "modelbin",
                 label = em("\\(\\text{ Analysis Prior Distribution}\\)"),
                 choices = list("\\( \\text{Beta} \\)" = 1, "\\( \\text{Moment} \\)" = 2),
                 selected = 1, inline = TRUE
               ),

               fluidRow(
                 column(4,
                        conditionalPanel("input.modelbin == 1",
                                         sliderInput("alphabin", "\\(\\alpha\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)),
                        conditionalPanel("input.modelbin == 2",
                                         sliderInput("sbin", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE))
                 ),

                 column(4,
                        conditionalPanel("input.modelbin == 1",
                                         sliderInput("betabin", "\\(\\beta\\)", min = 0.01, max = 100, value = 1, step = 0.01, ticks = FALSE))
                 )
               ),

               conditionalPanel("input.Modebin == 1 || input.Modebin == 2",
                                prettyRadioButtons(
                                  "priorbin", "\\( \\text{Design prior is the same as analysis prior: } \\)",
                                  choices = list("\\( \\text{Yes} \\)" = 1, "\\( \\text{No} \\)" = 2),
                                  selected = 1, inline = TRUE
                                ),

                                conditionalPanel("input.priorbin == 2",
                                                 prettyRadioButtons(
                                                   "modelbind", em("\\( \\text{Design prior distribution} \\)"),
                                                   choices = list(
                                                     "\\( \\text{Beta} \\)" = 1,
                                                     "\\( \\text{Moment} \\)" = 2,
                                                     "\\( \\text{Point} \\)" = 3
                                                   ),
                                                   selected = 1, inline = TRUE
                                                 ),

                                                 conditionalPanel("input.modelbind == 3",
                                                                  sliderInput("h0bind", "\\(p_1\\)", min = .01, max = .99, value = .5, step = .01, ticks = FALSE)),

                                                 fluidRow(
                                                   column(4,
                                                          conditionalPanel("input.modelbind == 1",
                                                                           sliderInput("alphabind", "\\(\\alpha\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)),

                                                          conditionalPanel("input.modelbind == 2",
                                                                           sliderInput("sbind", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE))
                                                   ),
                                                   column(4,
                                                          conditionalPanel("input.modelbind == 1",
                                                                           sliderInput("betabind", "\\(\\beta\\)", min = 0.01, max = 100, value = 1, step = 0.01, ticks = FALSE))
                                                   )
                                                 )
                                )
               ),

               conditionalPanel(
                 "input.modelbind == 3 && (input.h1bin == 2 || input.h1bin == 3 || input.h1bine == 2 || input.h1bine == 3)",
                 em(span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
               ),

               conditionalPanel("input.Modebin == 1",
                                em("\\(\\text{Controlling for the probability of}\\)"),
                                fluidRow(
                                  column(6, sliderInput("powerbin", "\\(\\text{True Positive Evidence:}\\)", min = .5, max = .99, value = .8, step = .01, ticks = FALSE)),
                                  column(6, sliderInput("FP_bin", "\\(\\text{False Positive Evidence:}\\)", min = .001, max = .05, value = .05, step = .001, ticks = FALSE))
                                )
               ),

               conditionalPanel("input.Modebin == 1 || input.Modebin == 2",
                                sliderInput("bbin", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)),

               conditionalPanel("input.Modebin == 2 || input.Modebin == 3",
                                fluidRow(
                                  column(6, numericInput("nbin", "\\(\\text{Sample Size:}\\)", value = 50)),
                                  column(6,
                                         conditionalPanel("input.Modebin == 3",
                                                          numericInput("xbin", "\\(\\text{Number of Success:}\\)", value = 25)))
                                )
               ),

               conditionalPanel("input.Modebin == 1 || input.Modebin == 2",
                                actionButton("runbin", label = "\\(\\text{Run}\\)"),
                                conditionalPanel("input.Modebin == 1",
                                                 em(span("\\(\\text{Note: Error when the required N > 10,000}\\)", style = "color: red;")))
               ),

               conditionalPanel("input.Modebin == 3",
                                actionButton("calbin", label = "\\(\\text{Calculate}\\)"),
                                htmlOutput("BFbin")),
               conditionalPanel(
                 condition = "input.Modebin == 1",
                 checkboxGroupInput(
                   "o_plot_bin",
                   label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
                   choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
                   selected = NULL
                 ),
                 downloadButton("export_bin", "Download result as PDF")
               )


             ),

             mainPanel(fluidRow(
               column(6, plotOutput("prior_bin")),
               column(6, htmlOutput("resultbin"))
             ),
             uiOutput("Optional_Plots_bin"))
           )
  ), tabPanel("\\(\\text{Two proportion}\\)",withMathJax(),
sidebarLayout(sidebarPanel(
  prettyRadioButtons(
    "Modep2", "\\(\\text{Select Mode}\\)",
    choices = list(
      "\\(\\text{Sample size determination}\\)" = 1,
      "\\(\\text{Fixed N}\\)" = 2,
      "\\(\\text{BF calculator}\\)" = 3
    ),
    selected = 1, inline = TRUE
  ) ,fluidRow(
    column(5,
           em("\\(\\mathcal{H}_0: p_0 = p_1 = p_2\\)"),
           br(),
           em("\\(p_0 \\sim \\text{Beta}(\\alpha_0, \\beta_0)\\)")




           ),




    column(5,
           em("\\(\\mathcal{H}_1: p_1 \\neq p_2\\)"),
           br(),
           em("\\(p_1 \\sim \\text{Beta}(\\alpha_1, \\beta_1)\\)"),
           br(),
           em("\\(p_2 \\sim \\text{Beta}(\\alpha_2, \\beta_2)\\)")
    )
  ),fluidRow(

  column(4,
  sliderInput("alpha0", "\\(\\alpha_0\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
  sliderInput("beta0", "\\(\\beta_0\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)


           ),
  column(4,
         sliderInput("alpha1", "\\(\\alpha_1\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
         sliderInput("beta1", "\\(\\beta_1\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)
           ),
  column(4,
         sliderInput("alpha2", "\\(\\alpha_2\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
         sliderInput("beta2", "\\(\\beta_2\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)

  )

  ),
  conditionalPanel("input.Modep2 == 1 | input.Modep2 == 2",
  prettyRadioButtons(
    "priorp2", "\\( \\text{Design prior is the same as analysis prior: } \\)",
    choices = list("\\( \\text{Yes} \\)" = 1, "\\( \\text{No} \\)" = 2),
    selected = 1, inline = TRUE
  )),

  conditionalPanel(
    condition = "input.priorp2 == 2",
  fluidRow(
    column(4,
           prettyRadioButtons(
             "model_p1",
             label = "\\( \\text{Model for } p_1: \\)",
             choices = list(
               "\\( \\text{Fixed } p_1 \\)" = 1,
               "\\( p_1 \\sim \\text{Beta}(\\alpha_{1d}, \\beta_{1d}) \\)" = 2
             ),
             selected = 1,
             inline = TRUE
           ),
           prettyRadioButtons(
             "model_p2",
             label = "\\( \\text{Model for } p_2: \\)",
             choices = list(
               "\\( \\text{Fixed } p_1 \\)" = 1,
               "\\( p_1 \\sim \\text{Beta}(\\alpha_{2d}, \\beta_{2d}) \\)" = 2
             ),
             selected = 1,
             inline = TRUE
           )



  ),

  column(4,
         conditionalPanel(
           condition = "input.model_p1 == 1",
           sliderInput("location1d", "\\(p_1 =\\)", min = 0.01, max = 0.99, value = 0.5, step = 0.01, ticks = FALSE)
         ),
         conditionalPanel(
           condition = "input.model_p1 == 2",
           sliderInput("alpha1d", "\\(\\alpha_{1d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
           sliderInput("beta1d", "\\(\\beta_{1d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)
         )
  )
  ,
  column(4,
         conditionalPanel(
           condition = "input.model_p2 == 1",
           sliderInput("location2d", "\\(p_2 =\\)", min = 0.01, max = 0.99, value = 0.5, step = 0.01, ticks = FALSE)
         ),
         conditionalPanel(
           condition = "input.model_p2 == 2",
           sliderInput("alpha2d", "\\(\\alpha_{2d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
           sliderInput("beta2d", "\\(\\beta_{2d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)
         )
  )



  )),

  conditionalPanel("input.Modep2 == 1",
                   em("\\(\\text{Controlling for the probability of}\\)"),
                   fluidRow(
                     column(12, sliderInput("powerp2", "\\(\\text{True Positive Evidence:}\\)", min = .5, max = .99, value = .8, step = .01, ticks = FALSE)),
                     #column(6, sliderInput("FP_p2", "\\(\\text{False Positive Evidence:}\\)", min = .001, max = .05, value = .05, step = .001, ticks = FALSE))
                   )
  ),
  conditionalPanel("input.Modep2 == 1 || input.Modep2 == 2",
                   sliderInput("bp2", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)),

  conditionalPanel("input.Modep2 == 2 ||input.Modep2 == 3 ",
                   em("\\( \\text{Sample size per group} \\)"),
                   fluidRow(
                     column(6, numericInput("n1p2", "\\(n_1\\)", value = 50)),
                     column(6, numericInput("n2p2", "\\(n_2\\)", value = 50))
                   )

                   ),
  conditionalPanel("input.Modep2 == 3",
                   em("\\( \\text{Number of success} \\)"),
                   fluidRow(
                     column(6, numericInput("x1p2", "\\(x_1\\)", value = 50)),
                     column(6, numericInput("x2p2", "\\(x_2\\)", value = 50))
                   )) ,

  conditionalPanel("input.Modep2 == 1 || input.Modep2 == 2",
                   actionButton("runp2", label = "\\(\\text{Run}\\)"),
                   conditionalPanel("input.Modep2 == 1", em(span("\\(\\text{Note: Error when required } N > 5,000\\)", style = "color: red;")))
  ),
  conditionalPanel("input.Modep2 == 3",
                   actionButton("calp2", label = "\\(\\text{Calculate}\\)"),
                   htmlOutput("BFp2")
  ),

  conditionalPanel(
    condition = "input.Modep2 == 1",
    checkboxGroupInput(
      "o_plot_p2",
      label = "\\(\\text{Additional Plots (computationally very intensive):}\\)",
      choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
      selected = NULL
    ),
    downloadButton("export_p2", "Download result as PDF")
  )
),mainPanel(

  fluidRow(
    column(6, plotOutput("prior_p0")),
    column(6, htmlOutput("resultp2"))
  ),
  fluidRow(
    column(6, plotOutput("prior_p1")),
    column(6, plotOutput("prior_p2"))
  ),
  uiOutput("Optional_Plots_p2")

))







)))

# Server logic
server <- function(input, output, session) {

  server_t1(input, output, session)
  server_t2(input, output, session)
  server_r(input, output, session)
  server_f(input, output, session)
  server_bin(input, output, session)
  server_p2(input, output, session)

}

# Run the application

#' Launch BayesPower Shiny App
#'
#' @export
BayesPower_testing <- function(){

pkgs <- c(
  "rootSolve", "shiny", "gsl", "shinyWidgets", "shinyjs",
  "kableExtra", "knitr", "Rcpp", "fontawesome", "BH",
  "bslib", "pracma", "profvis", "mombf", "ExtDist",
  "ggplot2", "patchwork"
)

invisible(lapply(pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}))



  # Run the application
  shinyApp(ui = ui, server = server)

}
