

ui <-
  shiny::navbarPage(id = "id",
                 "\\(\\text{BayesPower}_{1.0}\\)",
  shiny::navbarMenu(
    "\\(\\text{Standardized Mean Difference}\\)",

    shiny::tabPanel(
      "\\(\\text{One-sample/paired t-test}\\)",

      # Custom font for MathJax
      shiny::tags$style(shiny::HTML("
    body {
  font-family: sans-serif;
}
  ")),

      shiny::withMathJax(),

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          # Mode selection
          shinyWidgets::prettyRadioButtons(
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
          shiny::fluidRow(
            shiny::column(6,
                   shinyWidgets::prettyRadioButtons(
                     inputId = "h0t1",
                     label = shiny::em("\\(\\mathcal{H}_0:\\)"),
                     choices = list(
                       "\\(\\delta = 0\\)" = 1,
                       "\\(\\delta \\in \\{-\\epsilon, \\epsilon\\}\\)" = 2
                     ),
                     inline = TRUE,
                     selected = 1
                   )
            ),
            shiny::column(6,
                   shiny::conditionalPanel("input.h0t1 == 1",
                                    shinyWidgets::prettyRadioButtons(
                                      inputId = "h1t1",
                                      label = shiny::em("\\(\\mathcal{H}_1:\\)"),
                                      choices = list(
                                        "\\(\\delta ≠ 0\\)" = 1,
                                        "\\(\\delta > 0\\)" = 2,
                                        "\\(\\delta < 0\\)" = 3
                                      ),
                                      inline = TRUE,
                                      selected = 1
                                    )
                   ),
                   shiny::conditionalPanel("input.h0t1 == 2",
                                    shinyWidgets::prettyRadioButtons(
                                      inputId = "h1t1e",
                                      label = shiny::em("\\(\\mathcal{H}_1:\\)"),
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
          shiny::fluidRow(
            shiny::column(6,
                   shiny::conditionalPanel("input.h1t1e == 2 && input.h0t1 == 2",
                                    shiny::em("\\( -\\epsilon = 0 \\)")
                   ),
                   shiny::conditionalPanel("(input.h1t1e == 1 || input.h1t1e == 3) && input.h0t1 == 2",
                                    shiny::sliderInput("lbt1e", "\\( -\\epsilon \\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE)
                   )
            ),
            shiny::column(6,
                   shiny::conditionalPanel("input.h1t1e == 3 && input.h0t1 == 2",
                                    shiny::em("\\( \\epsilon = 0 \\)")
                   ),
                   shiny::conditionalPanel("(input.h1t1e == 1 || input.h1t1e == 2) && input.h0t1 == 2",
                                    shiny::sliderInput("ubt1e", "\\( \\epsilon \\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE)
                   )
            )
          ),

          # Analysis prior
          shinyWidgets::prettyRadioButtons(
            inputId = "modelt1",
            label = shiny::em("\\(\\text{Analysis Prior Distribution}\\)"),
            choices = list(
              "\\(\\text{Scaled t}\\)" = 1,
              "\\(\\text{Normal}\\)" = 2,
              "\\(\\text{Moment}\\)" = 3
            ),
            inline = TRUE,
            selected = 1
          ),

          # Location/scale/df
          shiny::fluidRow(
            shiny::column(4,
                   shiny::conditionalPanel("input.h0t1 == 1",
                                    shiny::sliderInput("lt1", "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
                   ),
                   shiny::conditionalPanel("input.h0t1 == 2", shiny::em("\\(\\text{Location = 0}\\)"))
            ),
            shiny::column(4,
                   shiny::sliderInput("st1", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 0.707, step = 0.001, ticks = FALSE)
            ),
            shiny::column(4,
                   shiny::conditionalPanel("input.modelt1 == 1",
                                    shiny::sliderInput("dft1", "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                   )
            )
          ),

          # Design prior if different
          shiny::conditionalPanel("input.Modet1 == 1 || input.Modet1 == 2",
                           shinyWidgets::prettyRadioButtons(
                             "prior",
                             "\\(\\text{Design prior is the same as analysis prior:}\\)",
                             choices = list("\\(\\text{Yes}\\)" = 1, "\\(\\text{No}\\)" = 2),
                             selected = 1,
                             inline = TRUE
                           ),
                           shiny::conditionalPanel("input.prior == 2",
                                            shinyWidgets::prettyRadioButtons(
                                              inputId = "modelt1d",
                                              label = shiny::em("\\(\\text{Design Prior Distribution}\\)"),
                                              choices = list(
                                                "\\(\\text{Scaled t}\\)" = 1,
                                                "\\(\\text{Normal}\\)" = 2,
                                                "\\(\\text{Moment}\\)" = 3,
                                                "\\(\\text{Point}\\)" = 4
                                              ),
                                              selected = 1,
                                              inline = TRUE
                                            ),
                                            shiny::fluidRow(
                                              shiny::column(4,
                                                     shiny::conditionalPanel("input.h0t1 == 1 || input.modelt1d == 4",
                                                                      shiny::sliderInput("lt1d", "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
                                                     ),
                                                     shiny::conditionalPanel("input.h0t1 == 2 && input.modelt1d != 4", shiny::em("\\(\\text{Location = 0}\\)"))
                                              ),
                                              shiny::column(4,
                                                     shiny::conditionalPanel("input.modelt1d != 4",
                                                                      shiny::sliderInput("st1d", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 0.707, step = 0.001, ticks = FALSE)
                                                     )
                                              ),
                                              shiny::column(4,
                                                     shiny::conditionalPanel("input.modelt1d == 1",
                                                                      shiny::sliderInput("dft1d", "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                                                     )
                                              )
                                            )
                           )
          ),

          # Note for point prior
          shiny::conditionalPanel("input.modelt1d == 4 && (input.h1t1 == 2 || input.h1t1 == 3 || input.h1t1e == 2 || input.h1t1e == 3)",
                           shiny::em(shiny::span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
          ),

          # Power / error control
          shiny::conditionalPanel("input.Modet1 == 1",
                           shiny::em("\\(\\text{Controlling for the probability of}\\)"),
                           shiny::fluidRow(
                             shiny::column(6,
                                    shiny::sliderInput("powert1", "\\(\\text{True Positive Evidence:}\\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)
                             ),
                             shiny::column(6,
                                    shiny::sliderInput("alphat1", "\\(\\text{False Positive Evidence:}\\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE)
                             )
                           )
          ),

          # Bound
          shiny::conditionalPanel("input.Modet1 == 1 || input.Modet1 == 2",
                           shiny::sliderInput("bt1", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)
          ),

          # Sample size input
          shiny::conditionalPanel("input.Modet1 == 2",
                           shiny::numericInput("nt1", "\\(\\text{Sample Size:}\\)", value = 50)
          ),

          # Run button + error message
          shiny::conditionalPanel("input.Modet1 == 1 || input.Modet1 == 2",
                           shiny::actionButton("runt1", label = "\\(\\text{Run}\\)"),
                           shiny::conditionalPanel("input.Modet1 == 1",
                                            shiny::em(shiny::span("\\(\\text{Note: Error when the required N > 10,000}\\)", style = "color: red;"))
                           )
          ),

          # BF calculator mode
          shiny::conditionalPanel("input.Modet1 == 3",
                           shiny::fluidRow(
                             shiny::column(6, shiny::numericInput("t1df", "\\(\\text{Degree of freedom:}\\)", value = 50)),
                             shiny::column(6, shiny::numericInput("t1tval", "\\(\\text{t-value:}\\)", value = 2))
                           ),
                           shiny::actionButton("cal1", label = "\\(\\text{Calculate}\\)"),
                           shiny::htmlOutput("BFt1")

          ),

          shiny::conditionalPanel(
            condition = "input.Modet1 == 1",
            shiny::checkboxGroupInput(
              "o_plot_t1",
              label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
              choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
              selected = NULL
            ),
            shiny::downloadButton("export_t1", "Download result as PDF")
          )
        ),

        # Main panel with result tabs
        shiny::mainPanel( shiny::fluidRow(
          shiny::column(6, shiny::plotOutput("priort1")),
          shiny::column(6, shiny::htmlOutput("resultt1"))),
          shiny::uiOutput("Optional_Plots_t1"))
      )
    )
    ,shiny::tabPanel(
      "\\(\\text{Independent samples t-test}\\)",
      shiny::withMathJax(),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          # Mode selection
          shinyWidgets::prettyRadioButtons(
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
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shinyWidgets::prettyRadioButtons(
                inputId = "h0t2",
                label = shiny::em("\\(\\mathcal{H}_0:\\)"),
                choices = list(
                  "\\(\\delta = 0\\)" = 1,
                  "\\(\\delta \\in \\{-\\epsilon, \\epsilon\\}\\)" = 2
                ),
                selected = 1,
                inline = TRUE
              )
            ),
            shiny::column(
              width = 6,
              shiny::conditionalPanel(
                condition = "input.h0t2 == 1",
                shinyWidgets::prettyRadioButtons(
                  inputId = "h1t2",
                  label = shiny::em("\\(\\mathcal{H}_1:\\)"),
                  choices = list(
                    "\\(\\delta ≠ 0\\)" = 1,
                    "\\(\\delta > 0\\)" = 2,
                    "\\(\\delta < 0\\)" = 3
                  ),
                  selected = 1,
                  inline = TRUE
                )
              ),
              shiny::conditionalPanel(
                condition = "input.h0t2 == 2",
                shinyWidgets::prettyRadioButtons(
                  inputId = "h1t2e",
                  label = shiny::em("\\(\\mathcal{H}_1:\\)"),
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
          shiny::fluidRow(
            shiny::column(
              width = 6,
              shiny::conditionalPanel("input.h1t2e == 2 && input.h0t2 == 2", shiny::em("\\(-\\epsilon = 0\\)")),
              shiny::conditionalPanel(
                condition = "(input.h1t2e == 1 || input.h1t2e == 3) && input.h0t2 == 2",
                shiny::sliderInput("lbt2e", label = "\\(-\\epsilon\\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE)
              )
            ),
            shiny::column(
              width = 6,
              shiny::conditionalPanel("input.h1t2e == 3 && input.h0t2 == 2", shiny::em("\\(\\epsilon = 0\\)")),
              shiny::conditionalPanel(
                condition = "(input.h1t2e == 1 || input.h1t2e == 2) && input.h0t2 == 2",
                shiny::sliderInput("ubt2e", label = "\\(\\epsilon\\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE)
              )
            )
          ),

          # Analysis prior
          shinyWidgets::prettyRadioButtons(
            inputId = "modelt2",
            label = shiny::em("\\(\\text{Analysis Prior Distribution}\\)"),
            choices = list(
              "\\(\\text{Scaled t}\\)" = 1,
              "\\(\\text{Normal}\\)" = 2,
              "\\(\\text{Moment}\\)" = 3
            ),
            selected = 1,
            inline = TRUE
          ),

          shiny::fluidRow(
            shiny::column(
              width = 4,
              shiny::conditionalPanel(
                "input.h0t2 == 1",
                shiny::sliderInput("lt2", label = "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
              ),
              shiny::conditionalPanel("input.h0t2 == 2", shiny::em("\\(\\text{Location = 0}\\)"))
            ),
            shiny::column(
              width = 4,
              shiny::sliderInput("st2", label = "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = .707, step = 0.001, ticks = FALSE)
            ),
            shiny::column(
              width = 4,
              shiny::conditionalPanel(
                "input.modelt2 == 1",
                shiny::sliderInput("dft2", label = "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
              )
            )
          ),

          # Design prior toggle
          shiny::conditionalPanel(
            "input.Modet2 == 1 || input.Modet2 == 2",
            shinyWidgets::prettyRadioButtons(
              "priort2",
              "\\(\\text{Design prior is the same as analysis prior:}\\)",
              choices = list("\\(\\text{Yes}\\)" = 1, "\\(\\text{No}\\)" = 2),
              selected = 1,
              inline = TRUE
            ),
            shiny::conditionalPanel(
              "input.priort2 == 2",
              shinyWidgets::prettyRadioButtons(
                inputId = "modelt2d",
                label = shiny::em("\\(\\text{Design prior distribution}\\)"),
                choices = list(
                  "\\(\\text{Scaled t}\\)" = 1,
                  "\\(\\text{Normal}\\)" = 2,
                  "\\(\\text{Moment}\\)" = 3,
                  "\\(\\text{Point}\\)" = 4
                ),
                selected = 1,
                inline = TRUE
              ),
              shiny::fluidRow(
                shiny::column(
                  width = 4,
                  shiny::conditionalPanel(
                    "input.h0t2 == 1 || input.modelt2d == 4",
                    shiny::sliderInput("lt2d", label = "\\(\\text{Location}\\)", min = -2, max = 2, value = 0, step = 0.01, ticks = FALSE)
                  ),
                  shiny::conditionalPanel("input.h0t2 == 2 && input.modelt2d != 4", shiny::em("\\(\\text{Location = 0}\\)"))
                ),
                shiny::column(
                  width = 4,
                  shiny::conditionalPanel(
                    "input.modelt2d == 1 || input.modelt2d == 2 || input.modelt2d == 3",
                    shiny::sliderInput("st2d", label = "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE)
                  )
                ),
                shiny::column(
                  width = 4,
                  shiny::conditionalPanel(
                    "input.modelt2d == 1",
                    shiny::sliderInput("dft2d", label = "\\(\\text{df}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                  )
                )
              )
            )
          ),

          # Hint for point prior direction
          shiny::conditionalPanel(
            "input.modelt2d == 4 && (input.h1t2 == 2 || input.h1t2 == 3 || input.h1t2e == 2 || input.h1t2e == 3)",
            shiny::em(shiny::span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
          ),

          # Controls for sample size determination
          shiny::conditionalPanel("input.Modet2 == 1",
                           shiny::em("\\(\\text{Controlling for the probability of}\\)"),
                           shiny::fluidRow(
                             shiny::column(6, shiny::sliderInput("powert2", "\\(\\text{True Positive Evidence:}\\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)),
                             shiny::column(6, shiny::sliderInput("alphat2", "\\(\\text{False Positive Evidence:}\\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE))
                           )
          ),

          # Shared parameters
          shiny::conditionalPanel("input.Modet2 == 1 || input.Modet2 == 2", shiny::sliderInput("bt2", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)),
          shiny::conditionalPanel("input.Modet2 == 1", shiny::sliderInput("rt2", "\\(N_2/N_1\\)", min = 1, max = 10, value = 1, ticks = FALSE)),

          # Fixed sample sizes
          shiny::conditionalPanel("input.Modet2 == 2||input.Modet2 == 3",
                           shiny::fluidRow(
                             shiny::column(4, shiny::numericInput("n1t2", "\\(N_1:\\)", value = 50)),
                             shiny::column(4, shiny::numericInput("n2t2", "\\(N_2:\\)", value = 50)),
                             shiny::conditionalPanel("input.Modet2 == 3",
                                                     shiny::column(4, shiny::numericInput("t2tval", "\\(\\text{t-value:}\\)", value = 2)))
                             ),
                           shiny::conditionalPanel("input.Modet2 == 3",
                                                   shiny::actionButton("cal1", label = "\\(\\text{Calculate}\\)"),
                                                   shiny::htmlOutput("BFt2"))



          ),

          # Action buttons
          shiny::conditionalPanel("input.Modet2 == 1 || input.Modet2 == 2",
                           shiny::actionButton("runt2", label = "\\(\\text{Run}\\)"),
                           shiny::conditionalPanel("input.Modet2 == 1", shiny::em(shiny::span("\\(\\text{Note: Error when required } N > 10,000\\)", style = "color: red;")))
          ),
          shiny::conditionalPanel(
            condition = "input.Modet2 == 1",
            shiny::checkboxGroupInput(
              "o_plot_t2",
              label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
              choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
              selected = NULL
            ),
            shiny::downloadButton("export_t2", "Download result as PDF")
          )





        ),

        shiny::mainPanel(shiny::fluidRow(
          shiny::column(6, shiny::plotOutput("priort2")),
          shiny::column(6, shiny::htmlOutput("resultt2"))),
          shiny::uiOutput("Optional_Plots_t2")





          )
      )
    )
  )
,
shiny::tabPanel("\\(\\text{Correlation}\\)", shiny::withMathJax(),
         shiny::sidebarLayout(
           shiny::sidebarPanel(
             shinyWidgets::prettyRadioButtons(
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
             shiny::fluidRow(
               shiny::column(6,
                      shinyWidgets::prettyRadioButtons(
                        inputId = "h0r",
                        label = shiny::em("\\(\\mathcal{H}_0:\\)"),
                        choices = list(
                          "\\(\\rho = \\rho_0\\)" = 1,
                          "\\(\\rho \\in (\\rho_0 - \\epsilon, \\rho_0 + \\epsilon)\\)" = 2
                        ),
                        selected = 1,
                        inline = TRUE
                      )
               ),
               shiny::column(6,
                      shiny::conditionalPanel(
                        condition = "input.h0r == 1",
                        shinyWidgets::prettyRadioButtons(
                          inputId = "h1r",
                          label = shiny::em("\\(\\mathcal{H}_1:\\)"),
                          choices = list(
                            "\\(\\rho ≠ \\rho_0\\)" = 1,
                            "\\(\\rho > \\rho_0\\)" = 2,
                            "\\(\\rho < \\rho_0\\)" = 3
                          ),
                          selected = 1,
                          inline = TRUE
                        )
                      ),
                      shiny::conditionalPanel(
                        condition = "input.h0r == 2",
                        shinyWidgets::prettyRadioButtons(
                          inputId = "h1re",
                          label = shiny::em("\\(\\mathcal{H}_1:\\)"),
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
             shiny::fluidRow(
               shiny::column(4,
                      shiny::sliderInput("h0pho", "\\(\\rho_0\\)", min = -.99, max = .99, value = 0, step = 0.01, ticks = FALSE)
               ),
               shiny::column(4,
                      shiny::conditionalPanel(
                        condition = "(input.h1re == 1 || input.h1re == 3) && input.h0r == 2",
                        shiny::sliderInput("lbre", "\\( -\\epsilon \\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE),
                        shiny::htmlOutput("r_lower")
                      )
               ),
               shiny::column(4,
                      shiny::conditionalPanel(
                        condition = "(input.h1re == 1 || input.h1re == 2) && input.h0r == 2",
                        shiny::sliderInput("ubre", "\\(\\epsilon \\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE),
                        shiny::htmlOutput("r_upper")
                      )
               )
             ),

             # Analysis prior
             shinyWidgets::prettyRadioButtons(
               inputId = "modelr",
               label = shiny::em("\\(\\text{ Analysis Prior Distribution}\\)"),
               choices = list(
                 "\\( \\text{Default Stretched Beta} \\)" = 1,
                 "\\( \\text{Stretched Beta} \\)" = 2,
                 "\\( \\text{Moment} \\)" = 3
               ),
               selected = 1,
               inline = TRUE
             ),

             shiny::fluidRow(
               shiny::column(4,
                      shiny::conditionalPanel("input.modelr == 1",
                                       shiny::sliderInput("kr", "\\(k \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                      ),
                      shiny::conditionalPanel("input.modelr == 3",
                                       shiny::sliderInput("sr", "\\(Scale \\)", min = 0.01, max = 1, value = 0.01, step = 0.01, ticks = FALSE)
                      )
               ),
               shiny::column(4,
                      shiny::conditionalPanel("input.modelr == 2",
                                       shiny::sliderInput("ralpha", "\\(\\alpha \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                      )
               ),
               shiny::column(4,
                      shiny::conditionalPanel("input.modelr == 2",
                                       shiny::sliderInput("rbeta", "\\(\\beta \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                      )
               )
             ),

             # Design prior
             shiny::conditionalPanel("input.Moder == 1 | input.Moder == 2",
                              shinyWidgets::prettyRadioButtons(
                                inputId = "priorr",
                                label = "\\( \\text{Design prior is the same as analysis prior: } \\)",
                                choices = list(
                                  "\\( \\text{Yes} \\)" = 1,
                                  "\\( \\text{No} \\)" = 2
                                ),
                                selected = 1,
                                inline = TRUE
                              ),
                              shiny::conditionalPanel("input.priorr == 2",
                                               shinyWidgets::prettyRadioButtons(
                                                 inputId = "modelrd",
                                                 label = shiny::em("\\( \\text{Design prior distribution} \\)"),
                                                 choices = list(
                                                   "\\( \\text{Default Stretched Beta} \\)" = 1,
                                                   "\\( \\text{Stretched Beta} \\)" = 2,
                                                   "\\( \\text{Moment} \\)" = 3,
                                                   "\\( \\text{Point} \\)" = 4
                                                 ),
                                                 selected = 1,
                                                 inline = TRUE
                                               ),
                                               shiny::conditionalPanel("input.modelrd == 4",
                                                                shiny::sliderInput("h0phod", "\\(\\rho_1\\)", min = -1, max = 1, value = 0.3, step = 0.01, ticks = FALSE)
                                               ),
                                               shiny::fluidRow(
                                                 shiny::column(4,
                                                        shiny::conditionalPanel("input.modelrd == 1",
                                                                         shiny::sliderInput("rkd", "\\( k \\)", min = 0.01, max = 10, value = 1, step = 0.01, ticks = FALSE)
                                                        ),
                                                        shiny::conditionalPanel("input.modelrd == 3",
                                                                         shiny::sliderInput("rsd", "\\( Scale \\)", min = 0.01, max = 1, value = 0.1, step = 0.01, ticks = FALSE)
                                                        )
                                                 ),
                                                 shiny::column(4,
                                                        shiny::conditionalPanel("input.modelrd == 2",
                                                                         shiny::sliderInput("ralphad", "\\(\\alpha \\)", min = 0.01, max = 10, value = 0.1, step = 0.01, ticks = FALSE)
                                                        )
                                                 ),
                                                 shiny::column(4,
                                                        shiny::conditionalPanel("input.modelrd == 2",
                                                                         shiny::sliderInput("rbetad", "\\(\\beta \\)", min = 0.01, max = 10, value = 0.1, step = 0.01, ticks = FALSE)
                                                        )
                                                 )
                                               )
                              ),
                              shiny::conditionalPanel("input.modelrd == 4 & (input.h1r == 2 | input.h1r == 3 | input.h1re == 2 | input.h1re == 3)",
                                               shiny::em(shiny::span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
                              )
             ),

             # Power planning
             shiny::conditionalPanel("input.Moder == 1",
                              shiny::em("\\(\\text{Controlling for the probability of}\\)"),
                              shiny::fluidRow(
                                shiny::column(6,
                                       shiny::sliderInput("powerr", "\\( \\text{True Positive Evidence:} \\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)
                                ),
                                shiny::column(6,
                                       shiny::conditionalPanel("input.h0r == 1 | input.h0r == 2",
                                                        shiny::sliderInput("alphapr", "\\( \\text{False Positive Evidence:} \\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE)
                                       )
                                )
                              )
             ),

             # Common for Mode 1 and 2
             shiny::conditionalPanel("input.Moder == 1 | input.Moder == 2",
                              shiny::sliderInput("br", "\\( \\text{Bound of compelling evidence:} \\)", min = 1, max = 20, value = 3, ticks = FALSE)
             ),

             shiny::conditionalPanel("input.Moder == 2",
                              shiny::numericInput("nr", "\\( \\text{Sample Size:} \\)", value = 50)
             ),

             shiny::conditionalPanel("input.Moder == 1 | input.Moder == 2",
                              shiny::actionButton("runr", label = "\\( \\text{Run} \\)"),
                              shiny::conditionalPanel("input.Moder == 1",
                                               shiny::em(shiny::span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;"))
                              )
             ),

             # BF calculator
             shiny::conditionalPanel("input.Moder == 3",
                              shiny::fluidRow(
                                shiny::column(6,
                                       shiny::numericInput("rdf", "\\( \\text{Sample size:} \\)", value = 50)
                                ),
                                shiny::column(6,
                                       shiny::numericInput("rval", "\\(\\text{Pearson's:} \\)", value = 0)
                                )
                              ),
                              shiny::actionButton("calr", label = "\\( \\text{Calculate} \\)"),
                              shiny::htmlOutput("BFrv")
             ),
             shiny::conditionalPanel(
               condition = "input.Moder == 1",
               shiny::checkboxGroupInput(
                 "o_plot_r",
                 label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
                 choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
                 selected = NULL
               ),
               shiny::downloadButton("export_r", "Download result as PDF")
             )
           ),

           # Main panel with output tabs
           shiny::mainPanel(
             shiny::fluidRow(
               shiny::column(6, shiny::plotOutput("prior_r")),
               shiny::column(6, shiny::htmlOutput("resultr"))),
             shiny::uiOutput("Optional_Plots_r")

           )
         )
)
,
shiny::tabPanel(shiny::em("\\(\\text{Regression}\\)"), shiny::withMathJax(),
         shiny::sidebarLayout(
           shiny::sidebarPanel(
             shiny::fluidRow(
               shiny::column(width = 4,
                      shinyWidgets::prettyRadioButtons(
                        inputId = "ANOREG",
                        label = shiny::em("\\(\\text{Type of Analysis}\\)"),
                        choices = list(
                          "\\(\\text{ANOVA}\\)" = 1,
                          "\\(\\text{Regression}\\)" = 2
                        ),
                        inline = TRUE,
                        selected = 1
                      )
               ),
               shiny::column(width = 6,
                      shinyWidgets::prettyRadioButtons(
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
             shiny::conditionalPanel(
               condition = "input.ANOREG == 2&(input.Modef ==1|input.Modef ==2)",
               shiny::fluidRow(
                 shiny::column(width = 6,
                        shiny::sliderInput("pf", "\\(p\\text{ predictor - reduced model :}\\)",
                                    min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                 ),
                 shiny::column(width = 6,
                        shiny::sliderInput("kf", "\\(k\\text{ predictor - full model :}\\)",
                                    min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                 )
               )
             ),

             # ANOVA model structure
             shiny::conditionalPanel(
               condition = "input.ANOREG == 1&(input.Modef ==1|input.Modef ==2)",
               shiny::fluidRow(
                 shiny::column(width = 6,
                        shinyWidgets::prettyRadioButtons(
                          inputId = "redf",
                          label = shiny::em("\\(\\text{Reduced model}\\)"),
                          choices = list(
                            "\\(\\text{Intercept}\\)" = 1,
                            "\\(\\text{One-factor }\\)" = 2,
                            "\\(\\text{Two-factor - main effect}\\)" = 3
                          ),
                          inline = TRUE,
                          selected = 1
                        )
                 ),
                 shiny::column(width = 6,
                        shiny::conditionalPanel(
                          condition = "input.redf == 1",
                          shinyWidgets::prettyRadioButtons("full1", "\\(\\text{Full Model}\\)",
                                             choices = list(
                                               "\\(\\text{One-factor}\\)" = 2,
                                               "\\(\\text{Two-factor - main effect}\\)" = 3,
                                               "\\(\\text{Two-factor - interaction}\\)" = 4
                                             ),
                                             selected = 2,
                                             inline = TRUE
                          )
                        ),
                        shiny::conditionalPanel(
                          condition = "input.redf == 2",
                          shinyWidgets::prettyRadioButtons("full2", "\\(\\text{Full Model}\\)",
                                             choices = list(
                                               "\\(\\text{Two-factor - main effect}\\)" = 3
                                             ),
                                             selected = 3,
                                             inline = TRUE
                          )
                        ),
                        shiny::conditionalPanel(
                          condition = "input.redf == 3",
                          shinyWidgets::prettyRadioButtons("full3", "\\(\\text{Full Model}\\)",
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
             shiny::conditionalPanel(
               condition = "input.ANOREG == 1&(input.Modef ==1|input.Modef ==2)",
               shiny::fluidRow(
                 shiny::column(width = 6,
                        shiny::sliderInput("f1", "\\(\\text{Factor 1 level:}\\)", min = 2, max = 10, value = 2, step = 1, ticks = FALSE)
                 ),
                 shiny::column(width = 6,
                        shiny::sliderInput("f2", "\\(\\text{Factor 2 level:}\\)", min = 2, max = 10, value = 2, step = 1, ticks = FALSE)
                 )
               )
             ),

             # Hypotheses
             shiny::fluidRow(
               shiny::column(width = 6,
                      shinyWidgets::prettyRadioButtons(
                        inputId = "h0f",
                        label = shiny::em("$$\\mathcal{H}_0:$$"),
                        choices = list(
                          "\\(\\lambda^2 = 0 \\)" = 1,
                          "\\(\\lambda^2 \\in \\{0, \\epsilon\\}\\)" = 2
                        ),
                        inline = TRUE,
                        selected = 1
                      )
               ),
               shiny::column(width = 6,
                      shiny::conditionalPanel(
                        condition = "input.h0f == 1",
                        shinyWidgets::prettyRadioButtons(
                          inputId = "h1f",
                          label = shiny::em("\\(\\mathcal{H}_1:\\)"),
                          choices = list("\\(\\lambda^2 > 0 \\)" = 1),
                          inline = TRUE,
                          selected = 1
                        )
                      ),
                      shiny::conditionalPanel(
                        condition = "input.h0f == 2",
                        shinyWidgets::prettyRadioButtons(
                          inputId = "h1fe",
                          label = shiny::em("\\(\\mathcal{H}_1:\\)"),
                          choices = list("\\(\\lambda^2 > \\epsilon\\)" = 1),
                          inline = TRUE,
                          selected = 1
                        ),
                        shiny::sliderInput("epsilinff", "\\(\\epsilon :\\)", min = 0.01, max = 0.25, value = 0.1, step = 0.01, ticks = FALSE)
                      )
               )
             ),

             # Prior
             shinyWidgets::prettyRadioButtons(
               inputId = "modelf",
               label = shiny::em("\\(\\text{Analysis Prior Distribution}\\)"),
               choices = list(
                 "\\( \\text{Effect size prior} \\)" = 1,
                 "\\( \\text{Moment prior (must df ≥ 3)} \\)" = 2
               ),
               inline = TRUE,
               selected = 1
             ),

             shiny::fluidRow(
               shiny::column(width = 4,
                      shiny::conditionalPanel(condition = "input.modelf == 1",
                                       shiny::sliderInput("rf", "\\(\\text{r scale:}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE)
                      )
               ),
               shiny::column(width = 4,
                      shiny::sliderInput("fsdf", "\\(\\mathcal{f}^2 :\\)", min = 0.01, max = 0.5, value = 0.1, step = 0.01, ticks = FALSE)
               ),
               shiny::column(width = 4,
                      shiny::sliderInput("dff", "\\(\\text{df :}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
               )
             ),

             # Design prior
             shiny::conditionalPanel(
               condition = "input.Modef == 1|input.Modef == 2",
               shinyWidgets::prettyRadioButtons(
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
             shiny::conditionalPanel(
               condition = "input.priorf == 2",
               shinyWidgets::prettyRadioButtons(
                 inputId = "modelfd",
                 label = shiny::em("\\(\\text{Design Prior Distribution}\\)"),
                 choices = list(
                   "\\( \\text{Effect size prior} \\)" = 1,
                   "\\( \\text{Moment prior (must df ≥ 3)} \\)" = 2,
                   "\\( \\text{Point} \\)" = 3
                 ),
                 inline = TRUE,
                 selected = 1
               ),
               shiny::fluidRow(
                 shiny::column(width = 4,
                        shiny::conditionalPanel(condition = "input.modelfd == 1",
                                         shiny::sliderInput("rfd", "\\(\\text{r scale:}\\)", min = 0, max = 3, value = 1, step = 0.01, ticks = FALSE)
                        ),
                        shiny::conditionalPanel(condition = "input.modelfd == 3",
                                         shiny::sliderInput("lfd", "\\(\\lambda^2:\\)", min = 0, max = 0.5, value = 0.1, step = 0.01, ticks = FALSE)
                        )
                 ),
                 shiny::column(width = 4,
                        shiny::conditionalPanel(condition = "input.modelfd == 1|input.modelfd == 2",
                                         shiny::sliderInput("fsdfd", "\\(\\mathcal{f}^2 :\\)", min = 0.01, max = 0.5, value = 0.1, step = 0.01, ticks = FALSE)
                        )
                 ),
                 shiny::column(width = 4,
                        shiny::conditionalPanel(condition = "input.modelfd == 1|input.modelfd == 2",
                                         shiny::sliderInput("dffd", "\\(\\text{df :}\\)", min = 1, max = 100, value = 1, step = 1, ticks = FALSE)
                        )
                 )
               )
             ),

             # Sample size determination controls
             shiny::conditionalPanel(condition = "input.Modef == 1",
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                       shiny::sliderInput("powerf", "\\( \\text{True Positive Evidence:} \\)", min = 0.5, max = 0.99, value = 0.8, step = 0.01, ticks = FALSE)
                                ),
                                shiny::column(width = 6,
                                       shiny::sliderInput("alphaf", "\\( \\text{False Positive Evidence:} \\)", min = 0.001, max = 0.05, value = 0.05, step = 0.001, ticks = FALSE)
                                )
                              )
             ),

             shiny::conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",
                              shiny::sliderInput("bff", "\\( \\text{Bound of compelling evidence:} \\)", min = 1, max = 20, value = 3, ticks = FALSE)
             ),
             shiny::conditionalPanel(condition = "input.Modef == 2",
                              shiny::numericInput("nf", "\\( \\text{Sample Size } N: \\)", value = 50)
             ),
             shiny::conditionalPanel(condition = "input.Modef == 1|input.Modef == 2",
                              shiny::actionButton("runf", label = "\\( \\text{Run} \\)"),
                              shiny::conditionalPanel(condition = "input.Modef == 1",
                                               shiny::em(shiny::span("\\(\\text{Note: Potential Error when the required N > 5,000} \\)", style = "color: red;"))
                              )
             ),

             # BF calculator mode
             shiny::conditionalPanel(condition = "input.Modef == 3",
                              shiny::fluidRow(
                                shiny::column(width = 4, shiny::numericInput("df1f", label = "\\( \\mathcal{df}_1: \\)", value = 1)),
                                shiny::column(width = 4, shiny::numericInput("df2f", label = "\\( \\mathcal{df}_2: \\)", value = 30)),
                                shiny::column(width = 4, shiny::numericInput("fval", label = "\\( f\\text{-value:} \\)", value = 1))
                              ),
                              shiny::actionButton("calf", label = "\\( \\text{Calculate} \\)"),
                              shiny::htmlOutput("BFcalf")
             ),

             shiny::em("\\( \\text{Recommended hyperparameters:} \\)"),
             shiny::htmlOutput("prior_suggest"),
             shiny::conditionalPanel(
               condition = "input.Modef == 1",
               shiny::checkboxGroupInput(
                 "o_plot_f",
                 label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
                 choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
                 selected = NULL
               ),
               shiny::downloadButton("export_f", "Download result as PDF")
             )
           ),

           # Main Panel
           shiny::mainPanel(shiny::fluidRow(
             shiny::column(6, shiny::plotOutput("priorff")),
             shiny::column(6, shiny::htmlOutput("resultf"))

           ),shiny::uiOutput("Optional_Plots_f")
           )
         )
)
,


shiny::navbarMenu(
  "\\(\\text{Proportion}\\)",

  shiny::tabPanel("\\(\\text{One proportion - binomial}\\)", shiny::withMathJax(),
           shiny::sidebarLayout(
             shiny::sidebarPanel(
               shinyWidgets::prettyRadioButtons(
                 "Modebin", "\\(\\text{Select Mode}\\)",
                 choices = list(
                   "\\(\\text{Sample size determination}\\)" = 1,
                   "\\(\\text{Fixed N}\\)" = 2,
                   "\\(\\text{BF calculator}\\)" = 3
                 ),
                 selected = 1, inline = TRUE
               ),

               shiny::fluidRow(
                 shiny::column(6,
                        shinyWidgets::prettyRadioButtons(
                          inputId = "h0bin", label = shiny::em("\\(\\mathcal{H}_0:\\)"),
                          choices = list(
                            "\\(p = p_0\\)" = 1,
                            "\\(p \\in (p_0 - \\epsilon, p_0 + \\epsilon)\\)" = 2
                          ),
                          selected = 1, inline = TRUE
                        )),

                 shiny::column(6,
                        shiny::conditionalPanel(
                          condition = "input.h0bin == 1",
                          shinyWidgets::prettyRadioButtons(
                            inputId = "h1bin", label = shiny::em("\\(\\mathcal{H}_1:\\)"),
                            choices = list(
                              "\\(p ≠ p_0\\)" = 1,
                              "\\(p > p_0\\)" = 2,
                              "\\(p < p_0\\)" = 3
                            ),
                            selected = 1, inline = TRUE
                          )
                        ),

                        shiny::conditionalPanel(
                          condition = "input.h0bin == 2",
                          shinyWidgets::prettyRadioButtons(
                            inputId = "h1bine", label = shiny::em("\\(\\mathcal{H}_1:\\)"),
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

               shiny::fluidRow(
                 shiny::column(4, shiny::sliderInput("h0prop", "\\(p_0\\)", min = .01, max = .99, value = .5, step = .01, ticks = FALSE)),

                 shiny::column(4,
                        shiny::conditionalPanel("input.h0bin == 2 && (input.h1bine == 1 || input.h1bine == 3)",
                                         shiny::sliderInput("lbbine", "\\(-\\epsilon\\)", min = -0.5, max = -0.01, value = -0.2, step = 0.01, ticks = FALSE),
                                         shiny::htmlOutput("bin_lower"))
                 ),

                 shiny::column(4,
                        shiny::conditionalPanel("input.h0bin == 2 && (input.h1bine == 1 || input.h1bine == 2)",
                                         shiny::sliderInput("ubbine", "\\(\\epsilon\\)", min = 0.01, max = 0.5, value = 0.2, step = 0.01, ticks = FALSE),
                                         shiny::htmlOutput("bin_upper"))
                 )
               ),

               shinyWidgets::prettyRadioButtons(
                 inputId = "modelbin",
                 label = shiny::em("\\(\\text{ Analysis Prior Distribution}\\)"),
                 choices = list("\\( \\text{Beta} \\)" = 1, "\\( \\text{Moment} \\)" = 2),
                 selected = 1, inline = TRUE
               ),

               shiny::fluidRow(
                 shiny::column(4,
                        shiny::conditionalPanel("input.modelbin == 1",
                                         shiny::sliderInput("alphabin", "\\(\\alpha\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)),
                        shiny::conditionalPanel("input.modelbin == 2",
                                         shiny::sliderInput("sbin", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE))
                 ),

                 shiny::column(4,
                        shiny::conditionalPanel("input.modelbin == 1",
                                         shiny::sliderInput("betabin", "\\(\\beta\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE))
                 )
               ),

               shiny::conditionalPanel("input.Modebin == 1 || input.Modebin == 2",
                                shinyWidgets::prettyRadioButtons(
                                  "priorbin", "\\( \\text{Design prior is the same as analysis prior: } \\)",
                                  choices = list("\\( \\text{Yes} \\)" = 1, "\\( \\text{No} \\)" = 2),
                                  selected = 1, inline = TRUE
                                ),

                                shiny::conditionalPanel("input.priorbin == 2",
                                                 shinyWidgets::prettyRadioButtons(
                                                   "modelbind", shiny::em("\\( \\text{Design prior distribution} \\)"),
                                                   choices = list(
                                                     "\\( \\text{Beta} \\)" = 1,
                                                     "\\( \\text{Moment} \\)" = 2,
                                                     "\\( \\text{Point} \\)" = 3
                                                   ),
                                                   selected = 1, inline = TRUE
                                                 ),

                                                 shiny::conditionalPanel("input.modelbind == 3",
                                                                  shiny::sliderInput("h0bind", "\\(p_1\\)", min = .01, max = .99, value = .5, step = .01, ticks = FALSE)),

                                                 shiny::fluidRow(
                                                   shiny::column(4,
                                                          shiny::conditionalPanel("input.modelbind == 1",
                                                                           shiny::sliderInput("alphabind", "\\(\\alpha\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)),

                                                          shiny::conditionalPanel("input.modelbind == 2",
                                                                           shiny::sliderInput("sbind", "\\(\\text{Scale}\\)", min = 0.01, max = 3, value = 1, step = 0.01, ticks = FALSE))
                                                   ),
                                                   shiny::column(4,
                                                          shiny::conditionalPanel("input.modelbind == 1",
                                                                           shiny::sliderInput("betabind", "\\(\\beta\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE))
                                                   )
                                                 )
                                )
               ),

               shiny::conditionalPanel(
                 "input.modelbind == 3 && (input.h1bin == 2 || input.h1bin == 3 || input.h1bine == 2 || input.h1bine == 3)",
                 shiny::em(shiny::span("\\(\\text{Note: The value should match the direction of } \\mathcal{H}_1\\)", style = "color: gray;"))
               ),

               shiny::conditionalPanel("input.Modebin == 1",
                                shiny::em("\\(\\text{Controlling for the probability of}\\)"),
                                shiny::fluidRow(
                                  shiny::column(6, shiny::sliderInput("powerbin", "\\(\\text{True Positive Evidence:}\\)", min = .5, max = .99, value = .8, step = .01, ticks = FALSE)),
                                  shiny::column(6, shiny::sliderInput("FP_bin", "\\(\\text{False Positive Evidence:}\\)", min = .001, max = .05, value = .05, step = .001, ticks = FALSE))
                                )
               ),

               shiny::conditionalPanel("input.Modebin == 1 || input.Modebin == 2",
                                shiny::sliderInput("bbin", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)),

               shiny::conditionalPanel("input.Modebin == 2 || input.Modebin == 3",
                                shiny::fluidRow(
                                  shiny::column(6, shiny::numericInput("nbin", "\\(\\text{Sample Size:}\\)", value = 50)),
                                  shiny::column(6,
                                         shiny::conditionalPanel("input.Modebin == 3",
                                                          shiny::numericInput("xbin", "\\(\\text{Number of Success:}\\)", value = 25)))
                                )
               ),

               shiny::conditionalPanel("input.Modebin == 1 || input.Modebin == 2",
                                shiny::actionButton("runbin", label = "\\(\\text{Run}\\)"),
                                shiny::conditionalPanel("input.Modebin == 1",
                                                 shiny::em(shiny::span("\\(\\text{Note: Error when the required N > 10,000}\\)", style = "color: red;")))
               ),

               shiny::conditionalPanel("input.Modebin == 3",
                                shiny::actionButton("calbin", label = "\\(\\text{Calculate}\\)"),
                                shiny::htmlOutput("BFbin")),
               shiny::conditionalPanel(
                 condition = "input.Modebin == 1",
                 shiny::checkboxGroupInput(
                   "o_plot_bin",
                   label = "\\(\\text{Additional Plots (computationally intensive):}\\)",
                   choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
                   selected = NULL
                 ),
                 shiny::downloadButton("export_bin", "Download result as PDF")
               )


             ),

             shiny::mainPanel(shiny::fluidRow(
               shiny::column(6, shiny::plotOutput("prior_bin")),
               shiny::column(6, shiny::htmlOutput("resultbin"))
             ),
             shiny::uiOutput("Optional_Plots_bin"))
           )
  ), shiny::tabPanel("\\(\\text{Two proportion}\\)",shiny::withMathJax(),
shiny::sidebarLayout(shiny::sidebarPanel(
  shinyWidgets::prettyRadioButtons(
    "Modep2", "\\(\\text{Select Mode}\\)",
    choices = list(
      "\\(\\text{Sample size determination}\\)" = 1,
      "\\(\\text{Fixed N}\\)" = 2,
      "\\(\\text{BF calculator}\\)" = 3
    ),
    selected = 1, inline = TRUE
  ) ,shiny::fluidRow(
    shiny::column(5,
           shiny::em("\\(\\mathcal{H}_0: p_0 = p_1 = p_2\\)"),
           shiny::br(),
           shiny::em("\\(p_0 \\sim \\text{Beta}(\\alpha_0, \\beta_0)\\)")




           ),




    shiny::column(5,
           shiny::em("\\(\\mathcal{H}_1: p_1 \\neq p_2\\)"),
           shiny::br(),
           shiny::em("\\(p_1 \\sim \\text{Beta}(\\alpha_1, \\beta_1)\\)"),
           shiny::br(),
           shiny::em("\\(p_2 \\sim \\text{Beta}(\\alpha_2, \\beta_2)\\)")
    )
  ),shiny::fluidRow(

  shiny::column(4,
  shiny::sliderInput("alpha0", "\\(\\alpha_0\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
  shiny::sliderInput("beta0", "\\(\\beta_0\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)


           ),
  shiny::column(4,
         shiny::sliderInput("alpha1", "\\(\\alpha_1\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
         shiny::sliderInput("beta1", "\\(\\beta_1\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)
           ),
  shiny::column(4,
         shiny::sliderInput("alpha2", "\\(\\alpha_2\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
         shiny::sliderInput("beta2", "\\(\\beta_2\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)

  )

  ),
  shiny::conditionalPanel("input.Modep2 == 1 | input.Modep2 == 2",
  shinyWidgets::prettyRadioButtons(
    "priorp2", "\\( \\text{Design prior is the same as analysis prior: } \\)",
    choices = list("\\( \\text{Yes} \\)" = 1, "\\( \\text{No} \\)" = 2),
    selected = 1, inline = TRUE
  )),

  shiny::conditionalPanel(
    condition = "input.priorp2 == 2",
  shiny::fluidRow(
    shiny::column(4,
           shinyWidgets::prettyRadioButtons(
             "model_p1",
             label = "\\( \\text{Model for } p_1: \\)",
             choices = list(
               "\\( \\text{Fixed } p_1 \\)" = 1,
               "\\( p_1 \\sim \\text{Beta}(\\alpha_{1d}, \\beta_{1d}) \\)" = 2
             ),
             selected = 1,
             inline = TRUE
           ),
           shinyWidgets::prettyRadioButtons(
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

  shiny::column(4,
         shiny::conditionalPanel(
           condition = "input.model_p1 == 1",
           shiny::sliderInput("location1d", "\\(p_1 =\\)", min = 0.01, max = 0.99, value = 0.5, step = 0.01, ticks = FALSE)
         ),
         shiny::conditionalPanel(
           condition = "input.model_p1 == 2",
           shiny::sliderInput("alpha1d", "\\(\\alpha_{1d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
           shiny::sliderInput("beta1d", "\\(\\beta_{1d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)
         )
  )
  ,
  shiny::column(4,
         shiny::conditionalPanel(
           condition = "input.model_p2 == 1",
           shiny::sliderInput("location2d", "\\(p_2 =\\)", min = 0.01, max = 0.99, value = 0.5, step = 0.01, ticks = FALSE)
         ),
         shiny::conditionalPanel(
           condition = "input.model_p2 == 2",
           shiny::sliderInput("alpha2d", "\\(\\alpha_{2d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE),
           shiny::sliderInput("beta2d", "\\(\\beta_{2d}\\)", min = 0.01, max = 100, value = 1, step = 1, ticks = FALSE)
         )
  )



  )),

  shiny::conditionalPanel("input.Modep2 == 1",
                   shiny::em("\\(\\text{Controlling for the probability of}\\)"),
                   shiny::fluidRow(
                     shiny::column(12, shiny::sliderInput("powerp2", "\\(\\text{True Positive Evidence:}\\)", min = .5, max = .99, value = .8, step = .01, ticks = FALSE)),
                     #shiny::column(6, shiny::sliderInput("FP_p2", "\\(\\text{False Positive Evidence:}\\)", min = .001, max = .05, value = .05, step = .001, ticks = FALSE))
                   )
  ),
  shiny::conditionalPanel("input.Modep2 == 1 || input.Modep2 == 2",
                   shiny::sliderInput("bp2", "\\(\\text{Bound of compelling evidence:}\\)", min = 1, max = 20, value = 3, ticks = FALSE)),

  shiny::conditionalPanel("input.Modep2 == 2 ||input.Modep2 == 3 ",
                   shiny::em("\\( \\text{Sample size per group} \\)"),
                   shiny::fluidRow(
                     shiny::column(6, shiny::numericInput("n1p2", "\\(n_1\\)", value = 50)),
                     shiny::column(6, shiny::numericInput("n2p2", "\\(n_2\\)", value = 50))
                   )

                   ),
  shiny::conditionalPanel("input.Modep2 == 3",
                   shiny::em("\\( \\text{Number of success} \\)"),
                   shiny::fluidRow(
                     shiny::column(6, shiny::numericInput("x1p2", "\\(x_1\\)", value = 50)),
                     shiny::column(6, shiny::numericInput("x2p2", "\\(x_2\\)", value = 50))
                   )) ,

  shiny::conditionalPanel("input.Modep2 == 1 || input.Modep2 == 2",
                   shiny::actionButton("runp2", label = "\\(\\text{Run}\\)"),
                   shiny::conditionalPanel("input.Modep2 == 1", shiny::em(shiny::span("\\(\\text{Note: Error when required } N > 5,000\\)", style = "color: red;")))
  ),
  shiny::conditionalPanel("input.Modep2 == 3",
                   shiny::actionButton("calp2", label = "\\(\\text{Calculate}\\)"),
                   shiny::htmlOutput("BFp2")
  ),

  shiny::conditionalPanel(
    condition = "input.Modep2 == 1",
    shiny::checkboxGroupInput(
      "o_plot_p2",
      label = "\\(\\text{Additional Plots (computationally very intensive):}\\)",
      choices = list("\\(\\text{Power Curve}\\)" = 1, "\\(\\text{Relationship between BF and data}\\)" = 2),
      selected = NULL
    ),
    shiny::downloadButton("export_p2", "Download result as PDF")
  )
),shiny::mainPanel(

  shiny::fluidRow(
    shiny::column(6, shiny::plotOutput("prior_p0")),
    shiny::column(6, shiny::htmlOutput("resultp2"))
  ),
  shiny::fluidRow(
    shiny::column(6, shiny::plotOutput("prior_p1")),
    shiny::column(6, shiny::plotOutput("prior_p2"))
  ),
  shiny::uiOutput("Optional_Plots_p2")

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

#' @export
BayesPower_testing <- function(){
  # Run the application
 shiny::shinyApp(ui = ui, server = server)

}
