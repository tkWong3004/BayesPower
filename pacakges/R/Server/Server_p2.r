input_p2 <- reactive({


  mode_bf <- switch(input$Modep2,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode

  a0 <- input$alpha0
  b0 <- input$beta0

  a1 <- input$alpha1
  b1 <- input$beta1

  a2 <- input$alpha2
  b2 <- input$beta2

  a1d <- input$alpha1d
  b1d <- input$beta1d

  a2d <- input$alpha2d
  b2d <- input$beta2d

  dp1 <- input$location1d
  dp2 <- input$location2d


  model_p1 <-switch(input$model_p1,
                    "1" = "Point",
                    "2" = "beta")
  model_p2 <-switch(input$model_p2,
                    "1" = "Point",
                    "2" = "beta")

  if (input$priorp2 == 1){
    model_p1 = model_p2 = "same"
  }
  de_an_prior<-input$priorp2
  D <- input$bp2

  n1 <- input$n1p2
  n2 <- input$n2p2

  x1 <- input$x1p2
  x2 <- input$x2p2

  target <- input$powerp2
  pc   <- "1" %in% input$o_plot_p2
  rela <- "2" %in% input$o_plot_p2

  if (mode_bf!=1){
    pc=rela=F
  }
 ############

  list(
    mode_bf = mode_bf,
    a0 = a0,
    b0 = b0,
    a1 = a1,
    b1 = b1,
    a2 = a2,
    b2 = b2,
    a1d =  a1d,
    b1d = b1d,
    a2d = a2d,
    b2d = b2d,
    model1 = model_p1,
    model2 = model_p2,
    dp1 = dp1,
    dp2 = dp2,
    D = D,
    n1 = round(n1),
    n2 = round(n2),
    k1 = round(x1),
    k2 = round(x2),
    target = target,
    r=1,
    pc=pc,
    rela=rela,
    de_an_prior=de_an_prior

  )
})



observeEvent(input$runp2, {
  p2  <- input_p2()
  dat <- tryCatch({pro_table_p2(p2$D,p2$target, p2$a0, p2$b0,
                      p2$a1, p2$b1, p2$a2, p2$b2, p2$r,
                      p2$model1,p2$a1d,p2$b1d,p2$dp1,
                      p2$model2,p2$a2d,p2$b2d,p2$dp2,
                      p2$mode_bf,p2$n1,p2$n2)},
    error = function(e) {
                        "Error"
                      })
  if (any(dat != "Error")){
    table <- dat[[1]]
    grid  <- dat[[2]]

  }



  output$prior_p0 <- renderPlot({
    p2_prior_plot(p2$a0,p2$b0,1,1,0,"same",0)
  })

  output$prior_p1 <- renderPlot({
    p2_prior_plot(p2$a1,p2$b1,p2$a1d,p2$b1d,p2$dp1,p2$model1,1)
  })
  output$prior_p2 <- renderPlot({
    p2_prior_plot(p2$a2,p2$b2,p2$a2d,p2$b2d,p2$dp2,p2$model2,2)
  })

  output$resultp2 <- renderUI({
    if (identical(dat, "Error")){
      table_html <- "Required sample size is more than the upper limit."
    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('$$', '
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    \\text{p(BF}_{10} > ', p2$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(table[1,1], 3), ' \\\\
    \\text{p(BF}_{01} > ', p2$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(table[1,3], 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    \\text{p(BF}_{01} > ', p2$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(table[1,2], 3), ' \\\\
    \\text{p(BF}_{10} > ', p2$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(table[1,4], 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N}_1 & ', table[1,5], ' \\\\
    \\text{N}_2 & ', table[1,6], ' \\\\
    \\end{array}
  ', '$$')
    }
    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em(table_html)
      )
    )
  })
  # Save plots based on conditions

  if (p2$pc) {
    Power_p2(p2$D, table[1,5], p2$a0, p2$b0, p2$a1, p2$b1, p2$a2,
             p2$b2, table[1,6] / table[1,5], p2$model1, p2$a1d, p2$b1d, p2$dp1,
             p2$model2, p2$a2d, p2$b2d, p2$dp2)
    pc_p2 <- recordPlot()
  } else pc_p2 <- NA

  # Render plots OUTSIDE renderUI
  output$PCp2 <- renderPlot({
    pc_p2
  })



  if (p2$rela) {

    rela_p2 <-heatmap_p2(grid, p2$D)
  } else rela_p2 <- NA


  output$bfrp2 <- renderPlot({
    print(rela_p2)
  })

  # Render the UI with proper plotOutput
  output$Optional_Plots_p2 <- renderUI({
    tagList(
      if (p2$pc) {
        tagList(
          withMathJax(em("$$\\text{Power Curve}$$")),
          plotOutput("PCp2")
        )
      },
      if (p2$rela) {
        tagList(
          withMathJax(em("$$\\text{Relationship between BF and data}$$")),
          plotOutput("bfrp2")
        )
      }
    )
  })

  output$export_p2 <- downloadHandler(
    filename = function() {
      "BayesPower-report.pdf"
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report_2p.Rmd")
      file.copy("R/pdf/report_2p.Rmd", tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,
        output_file = file,
        params = list(p2 = p2, dat = dat,pc_p2=pc_p2,rela_p2=rela_p2),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )

})


observeEvent(input$calp2, {
  p2 = input_p2()
  BF10 <- BF10_p2(p2$a0, p2$b0, p2$a1, p2$b1, p2$a2, p2$b2,p2$n1,p2$n2,p2$k1,p2$k2)

  output$BFp2 <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0(
      'n_1 = ', p2$n1, ', ',
      'n_2 = ', p2$n2, ', ',
      'x_1 = ', p2$k1, ', ',
      'x_2 = ', p2$k2, ', ',
      '\\textit{BF}_{10} = ', round(BF10, 4)
    )





    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em('$$', table_html, '$$')
      )
    )
  })

})
