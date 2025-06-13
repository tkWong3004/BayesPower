input_t2 <- reactive({
  mode_bf <- switch(input$Modet2,
                    "1" = 1,
                    "2" = 0,
                    "3" = 3)# mode
  interval <- input$h0t2 # point null or interval

  e <- switch(input$h1t2e,        # bound for interval test
              "1" = c(input$lbt2e, input$ubt2e),
              "2" = input$ubt2e,
              "3" = input$lbt2e)
  inter <- switch(interval,
                  "1" = input$h1t2,
                  "2" = input$h1t2e)

  hypothesis <- switch(interval,
                       "1" =   switch(input$h1t2,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1t2e,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))



  model <- switch(input$modelt2,
                  "1" = "t-distribution",
                  "2" = "Normal",
                  "3" = "NLP")

  location <- switch(input$h0t2,
                     "1" =input$lt2,
                     "2" = 0,
                     "3" = 0,
                     "4" = 0)
  scale <- input$st2
  dff <- input$dft2
  de_an_prior <- switch(input$priort2,
                        "1" = 1,
                        "2" = 0)
  model_d <- switch(input$modelt2d,
                    "1" = "t-distribution",
                    "2" = "Normal",
                    "3" = "NLP",
                    "4" = "Point")
  location_d <- switch(interval,
                       "1" = input$lt2d,
                       "2" = 0)
  scale_d <- input$st2d
  dff_d <- input$dft2d
  D <- input$bt2
  type <- input$typet2
  target <- input$powert2
  alpha <- input$alphat2
  tval <- input$t2tval
  r <- switch(input$Modet2,
              "1" = input$rt2,
              "2" = input$n2t2/input$n1t2,
              "3" = input$rt2)
  N1 = input$n1t2
  N2 = input$n2t2
  df = input$t2df
  pc   <- "1" %in% input$o_plot_t2
  rela <- "2" %in% input$o_plot_t2
  if (mode_bf!=1){
    pc=rela=F
  }
  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    interval = interval,
    hypothesis = hypothesis,
    e = e,
    model = model,
    location = location,
    scale = scale,
    dff = dff,
    de_an_prior = de_an_prior,
    model_d = model_d,
    location_d = location_d,
    scale_d = scale_d,
    dff_d = dff_d,
    type = type,
    D = D,
    target = target,
    alpha = alpha,
    tval = tval,
    r = r,
    N1=N1,
    N2=N2,
    df = df,
    pc = pc,
    rela = rela
  )
})

observeEvent(input$runt2, {
  t2 = input_t2()
  dat <- tryCatch({
    suppressWarnings(switch(t2$interval,
                            "1" = t2_Table(t2$D, t2$r, t2$target, t2$model, t2$location, t2$scale, t2$dff, t2$hypothesis,
                                           t2$model_d, t2$location_d, t2$scale_d, t2$dff_d, t2$de_an_prior, t2$N1, t2$N2, t2$mode_bf, t2$alpha),
                            "2" = t2e_table(t2$D, t2$r, t2$target, t2$model, t2$scale, t2$dff, t2$hypothesis, t2$e,
                                            t2$model_d, t2$scale_d, t2$dff_d, t2$de_an_prior, t2$mode_bf, t2$location, t2$N1, t2$N2, t2$alpha)
    ))
  }, error = function(e) {
    "Error"
  })


  output$priort2 <- renderPlot({
    suppressWarnings(switch(t2$interval,
           "1"=
             t1_prior_plot(
               D = t2$D,                  # Access 'D' explicitly
               target = t2$target,        # Access 'target' explicitly
               model = t2$model,          # Access 'model' explicitly
               location = t2$location,    # Access 'location' explicitly
               scale = t2$scale,          # Access 'scale' explicitly
               dff = t2$dff,              # Access 'dff' explicitly
               hypothesis = t2$hypothesis,  # Access 'hypothesis' explicitly
               model_d = t2$model_d,        # Access 'model_d' explicitly
               location_d = t2$location_d,  # Access 'location_d' explicitly
               scale_d = t2$scale_d,        # Access 'scale_d' explicitly
               dff_d = t2$dff_d,            # Access 'dff_d' explicitly
               de_an_prior = t2$de_an_prior   # Access 'de_an_prior' explicitly
             ), "2" =
             t1e_prior_plot(t2$model,
                            t2$scale,
                            t2$dff ,
                            t2$hypothesis,
                            t2$e,
                            t2$de_an_prior,
                            t2$model_d,
                            t2$scale_d,
                            t2$dff_d,
                            t2$location )

    ))

  })



  output$resultt2 <- renderUI({
    # Create the LaTeX formatted strings for the table
    if (identical(dat, "Error")){
      table_html <- "Required sample size is more than the upper limit."
    }else{
    table_html <- paste0("$$",'
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    \\text{p(BF}_{10} > ', t2$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[1,1], 3), ' \\\\
    \\text{p(BF}_{01} > ', t2$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[1,3], 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    \\text{p(BF}_{01} > ', t2$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[1,2], 3), ' \\\\
    \\text{p(BF}_{10} > ', t2$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[1,4], 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N}_1 & ', dat[1,5], ' \\\\
    \\text{N}_2 & ', dat[1,6], ' \\\\
    \\end{array}
  ',"$$")}

    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em(table_html)
      )
    )
  })


  if (t2$pc) {
    suppressWarnings(switch(t2$interval,
                            "1" = Power_t2(t2$D,t2$model,t2$location,t2$scale,t2$dff, t2$hypothesis,
                                           t2$model_d,t2$location_d,t2$scale_d,t2$dff_d, t2$de_an_prior,dat[1,5],dat[1,6]/dat[1,5]),
                            "2" = Power_t2e(t2$D,t2$model,t2$location,t2$scale,t2$dff, t2$hypothesis,
                                            t2$model_d,t2$location_d,t2$scale_d,t2$dff_d, t2$de_an_prior,dat[1,5],dat[1,6]/dat[1,5])))
    pc_t2 <- recordPlot()
  } else pc_t2 <- NA

  if (t2$rela) {
    suppressWarnings(switch(t2$interval,
                            "1"= t2_bf10(t2$D ,dat[1,5],t2$r, t2$target,t2$model ,t2$location ,t2$scale,t2$dff  , t2$hypothesis ),
                            "2" =t2e_BF (t2$D,dat[1,5],t2$r,t2$model ,t2$scale,t2$dff , t2$hypothesis ,t2$e) ))
    rela_t2 <- recordPlot()
  } else rela_t2 <- NA

  output$Optional_Plots_t2 <- renderUI({
    tagList(
      if (t2$pc) {
        tagList(
          withMathJax(em("$$\\text{Power Curve}$$")),
          output$PCt2 <- renderPlot({

            pc_t2

          })
        )
      },
      if (t2$rela) {
        tagList(
          withMathJax(em("$$\\text{Relationship between BF and data}$$")),
          output$bfrt2 <- renderPlot({

            rela_t2

          })
        )
      }
    )
  })


  output$export_t2 <- downloadHandler(
    filename = function() {
      "BayesPower-report.pdf"
    },
    content = function(file) {
      tempReport <- file.path(tempdir(), "report_t2.Rmd")
      file.copy("R/pdf/report_t2.Rmd", tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,
        output_file = file,
        params = list(t2 = t2, dat = dat,pc_t2=pc_t2,rela_t2=rela_t2),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )







})

observeEvent(input$cal1, {
  t2 = input_t2()

  N = t2$df +2
  N1 = N/(1+t2$r)
  BF10 =  suppressWarnings(t2_BF10(t2$tval,N1,t2$r,t2$model ,t2$location,t2$scale,t2$dff , t2$hypothesis ))

  output$BFt2 <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{t}(', t2$df, ') = ',t2$tval,', \\textit{BF}_{10} = ', round(BF10, 4), '
')


    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em('$$', table_html, '$$')
      )
    )
  })

})




