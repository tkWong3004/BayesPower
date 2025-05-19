input_t1 <- reactive({
  mode_bf <- switch(input$Modet1,
                    "1" = 1,
                    "2" = 0)# mode
  N <-  switch(input$Modet1,
               "1" = 2,
               "2" = input$nt1,
               "3" = input$t1df)
               
  
  interval <- input$h0t1 # point null or interval
  
  e <- switch(input$h1t1e,        # bound for interval test
              "1" = c(input$lbt1e, input$ubt1e),
              "2" = input$ubt1e,
              "3" = input$lbt1e)
  inter <- switch(interval,       
                  "1" = input$h1t1,
                  "2" = input$h1t1e)
  
  hypothesis <- switch(interval,
                       "1" =   switch(input$h1t1,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1t1e,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))
  
  model <- switch(input$modelt1,
                  "1" = "t-distribution",
                  "2" = "Normal",
                  "3" = "NLP")
  
  location <- switch(input$h0t1,
                     "1" =input$lt1,
                     "2" = 0)
  scale <- input$st1
  dff <- input$dft1
  de_an_prior <- switch(input$prior,
                        "1" = 1,
                        "2" = 0)
  model_d <- switch(input$modelt1d,
                    "1" = "t-distribution",
                    "2" = "Normal",
                    "3" = "NLP",
                    "4" = "Point")
  location_d <- input$lt1d
  
  scale_d <- input$st1d
  dff_d <- input$dft1d
  D <- input$bt1
  type <- input$typet1
  target <- input$powert1
  alpha <- input$alphat1
 
  tval <- input$t1tval
  
  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    interval = interval,
    hypothesis = hypothesis ,
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
    alpha = alpha ,
    N = N,
    tval = tval
  )
})

observeEvent(input$runt1, {
  x = input_t1()
  
  dat = suppressWarnings(switch(x$interval, "1" =  t1_Table(x$D,x$target,x$model,x$location,x$scale,x$dff, x$hypothesis,
                  x$model_d,x$location_d,x$scale_d,x$dff_d, x$de_an_prior,x$N, x$mode_bf ,
                  x$alpha),"2" = t1e_table(x$D,x$target,x$model,x$scale,x$dff, x$hypothesis,x$e ,
                                           x$model_d,x$scale_d,x$dff_d, x$de_an_prior,x$N,x$mode_bf,x$location_d ,x$alpha)))
  
  output$priort1 <- renderPlot({
    suppressWarnings(switch(x$interval,
           "1"= t1_prior_plot(
      D = x$D,                  # Access 'D' explicitly
      target = x$target,        # Access 'target' explicitly
      model = x$model,          # Access 'model' explicitly
      location = x$location,    # Access 'location' explicitly
      scale = x$scale,          # Access 'scale' explicitly
      dff = x$dff,              # Access 'dff' explicitly
      hypothesis = x$hypothesis,  # Access 'hypothesis' explicitly
      model_d = x$model_d,        # Access 'model_d' explicitly
      location_d = x$location_d,  # Access 'location_d' explicitly
      scale_d = x$scale_d,        # Access 'scale_d' explicitly
      dff_d = x$dff_d,            # Access 'dff_d' explicitly
      de_an_prior = x$de_an_prior   # Access 'de_an_prior' explicitly
    ),
    "2" = t1e_prior_plot(x$model,
                   x$scale,
                   x$dff ,
                   x$hypothesis,
                   x$e,
                   x$de_an_prior,
                   x$model_d,
                   x$scale_d,
                   x$dff_d,
                   x$location )
    
  ))
    
  })
  
  output$bfrt1 <- renderPlot({
    
    suppressWarnings(switch(x$interval,
           "1"=
    bf10_t1(
      D = x$D,                  # Access 'D' explicitly
      df = dat[1,5],               # Access 'dff' (degrees of freedom) for 'df'
      target = x$target,        # Access 'target' explicitly
      model = x$model,          # Access 'model' explicitly
      location = x$location,    # Access 'location' explicitly
      scale = x$scale,          # Access 'scale' explicitly
      dff = x$dff,              # Access 'dff' explicitly again, if required
      hypothesis = x$hypothesis # Access 'hypothesis' explicitly
    ), "2"= te1_BF (x$D,dat[1,5],x$model ,x$scale,x$dff , x$hypothesis ,x$e)))
    
  })
  
  output$resultt1 <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    \\text{p(BF}_{10} > ', x$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[1], 3), ' \\\\
    \\text{p(BF}_{01} > ', x$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[3], 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    \\text{p(BF}_{01} > ', x$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[2], 3), ' \\\\
    \\text{p(BF}_{10} > ', x$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[4], 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
    \\end{array}
  ')
    
    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em('$$', table_html, '$$')
      )
    )
  })
  
  
  output$PCt1 <- renderPlot({
    suppressWarnings(switch(x$interval,
           "1"=
    Power_t1(x$D,x$model,x$location,x$scale,x$dff, x$hypothesis,
                       x$model_d,x$location_d,x$scale_d,x$dff_d, x$de_an_prior,dat[1,5]),
    "2" = Power_t1e(x$D,x$model,x$location,x$scale,x$dff, x$hypothesis,
                    x$model_d,x$location_d,x$scale_d,x$dff_d, x$de_an_prior,dat[1,5],x$e)))
  })
  
})

observeEvent(input$cal1, {
  x = input_t1()
  BF10 <- switch(x$interval,
                 "1" = t1_BF10(x$tval,x$N,x$model ,x$location,x$scale,x$dff , x$hypothesis ),
                 "2" = t1e_BF10(x$tval,x$N,x$model,x$scale,x$dff , x$hypothesis,x$e ))
  
  output$BFt1 <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{t}(', x$N , ') = ',x$tval,', \\textit{BF}_{10} = ', round(BF10, 4), '
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




