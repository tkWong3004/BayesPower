input_f <- reactive({
  mode_bf <- switch(input$Modef,
                    "1" = 1,
                    "2" = 0)

  if (input$ANOREG == 2){
    p <- input$pf
    k <- input$kf
  }else{
    p <-switch(input$redf,
               "1" = 1,
               "2" = input$f1-1+1,
               "3" = input$f1-1 +input$f2-1 +1
    )
    full_model <-switch(input$redf,
                        "1"=input$full1,
                        "2"=input$full2,
                        "3"=input$full3)
    k <-switch(full_model,
               "2" = input$f1-1+1,
               "3" = input$f1-1 +input$f2-1 +1,
               "4" = input$f1-1 +input$f2-1 +1 + (input$f1-1)*(input$f2-1)
    )
  }

  

  inter <- input$h0f
  
  e <- switch(inter,
              "1" = input$epsilinff,
              "2" = input$epsilinff)

  
  model <- switch(input$modelf,
                  "1" = "tdis",
                  "2" = "Moment")
  
  rscale <- input$rf
  f_m <- sqrt(input$fsdf)
  dff <- input$dff
  de_an_prior <- switch(input$priorf,
                        "1" = 1,
                        "2" = 0)
  
  model_d <- switch(input$modelfd,
                    "1" = "tdis",
                    "2" = "Moment",
                    "3" = "Point")
  
  rscale_d <- input$rfd
  f_m_d <- sqrt(input$fsdfd)
  
  if ( input$modelfd == "3"){
    f_m_d <-sqrt(input$lfd)
  }
  
  dff_d <- input$dffd
  target <- input$powerf
  alpha <- input$alphaf
  N <- input$nf
  D <- input$bff
  fval <- input$fval
  df1 <- input$df1f
  df2 <- input$df2f
  q = k -p
  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    p = p,
    k = k,
    q = q,
    inter=inter,
    e=e,
    model=model,
    rscale=rscale,
    f_m=f_m,
    dff=dff,
    de_an_prior=de_an_prior,
    model_d=model_d,
    rscale_d=rscale_d,
    f_m_d=f_m_d,
    dff_d=dff_d,
    target=target,
    alpha=alpha,
    N=N,
    D=D,
    fval = fval,
    df1=df1,
    df2=df2
  )
})


output$prior_suggest <- renderUI({
  ff = input_f()
  if (ff$model == "tdis"){
  
  table_html <- paste0('
\\textit{df} = ', 3, ', \\textit{r} = \\sqrt{\\frac{df - 2}{dfq}} \\times f = ',round(sqrt((3 - 2) / 3*ff$q) * sqrt(ff$f_m),2),'
')
  }else{
  
    table_html <- paste0('
\\textit{df = 5+(q-1)} = ', 5+ff$q-1,'
')
  
  
}
  # Render the table using MathJax
  tagList(
    # Render the table using MathJax
    withMathJax(
      em('$$', table_html, '$$')
    )
  )
})


observeEvent(input$runf, {
  ff = input_f()
  
  dat = switch(ff$inter,
               "1" = f_table(ff$D,ff$target,ff$p,ff$k,ff$dff,ff$rscale,ff$f_m,ff$model,
                ff$dff_d,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,ff$N, ff$mode_bf,ff$alpha ),
               "2" = fe_table(ff$D,ff$target,ff$p,ff$k,ff$dff,ff$rscale,ff$f_m,ff$model,
                              ff$dff_d,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,ff$N, ff$mode_bf,ff$e ,ff$alpha))
  output$priorff <- renderPlot({
    
    switch(ff$inter,
           "1" =prior_plot_f(ff$q,ff$dff,ff$rscale,ff$f_m,ff$model,ff$dff_d
                 ,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior),
           "2" = prior_plot_fe(ff$q,ff$dff,ff$rscale,ff$f_m,ff$model,ff$dff_d
                               ,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,ff$e))
    
  })
  
  output$resultf <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    \\text{p(BF}_{10} > ', ff$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[1], 3), ' \\\\
    \\text{p(BF}_{01} > ', ff$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[3], 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    \\text{p(BF}_{01} > ', ff$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[2], 3), ' \\\\
    \\text{p(BF}_{10} > ', ff$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[4], 3), ' \\\\
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

  
 
  output$bff <- renderPlot({
    
    switch(ff$inter,
           "1" = bf10_f(ff$D,dat[1,5],ff$k,ff$p,ff$dff,ff$rscale,ff$f_m,ff$model),
           "2" = bf10_fe(ff$D,dat[1,5],ff$k,ff$p,ff$dff,ff$rscale,ff$f_m,ff$model,ff$e))
    
  })
  
  
  output$PCf <- renderPlot({
    switch(ff$inter,
    "1" = Power_f(ff$D,ff$k,ff$p,ff$dff,ff$rscale,ff$f_m,ff$model,ff$k_d,ff$p_d,
            ff$dff_d,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,dat[1,5]),
    "2" = Power_fe(ff$D,ff$k,ff$p,ff$dff,ff$rscale,
                   ff$f_m,ff$model,ff$k_d,ff$p_d,ff$dff_d,ff$rscale_d,ff$f_m_d,ff$model_d,ff$de_an_prior,dat[1,5],ff$e))
    
  })
  
})






observeEvent(input$calf, {
  ff = input_f()
  m = ff$df2+ff$df1
  BF10 <- F_BF(ff$fval,ff$df1,m,ff$dff,ff$rscale,ff$f_m,ff$model)
  
  output$BFcalf <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{F(',ff$df1,',',ff$df2,')} = ',ff$fval,', \\textit{BF}_{10} = ', round(BF10, 4), '
')
    
    
    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em(table_html)
      )
    )
  })
  
})