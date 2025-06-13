#' @export
server_r<- function(input, output, session) {
input_r <- reactive({
  mode_bf <- switch(input$Moder,
                    "1" = 1,
                    "2" = 2,
                    "3" = 3)# mode


  interval <- input$h0r # point null or interval
  h0 <- input$h0pho
  lbre <- input$lbre
  ubre <- input$ubre

  if ((h0+lbre)<(-1)){
    lbre = lbre+-1-(h0+lbre)

  }

  if ((h0+ubre)>(+1)){
    ubre = ubre+1-(h0+ubre)

  }



  e <- switch(input$h1re,        # bound for interval test
              "1" = c(lbre, ubre),
              "2" = ubre,
              "3" = lbre)
  inter <- switch(interval,
                  "1" = input$h1r,
                  "2" = input$h1re)


  hypothesis <- switch(interval,
                       "1" =   switch(input$h1r,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1re,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))





  model <- switch(input$modelr,
                  "1" = "d_beta",
                  "2" = "beta",
                  "3" = "NLP")
  k <- input$kr
  scale <- input$sr
  alpha <- input$ralpha
  beta <- input$rbeta
  de_an_prior <- switch(input$priorr,
                        "1" = 1,
                        "2" = 0)
  model_d <- switch(input$modelrd,
                    "1" = "d_beta",
                    "2" = "beta",
                    "3" = "NLP",
                    "4" = "Point")
  location_d <- input$h0phod
  k_d <- input$rkd
  scale_d <- input$rsd
  alpha_d <- input$ralphad
  beta_d<- input$rbetad
  target <- input$powerr
  FP <- input$alphapr
  D <- input$br
  N <-  switch(input$Moder,
               "1" = 2,
               "2" = input$nr,
               "3" = input$rdf+1)
  rval <- input$rval
  pc   <- 1 %in% input$o_plot_r
  rela <- 2 %in% input$o_plot_r
  if (mode_bf!=1){
    pc=rela=F
  }
  ###########
  location <- h0
  dff <- 1

  dff_d <- 1



 ############


  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    interval = interval,
    e = e,
    lbre = lbre,
    ubre = ubre,
    inter = inter,
    hypothesis = hypothesis,
    h0 = h0,
    model = model,
    k =k,
    scale = scale,
    alpha = alpha,
    beta = beta,
    de_an_prior = de_an_prior,
    model_d =model_d,
    location_d = location_d,
    k_d = k_d,
    scale_d = scale_d,
    alpha_d = alpha_d,
    beta_d = beta_d,
    target = target,
    FP = FP,
    D = D,
    N = N,
    rval = rval,
    location = location,
    dff = dff ,
    dff_d = dff_d,
    pc = pc,
    rela = rela

  )
})

output$r_lower<-renderUI({
 rr = input_r()


  table_html <-  paste0('
                        \\rho_0 - \\epsilon = ', rr$h0+rr$lbre,'')

  tagList(
    # Render the table using MathJax
    withMathJax(
      em('$$', table_html, '$$')
    )
  )

})


output$r_upper<-renderUI({
  rr = input_r()

  table_html <-  paste0('
                        \\rho_0 - \\epsilon = ', rr$h0+rr$ubre,'')

  tagList(
    # Render the table using MathJax
    withMathJax(
      em('$$', table_html, '$$')
    )
  )


})

observeEvent(input$runr, {
  rr = input_r()

  dat <- tryCatch({
    switch(rr$interval,

    "1" = r_table(rr$D,rr$target,rr$model,rr$k,
            rr$alpha, rr$beta,rr$h0,rr$location,
            rr$scale,rr$dff, rr$hypothesis ,rr$model_d,
            rr$location_d,rr$k_d, rr$alpha_d, rr$beta_d,
            rr$scale_d,rr$dff_d,rr$de_an_prior,rr$N,
            rr$mode_bf,rr$FP ),
    "2" = re_table(rr$D,rr$target,rr$model,rr$k,
                  rr$alpha, rr$beta,rr$h0,rr$location,
                  rr$scale,rr$dff, rr$hypothesis ,rr$model_d,
                  rr$location_d,rr$k_d, rr$alpha_d, rr$beta_d,
                  rr$scale_d,rr$dff_d,rr$de_an_prior,rr$N,
                  rr$mode_bf,rr$FP,rr$e ))
  }, error = function(e) {
    "Error"
  })

  output$prior_r <- renderPlot({

    switch(rr$interval,
           "1" = r_prior_plot(rr$k, rr$alpha, rr$beta,
                              rr$h0,rr$location,rr$scale,
                              rr$dff,rr$model,rr$de_an_prior,
                              rr$k_d, rr$alpha_d, rr$beta_d,
                              rr$location_d,rr$scale_d,rr$dff_d,
                              rr$model_d,rr$hypothesis),
           "2" = re_prior_plot(rr$k, rr$alpha, rr$beta,
                              rr$h0,rr$location,rr$scale,
                              rr$dff,rr$model,rr$de_an_prior,
                              rr$k_d, rr$alpha_d, rr$beta_d,
                              rr$location_d,rr$scale_d,rr$dff_d,
                              rr$model_d,rr$hypothesis,rr$e))



  })

  output$resultr <- renderUI({
    if (identical(dat, "Error")){
      table_html <- "Required sample size is more than the upper limit."
    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('$$','
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    \\text{p(BF}_{10} > ', rr$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[1], 3), ' \\\\
    \\text{p(BF}_{01} > ', rr$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[3], 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    \\text{p(BF}_{01} > ', rr$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[2], 3), ' \\\\
    \\text{p(BF}_{10} > ', rr$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[4], 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
    \\end{array}
  ','$$')
    }
    # Render the table using MathJax
    tagList(
      # Render the table using MathJax
      withMathJax(
        em(table_html)
      )
    )
  })




  if (rr$pc) {
    switch(rr$interval,
                     "1" = Power_r(rr$D,rr$k, rr$alpha, rr$beta,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model,
                                   rr$k_d, rr$alpha_d, rr$beta_d,rr$location_d,rr$scale_d,rr$dff_d,rr$model_d, rr$de_an_prior,dat[1,5]),
                     "2" = Power_re(rr$D,rr$k, rr$alpha, rr$beta,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model,
                                    rr$k_d, rr$alpha_d, rr$beta_d,rr$location_d,rr$scale_d,rr$dff_d,rr$model_d, rr$de_an_prior,dat[1,5],rr$e))

    pc_r <- recordPlot()
  }else{
    pc_r <-NA
    }

  if (rr$rela) {
    switch(rr$interval,
           "1" = r_bf10_p(rr$D,dat[1,5],rr$k,rr$alpha, rr$beta,rr$h0,
                          rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model),
           "2" = re_bf10_p(rr$D,dat[1,5],rr$k,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model,rr$e))


    rela_r <- recordPlot()
  } else{
    rela_r <- NA
  }



  output$Optional_Plots_r <- renderUI({
    tagList(
      if (rr$pc) {
        tagList(
          withMathJax(em("$$\\text{Power Curve}$$")),
          output$PCr <- renderPlot({

            pc_r

          })
        )
      },
      if (rr$rela) {
        tagList(
          withMathJax(em("$$\\text{Relationship between BF and data}$$")),
          output$bfrr <- renderPlot({

            rela_r

          })
        )
      }
    )
  })


  output$export_r <- downloadHandler(
    filename = function() {
      "BayesPower-report.pdf"
    },
    content = function(file) {
      template_path <- system.file("report_templates", "report_r.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_r.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,
        output_file = file,
        params = list(rr = rr, dat = dat,pc_r=pc_r,rela_r=rela_r),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )



})

observeEvent(input$calr, {
  rr = input_r()
  BF10 <- switch(rr$interval ,
                 "1" = r_BF10(rr$rval,rr$N,rr$k, rr$alpha, rr$beta,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model),
                 "2" = re_BF10(rr$rval,rr$N,rr$k, rr$alpha, rr$beta,rr$h0,rr$hypothesis,rr$location,rr$scale,rr$dff,rr$model,rr$e))

  output$BFrv <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    \\textit{r}(', rr$N-1 , ') = ',rr$rval,', \\textit{BF}_{10} = ', round(BF10, 4), '
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
}



