#' @export
server_bin<- function(input, output, session) {
input_bin <- reactive({
  mode_bf <- switch(input$Modebin,
                    "1" = 1,
                    "2" = 0,
                    "3" = 0)# mode



  interval <- input$h0bin # point null or interval

  hypothesis <- switch(interval,
                       "1" =   switch(input$h1bin,        # direction of the test
                                      "1" = "!=",
                                      "2" =  ">",
                                      "3" =  "<"),
                       "2" = switch(input$h1bine,        # direction of the test
                                    "1" = "!=",
                                    "2" =  ">",
                                    "3" =  "<"))
  location <- input$h0prop
  lbbin <- input$lbbine
  ubbin <- input$ubbine

  if ((location+lbbin)<(0)){
    lbbin = lbbin+-1-(location+lbbin)

  }

  if ((location+ubbin)>(+1)){
    ubbin = ubbin+1-(location+ubbin)

  }


  e <- switch(input$h1bine,        # bound for interval test
              "1" = c(lbbin, ubbin),
              "2" = ubbin,
              "3" = lbbin)

  inter <- switch(interval,
                  "1" = input$h1bin,
                  "2" = input$h1bine)


  model <- switch(input$modelbin,
                  "1" = "beta",
                  "2" = "Moment")
  alpha <- input$alphabin
  beta <- input$betabin
  scale <- input$sbin
  de_an_prior <- switch(input$priorbin,
                        "1" = 1,
                        "2" = 0)
  alpha_d <- input$alphabind
  beta_d <- input$betabind
  scale_d <- input$sbind
  model_d <- switch(input$modelbind,
                  "1" = "beta",
                  "2" = "Moment",
                  "3" = "Point")
  location_d <- input$h0bind
  target <- input$powerbin
  FP <- input$FP_bin
  D <- input$bbin
  N <- input$nbin
  Suc <- input$xbin
  pc   <- "1" %in% input$o_plot_bin
  rela <- "2" %in% input$o_plot_bin

  if (mode_bf!=1){
    pc=rela=F
  }
 ############


  # Add all variables to the final list
  list(
    mode_bf = mode_bf,
    interval = interval,
    hypothesis =hypothesis,
    location = location,
    e = e,
    lbbin = lbbin,
    ubbin = ubbin,
    inter = inter,
    model = model,
    alpha = alpha,
    beta = beta,
    scale = scale,
    de_an_prior = de_an_prior,
    alpha_d = alpha_d,
    beta_d = beta_d,
    scale_d = scale_d,
    location_d = location_d,
    model_d = model_d,
    target = target,
    FP = FP,
    D = D,
    N = N,
    Suc =Suc,
    pc=pc,
    rela=rela

  )
})



output$bin_lower<-renderUI({
  bin = input_bin()


  table_html <-  paste0('
                        p_0 - \\epsilon = ', bin$location+bin$lbbin,'')

  tagList(
    # Render the table using MathJax
    withMathJax(
      em('$$', table_html, '$$')
    )
  )

})


output$bin_upper<-renderUI({
 bin = input_bin()

  table_html <-  paste0('
                        \\rho_0 - \\epsilon = ', bin$location+bin$ubbin,'')

  tagList(
    # Render the table using MathJax
    withMathJax(
      em('$$', table_html, '$$')
    )
  )


})




observeEvent(input$runbin, {
  bin = input_bin()

  dat <- tryCatch({switch(bin$interval,
                "1" = {bin_table(bin$D,bin$target,bin$alpha,bin$beta,bin$location,
                  bin$scale,bin$model,bin$hypothesis,
                  bin$alpha_d,bin$beta_d,bin$location_d,bin$scale_d,
                  bin$model_d,bin$de_an_prior,bin$N, bin$mode_bf,bin$FP)},
                "2" = {
                  bin_e_table(bin$D,bin$target,bin$alpha,bin$beta,bin$location,
                                 bin$scale,bin$model,bin$hypothesis,
                                 bin$alpha_d,bin$beta_d,bin$location_d,bin$scale_d,
                                 bin$model_d,bin$de_an_prior,bin$N, bin$mode_bf,bin$FP,bin$e)
                  })}, error = function(e) {
                    "Error"
                  })

  output$prior_bin <- renderPlot({
    switch(bin$interval,
           "1" = {bin_prior_plot(bin$alpha,bin$beta,bin$location,bin$scale,bin$model,
                                 bin$alpha_d,bin$beta_d,bin$location_d,
                                 bin$scale_d,bin$model_d,bin$hypothesis,
                                 bin$de_an_prior)},
           "2" =bin_e_prior_plot (bin$alpha,bin$beta,bin$location,bin$scale,
                                  bin$model,bin$alpha_d,bin$beta_d,bin$location_d,
                                  bin$scale_d,bin$model_d,
                                  bin$hypothesis,bin$de_an_prior,bin$e))



  })

  output$resultbin <- renderUI({
    if (identical(dat, "Error")){
      table_html <- "Required sample size is more than the upper limit."
    }else{
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('$$', '
    \\begin{array}{l c}
    \\textbf{Probability of Compelling Evidence} & \\\\
    \\hline
    \\text{p(BF}_{10} > ', bin$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[1], 3), ' \\\\
    \\text{p(BF}_{01} > ', bin$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[3], 3), ' \\\\
    \\textbf{Probability of Misleading Evidence} & \\\\
    \\hline
    \\text{p(BF}_{01} > ', bin$D, '\\, | \\, \\mathcal{H}_1)\\ & ', round(dat[2], 3), ' \\\\
    \\text{p(BF}_{10} > ', bin$D, '\\, | \\, \\mathcal{H}_0)\\ & ', round(dat[4], 3), ' \\\\
    \\textbf{Required Sample Size} & \\\\
    \\hline
    \\text{N} & ', dat[5], ' \\\\
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

  if (bin$pc) {
    switch(bin$interval,
           "1" = Power_bin(bin$D,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis,
                           bin$alpha_d,bin$beta_d,bin$location_d,
                           bin$scale_d,bin$model_d, bin$de_an_prior,dat[1,5]),
           "2" = Power_e_bin(bin$D,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis,
                             bin$alpha_d,bin$beta_d,bin$location_d,
                             bin$scale_d,bin$model_d, bin$de_an_prior,dat[1,5],bin$e))

    pc_bin <- recordPlot()
  } else pc_bin <- NA


  if (bin$rela) {
    switch(bin$interval,
           "1" =bin_bf10(bin$D,dat[1,5],bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis),
           "2" =bin_e_bf10(bin$D,dat[1,5],bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis,bin$e) )

    rela_bin <- recordPlot()
  } else rela_bin <- NA


  output$Optional_Plots_bin <- renderUI({
    tagList(
      if (bin$pc) {
        tagList(
          withMathJax(em("$$\\text{Power Curve}$$")),
          output$PCbin <- renderPlot({

            pc_bin

          })
        )
      },
      if (bin$rela) {
        tagList(
          withMathJax(em("$$\\text{Relationship between BF and data}$$")),
          output$bfrbin <- renderPlot({

            rela_bin

          })
        )
      }
    )
  })



  output$export_bin <- downloadHandler(
    filename = function() {
      "BayesPower-report.pdf"
    },
    content = function(file) {
      template_path <- system.file("report_templates", "report_bin.Rmd", package = "BayesPower")

      tempReport <- file.path(tempdir(), "report_bin.Rmd")
      file.copy(template_path, tempReport, overwrite = TRUE)

      rmarkdown::render(
        input = tempReport,
        output_file = file,
        params = list(bin = bin, dat = dat,pc_bin=pc_bin,rela_bin=rela_bin),  # ✅ pass to `params`
        envir = new.env(parent = globalenv())  # environment still required
      )
    }
  )




})


observeEvent(input$calbin, {
  bin = input_bin()
  BF10 <- switch(bin$interval,
                 "1" = bin_BF(bin$Suc,bin$N,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis),
                 "2" = bin_e_BF(bin$Suc,bin$N,bin$alpha,bin$beta,bin$location,bin$scale,bin$model,bin$hypothesis,bin$e))

  output$BFbin <- renderUI({
    # Create the LaTeX formatted strings for the table
    table_html <- paste0('
    N = ', bin$N, ', x = ', bin$Suc, '; \\textit{BF}_{10} = ', round(BF10, 4), '
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
