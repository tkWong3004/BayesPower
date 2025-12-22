setwd("~/MyGitHub/BayesPower/package/R")

# Questions:
# - Parameter e: Why ">" requires e>0 and "<" requires e<0?
# - Why dff = 0 when checking on prior_analysis?

#devtools::install_github(repo = "tkWong3004/BayesPower", subdir = "package")
devtools::document()
devtools::install()
library(BayesPower)

BFpower.ttest.OneSample(
  alternative = "two.sided", 
  e = NULL, 
  prior_analysis = "t", location = 0, scale = 0.707, dff = 10,
  prior_design = "normal", location_d = 0, scale_d = 0.707, dff_d = 1,
  type_rate = "positive", true_rate = 0.8,  false_rate = 0.05,
  D = 10, 
  plot_power = FALSE,  
  plot_rel = FALSE
)

BFpower.ttest.TwoSamples(
  alternative = "two.sided", 
  e = NULL, 
  prior_analysis = "t", location = 0, scale = 0.707, dff = 10,
  prior_design = "normal", location_d = 0, scale_d = 0.707, dff_d = 1,
  type_rate = "positive", true_rate = 0.8,  false_rate = 0.05,
  D = 10, 
  plot_power = FALSE,  
  plot_rel = FALSE, 
  r = 1
)

BFpower.ttest.TwoSamples(
  alternative = "greater", 
  e = .2, 
  prior_analysis = "t", location = 0, scale = 0.707, dff = 10,
  prior_design = "normal", location_d = 0, scale_d = 0.707, dff_d = 1,
  type_rate = "positive", true_rate = 0.8,  false_rate = 0.05,
  D = 10, 
  plot_power = FALSE,  
  plot_rel = FALSE, 
  N1             = 30,
  N2             = 25
)














