launch_app_addin <- function() {
  choice <- rstudioapi::showQuestion(
    "Launch options",
    "Where would you like to open the Shiny app?",
    ok = "Browser",
    cancel = "RStudio Viewer"
  )
  options(shiny.launch.browser = choice)
  callr::r_bg(
    func = function() BayesPower::BayesPower_BayesFactor(),
    package = TRUE
  )
}
