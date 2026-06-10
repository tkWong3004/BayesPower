# This function creates an Addin that launches the Shiny app from RStudio's
# Addins menu on the top (under the menus bar).

launch_app_addin <- function() {
  # Launch the app in a background Rscript process on a fixed port, without
  # opening a browser:
  port <- httpuv::randomPort()
  system2(
    "Rscript",
    args = c("-e", shQuote(paste0(
      "options(shiny.port = ", port, ", shiny.launch.browser = FALSE);",
      "BayesPower::BayesPower_BayesFactor()"
    ))),
    wait = FALSE
  )

  # Wait until the loading the app above is ready:
  repeat {
    conn <- tryCatch(
      suppressWarnings(socketConnection("127.0.0.1", port, timeout = 1)),
      error = function(e) NULL
    )
    if (!is.null(conn)) { close(conn); break }
    Sys.sleep(0.2)
  }

  # Launch the app in the browser:
  utils::browseURL(paste0("http://127.0.0.1:", port))
}


