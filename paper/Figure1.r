
# This R script provides the codes for reproducing the Figure 1.

library(BayesPower)
library(patchwork)
# Input
threshold      <- 3
n              <- 100
df             <- n-1
prior_analysis <- "t-distribution"
location       <- 0
scale          <- .707
dff            <- 1
alternative    <- "two.sided"
prior_design   <- "Point"
location_d     <- 0.2


# Panel 1
tt <- seq(-5, 5, 0.2)

# Compute BF10 and bounds
BF10   <- BayesPower:::t1_BF10(tt, df, prior_analysis, location, scale, dff, alternative)
t.BF10 <- BayesPower:::t1_BF10_bound(threshold, df, prior_analysis, location, scale, dff, alternative)


# BF10 title
main.bf10 <- if (length(t.BF10) == 1) {
  bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when t = " ~ .(round(t.BF10, 2))))
} else {
  bquote(bold("BF"[10] ~ "=" ~ .(threshold) ~ " when t = " ~ .(round(t.BF10[1], 2)) ~
                " or " ~ .(round(t.BF10[2], 2))))
}


df_bf10 <- data.frame(tt = tt, BF = BF10)

clean_theme <- ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 14, face = "bold"),
    axis.text  = ggplot2::element_text(size = 12),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
  )

## ---------- BF10 ----------
x_breaks_10 <- sort(unique(c(-5, 5, round(t.BF10, 2))))

p1 <- ggplot2::ggplot(df_bf10, ggplot2::aes(tt, BF)) +
  ggplot2::geom_line(linewidth = 1.2, color = "black") +
  ggplot2::geom_vline(xintercept = t.BF10, linetype = "dashed") +
  ggplot2::scale_y_log10() +
  ggplot2::scale_x_continuous(limits = c(-5, 5), breaks = x_breaks_10) +
  ggplot2::labs(
    x = "t-value",
    y = expression("BF"[10] * " (log scale)"),
    title = main.bf10
  ) +
  clean_theme



# Panel 2

# Use a finer grid and explicitly include the BF cutoffs
tt2 <- sort(unique(c(
  seq(-5, 5, length.out = 1000),
  t.BF10
)))

h1_design <- dt(
  tt2,
  df,
  ncp = sqrt(n) * .2
)

df_h1 <- data.frame(
  tt = tt2,
  density = h1_design
)

lower <- min(t.BF10)
upper <- max(t.BF10)

p2 <- ggplot2::ggplot(
  df_h1,
  ggplot2::aes(x = tt, y = density)
) +

  # Left tail
  ggplot2::geom_ribbon(
    data = subset(df_h1, tt <= lower),
    ggplot2::aes(
      ymin = 0,
      ymax = density
    ),
    fill = "lightcoral",
    alpha = 0.4
  ) +

  # Right tail
  ggplot2::geom_ribbon(
    data = subset(df_h1, tt >= upper),
    ggplot2::aes(
      ymin = 0,
      ymax = density
    ),
    fill = "lightcoral",
    alpha = 0.4
  ) +

  # Density curve
  ggplot2::geom_line(
    linewidth = 1.2,
    color = "black"
  ) +

  # BF boundaries
  ggplot2::geom_vline(
    xintercept = t.BF10,
    linetype = "dashed"
  ) +

  ggplot2::scale_x_continuous(
    limits = c(-5, 5),
    breaks = x_breaks_10
  ) +

  ggplot2::labs(
    x = "t-value",
    y = "Density",
    title = expression(bold(delta == 0.2))
  ) +

  clean_theme
# Panel 3

# df range
dfs <- seq(2, ceiling(336 * 1.2), length.out = 31)

# Power curve
TPE <- sapply(dfs, function(df_i) {
  t10 <- BayesPower:::t1_BF10_bound(
    threshold, df_i,
    prior_analysis, location, scale, dff, alternative
  )

  BayesPower:::t1_TPE(
    t10, df_i,
    prior_design, location_d,
    scale = 0, dff = 0, alternative
  )
})

df1 <- data.frame(
  SampleSize = dfs + 1,
  Probability = TPE
)

# Helper function: power for a given total sample size
get_power <- function(N) {

  df_i <- N - 1

  t10 <- BayesPower:::t1_BF10_bound(
    threshold, df_i,
    prior_analysis, location, scale, dff, alternative
  )

  BayesPower:::t1_TPE(
    t10, df_i,
    prior_design, location_d,
    scale = 0, dff = 0, alternative
  )
}

# Power at n and 336
power_n   <- get_power(n)
power_336 <- get_power(336)


p3 <- ggplot2::ggplot(
  df1,
  ggplot2::aes(x = SampleSize, y = Probability)
) +

  ggplot2::geom_line(
    linewidth = 1.2,
    color = "black"
  ) +

  # Vertical dotted lines
  ggplot2::geom_segment(
    x = n, xend = n,
    y = 0, yend = power_n,
    linetype = "dotted",
    linewidth = 0.8
  ) +

  ggplot2::geom_segment(
    x = 336, xend = 336,
    y = 0, yend = power_336,
    linetype = "dotted",
    linewidth = 0.8
  ) +

  # Points
  ggplot2::geom_point(
    x = n, y = power_n,
    size = 2.5
  ) +

  ggplot2::geom_point(
    x = 336, y = power_336,
    size = 2.5
  ) +

  # Horizontal arrow from n to 336
  ggplot2::geom_segment(
    x = n,
    y = power_n,
    xend = 336,
    yend = power_n,
    arrow = grid::arrow(
      length = grid::unit(0.25, "cm"),
      type = "closed"
    ),
    linewidth = 0.8
  ) +

  # Horizontal reference line at .80
  ggplot2::geom_hline(
    yintercept = 0.80,
    linetype = "dotted",
    linewidth = 0.8
  ) +

  ggplot2::scale_x_continuous(
    breaks = sort(unique(c(
      pretty(df1$SampleSize),
      n,
      336
    )))
  ) +

  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = sort(unique(c(
      seq(0, 1, by = 0.2),
      0.80
    )))
  ) +

  ggplot2::labs(
    x = "Total sample size",
    y = "Probability",
    title = bquote(
      bold("Power curve for BF"[10] ~ ">" ~ .(threshold))
    )
  ) +

  clean_theme


# conbime plot

combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3)

combined_plot

ggplot2::ggsave(
  filename = "Procedure.pdf",
  plot = combined_plot,
  width = 15,
  height = 5,
  dpi = 300
)

