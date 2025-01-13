# Plot PCA
helper.plot.PCA <- function(df, show.legend = FALSE) {
   pl <- ggplot2::ggplot(
      data = df,
      ggplot2::aes(
         x = .data[[names(df)[1]]],  # First column
         y = .data[[names(df)[2]]],  # Second column
         color = .data[[names(df)[3]]]  # Third column
      )
   ) +
      ggplot2::geom_point(show.legend = show.legend) +
      ggplot2::theme_bw()
   
   return(pl)
}
