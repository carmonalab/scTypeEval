# Plot PCA
helper.plot.scatter <- function(df,
                                show.legend = FALSE,
                                label = TRUE) {
   pl <- ggplot2::ggplot() +
      ggplot2::geom_point(
         data = df,
         ggplot2::aes(
            x = .data[[names(df)[1]]],  # First column
            y = .data[[names(df)[2]]],  # Second column
            color = .data[[names(df)[3]]]  # Third column
         ),
         show.legend = show.legend) +
      ggplot2::theme_bw()
   
   if(label){
      # Calculate centroids
      cnt <- df |>
         dplyr::group_by(.data[[names(df)[3]]]) |>
         dplyr::summarize(
            x = mean(.data[[names(df)[1]]]),
            y = mean(.data[[names(df)[2]]])
         )
      
      # Add labels to the plot
      pl <- pl +
         ggrepel::geom_label_repel(
            data = cnt,
            ggplot2::aes(
               x = x,
               y = y,
               label = .data[[names(df)[3]]], # Use group name as label
               color = .data[[names(df)[3]]]  # Match colors
            ),
            alpha = 0.8,
            show.legend = FALSE # Avoid duplicate legends for labels
         )
   }
   
   return(pl)
}
