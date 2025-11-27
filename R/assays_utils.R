fit.ReferenceLine <- function(x, y) {
   # Remove NA values
   valid_idx <- !is.na(x) & !is.na(y)
   x <- x[valid_idx]
   y <- y[valid_idx]
   
   # Check for sufficient data (at least 2 points required)
   if (length(x) < 2 || length(y) < 2) {
      warning("Not enough points for computing fit.ReferenceLine")
      return(list("r.squared" = NA, "p.value" = NA))
   }
   # Compute residuals for y relative to the reference line
   residuals <- y - x
   # Compute R-squared
   ## Root Mean Squared Error (RMSE)
   ss_residual <- sqrt(mean(residuals^2))  
   # Use ANOVA to compare the full model to the null hypothesis (y = x)
   # Full model: y ~ x (allowing slope and intercept to vary)
   model_full <- lm(y ~ x)
   #r_squared <- summary(model_full)$r.squared
   # Null model: force y = x (slope = 1 and intercept = 0)
   model_null <- lm(y ~ offset(x))
   # Perform ANOVA comparison
   anova_result <- anova(model_null, model_full)
   p_value_anova <- anova_result$`Pr(>F)`[2]  # p-value for comparison
   
   if(is.na(p_value_anova)){p_value_anova <- 1}
   
   # Return results
   ret <- list(
      "r.squared" = ss_residual,         # R-squared value
      "p.value" = p_value_anova  # p-value from ANOVA
   )
   return(ret)
}



fit.Constant <- function(x, y) {
   # remove NA
   y <- y[!is.na(y)]
   
   # if fewer than 2 values, return NA
   if (length(y) < 2) {
      warning("Not enough points for computing fit.Constant")
      return(list("1-rss" = NA,
                  "p.value" = NA))
   }
   
   # Calculate mean of y
   y_mean <- mean(y)
   # Extract residual
   residuals <- y - y_mean
   # Compute Root Mean Squared Error (RMSE)
   rss <- sqrt(mean(residuals^2))
   
   # Compute p-value for residuals (testing if mean residual = 0)
   p_value <- t.test(residuals, mu = 0)$p.value
   if(is.na(p_value)){p_value <- 1}
   # Return results
   ret <- list("1-rss" = rss,
               "p.value" = p_value)
   return(ret)
}

monotonicity_score <- function(x) {
   if (length(x) < 2) return(1)
   diffs <- diff(x)
   sum(diffs > 0) / (length(x) - 1)
}

consistency_drop <- function(x, # consistency score
                             y # it has to be 1 (no split) or 0.5 (split)
                             ){
   whole <- x[y == 1]
   split <- x[y == 0.5]
   
   diff <- whole - split
   score <- 1 - diff
   score <- ifelse(score<0, 0, score)
   return(score)
}



