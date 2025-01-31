fit.ReferenceLine <- function(x, y) {
   # Compute residuals for y relative to the reference line
   residuals <- y - x
   # Compute R-squared
   ss_residual <- sum(residuals^2)  # Residual sum of squares
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
   # Calculate mean of y
   y_mean <- mean(y)
   # Extract residual
   residuals <- y - y_mean
   # Compute sum of squares (RSS)
   rss <- sum(residuals^2)
   
   # Compute p-value for residuals (testing if mean residual = 0)
   p_value <- t.test(residuals, mu = 0)$p.value
   if(is.na(p_value)){p_value <- 1}
   # Return results
   ret <- list("1-rss" = rss,
               "p.value" = p_value)
   return(ret)
}



