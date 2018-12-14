CompareProportionGroups <- function(n1, n2, f) {
  n <- cbind(n1, n2); # Count matrix of sucesses and fails
  f <- as.factor(f);  # Group factor
  
  fit <- glm(n ~ f, family=binomial);
  anova(fit, test='Chisq');
} 