################# This function is written by Koissi Savi and Kossi Doulabe





rm(list=ls()) # Cleaning the environment of R
# Abbott correction function
Abbott <- function(d, x, t){ 
  ab <- 1 - ((t - x)/(t[d == 0]-x[d == 0]))
  ab <- ifelse(ab < 0, yes = 0, no = ab)
  ab
} 

# function for regression model with family spec
reg <- function(Abb, dose, link){
  instruc <- paste0(
    "glm(Abb ~ ", 
    paste(paste0("dose, family = stats::quasibinomial(link =", link,")"), collapse = ""),
    ")"
  )
  eval(parse(text = instruc))
}


# The function lc, main function embedding the determination of lethal concentration, the CI and the #goodness of fit of the selected model used for the determination of the lethal concentrations
lc <- function(d, x, t, prob = c(0.50,0.90,0.95)){ 
  ab <- Abbott(d, x, t) 
  mynewdat <- data.frame(d, x, t, ab)
  mynewdat <- mynewdat[!d == 0, ] 
  # This part removed the control which is no more useful for the determination of l.c 
  mod <- vector(mode = 'list', length = 3) 
  links <- c("probit", "logit", "cloglog")
  for(l in seq_along(links)){
    if( all(d < 1) ){
      mod[[l]] <- reg(Abb = mynewdat$ab, dose = mynewdat$d, link = links[l]) #to avoid negative values of doses in dose.p
    } else {
      mod[[l]] <- reg(Abb = mynewdat$ab, dose = log(mynewdat$d), link = links[l])
    } 
  }
  b <- which.min(sapply(mod, FUN = deviance))
  Resum <- summary( mod[[b]] ) 
  ld <- MASS::dose.p(mod[[b]], p = prob) # Computation of lethal concentrations 50, 90 and 95 but may also 10, 20, 80 etc. 
  ld.ci <- ld + attr(ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow = 1) # Defines the confidence intervals 
  ld.est <- round(cbind(exp(as.vector(ld)), as.vector(attr(ld, "SE")), exp(ld.ci[,1]), exp(ld.ci[,2])), 3)
  dimnames(ld.est)[[2]] <- c("LD", "SE", "LCL","UCL") # Return the lethal concentrations and their  confidence intervals 
  #Goodness of fit of the model 
  x1=seq(min(d),max(d), 0.001)
  lpredmod <- predict(mod[[b]], data.frame(d=x1), type="response")
  plot(d,ab, pch=16, ylim=c(0,1),xlab="Concentration ", ylab = "Mortality rate")
  #lines(x1,lpredmod)
  comparison <- pchisq(deviance(mod[[b]]), mod[[b]]$df.residual, lower.tail = FALSE) 
  # Comparison of deviance
  R2_Naglekerke <- round((1 - exp((mod[[b]]$dev - mod[[b]]$null)/sum(t)))/(1 - exp(- mod[[b]]$null/sum(t))), 3) 
  #Determination of  Naglekerke R square 
  out <- list(Resum = Resum, ld.est = ld.est, comparison = comparison, R2_Naglekerke = R2_Naglekerke) 
  out 
  
} 



lethal <- lc(d , x , t , prob = c(0.50,0.90,0.95, 0.99))
lethal
