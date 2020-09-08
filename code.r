#The function lc() in the open source of R
rm(list=ls()) # Clears the memory of R
library(MASS) # Load the package useful for the determination of lethal concentrations
lc<-function(d,x,t){ # This function has three main entries that are the concentrations of effluents (d), the number of dead (x) after   treatment and the total number of species (t) that received effluent concentrations (lc)
  options(warn=-1) # This is to avoid warming alerts
  # The following function will compute the corrected mortality rates
  Abbott <- function(d,x,t){
    ab<-(1-((t-x)/(t[d==0]-x[d==0])))
    ifelse(ab < 0,0,ab[])}
  Abb<-Abbott(d,x,t)
  mynewdat<-data.frame(d,x,t,Abb); (mynewdat = mynewdat[!d==0,]) # This part removed the control which is no more useful for
  #the determination of l.c
  mod<-list()
  #The following help to select the binomial family link
  family<-function (i){
    if (i==1){family=binomial(link = "probit")}
    else{
      if (i==2){family=binomial(link = "logit")}
      else{
        if (i==3){family=binomial(link = "cloglog")}
      }}
    #return(family)
  }
  for (i in 1:3){
    mod[[i]]<-glm(Abb ~ log(d),family=family(i),data=mynewdat)
  }
  b<-which.min(c(deviance(mod[[1]]),deviance(mod[[2]]),deviance(mod[[3]]))) # Compare and select the best link
  cat("The best model is the model", b,"\n") # Gives the output of the previous comparison
  out=list()
  out$Resum<-summary(mod[[b]])
  ld<-dose.p(mod[[b]],p=c(0.50,0.90,0.95)) # Computation of lethal concentrations 50, 90 and 95 but may also 10, 20, 80 etc.
  ld.ci <- as.vector(ld)+ attr(ld, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1) # Defines the confidence intervals
  out$ld.est <- round(exp((cbind(ld, attr(ld, "SE"), ld.ci[,1], ld.ci[,2]))),3)
  dimnames(out$ld.est)[[2]] <- c("LD", "SE", "LCL","UCL") # Return the lethal concentrations and their confidence intervals
  #The following part plots the curve
  x1=seq(min(d),max(d), 0.001)
  lpredmod <- predict(mod[[b]], data.frame(d=x1), type="response")
  plot(d,Abb, pch=16, ylim=c(0,1),xlab="Concentration ", ylab = "Mortality rate")
  lines(x1,lpredmod)
  #Goodness of fit of the model
  out$comparison<-pchisq(deviance(mod[[b]]), mod[[b]]$df.residual, lower=FALSE) # Comparison of deviance
  out$R2_Naglekerke<-round((1-exp((mod[[b]]$dev-mod[[b]]$null)/sum(t)))/(1-exp(-mod[[b]]$null/sum(t))),3) #Determination of
 # Naglekerke R square
  return(out)
}
