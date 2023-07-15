################################
### Gini coefficient ###########

#  Trapezoidal approcimation
auc_trapez_approx <- function(x, y){
  dy <- c(diff(y), 0)
  dx <- c(diff(x), 0)
  sum(y*dx) + sum(dx*dy)/2
}

#Obtain the lorenz curve (gini)
# WARNING: pdis, gdis = diff(physical_position), diff(genetic_position)
# Modified fucntion from Kaur and Rockman, 2014  https://doi.org/10.1534/genetics.113.158857 (credit them)

lorenz_curve <- function(pdis,gdis){
  #Sort the intervals by their recombination rate	
  rdata = cbind(pdis , gdis) 
  rdata=rdata[order(rdata[,2]/rdata[,1]),]
  
  #Rescale the intervals so that distances are rendered as proportions
  rdata[,1] = rdata[,1]/sum(rdata[,1])
  rdata[,2] = rdata[,2]/sum(rdata[,2])
  
  #Collect the cumulative distributions
  ginidata = cbind(cumsum(rdata[,1]) , cumsum(rdata[,2]))
  ginidata=as.data.frame(ginidata)
  names(ginidata) = c("cum.pos", "cum.gen")
  return(ginidata)
}


#Obtain gini index from lorenz curve
ginifromlorenz = function(lo){
  auc_approx <- function(x, y){
    dy <- c(diff(y), 0)
    dx <- c(diff(x), 0)
    sum(x*dy) + sum(dx*dy)/2
  }
  
  gini.coef = 1-2*auc_approx(lo$cum.gen,lo$cum.pos)
  names(gini.coef) <- "Gini.coefficient"
  return(gini.coef)
}


# Get gini coefficient from Marey maps
mapgini <- function(pdis,gdis){
  locurve=lorenz_curve(pdis,gdis)
  ginifromlorenz(locurve)
}








#################################
#################################
#Relative area under the curve

# x = physical distance; y = genetic disatnce 
rel_auc = function(x,y){
  rdata = cbind(x , y) 
  #Rescale the intervals so that distances are rendered as proportions
  rdata[,1] = rdata[,1]/max(rdata[,1])
  rdata[,2] = rdata[,2]/max(rdata[,2])
  
  #plot(rdata[,1],rdata[,2])
  
  relauc = auc_trapez_approx(rdata[,1],rdata[,2]) # Area under the curve
  return(relauc)
}