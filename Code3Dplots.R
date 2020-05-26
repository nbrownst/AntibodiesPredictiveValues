#program to generate overall 3/4D plot
#written by Naomi Brownstein, May 2020

############################################################
#install/load necessary packages
############################################################

## Use BiocManager to install limma
BiocManager::install("limma")

#function for PPV and NPV
library(MKmisc)

#allows the 3D interactive plot
library(plotly)



############################################################
#functions
############################################################

#function to calculate NPV from a vector of length 3 (sens, spec, and prev)
getPPVandNPV=function(x){
  #predValues(x$sens.read,x$spec.ready,x$prev.ready)
  predValues(x[1],x[2],x[3])
}

#full function to get NPV and PPV from a set of vectors of sens, spec, p
#note each should be a vector, but they can have different lengths
#the function looks at all possibilities for each element of each vector
fullPPVNPVresults=function(sens.input,spec.input,prev.input){
  #from here on it calculates PPV and NPV
  numsensparms=length(sens.input)
  numspecparms=length(spec.input)
  numprevparms=length(prev.input)
  
  sens.repeat=rep(sens.input,numspecparms*numprevparms)
  sens.ready=sens.repeat
  
  spec.repeat=rep(spec.input,numsensparms)
  ordered.spec1=spec.repeat[order(spec.repeat)]
  ordered.spec2=rep(ordered.spec1,numprevparms)
  spec.ready=ordered.spec2
  
  prev.repeat=rep(prev.input,numsensparms*numspecparms)
  ordered.prev=prev.repeat[order(prev.repeat)]
  prev.ready=ordered.prev
  
  #length(rep(sens.input,length(spec.input)))
  parms=cbind(sens.ready,spec.ready,prev.ready)
  
  
  results=t(apply(parms,1,getPPVandNPV))
  final=cbind(results,parms)
  colnames(final)=cbind("PPV","NPV","sens","spec","p")
  #now we have the results of NPV and PPV with other variables contained in the 
  #object called final
  final.dat=as.data.frame(final)
  return(final.dat)
}

#plot function: takes in results and outputs the PPV figure 
makePPVplot=function(myresults){
  fig <- plot_ly(myresults, x = ~spec, y = ~p, z = ~PPV,
                 marker = list(color = ~sens, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'Specificity'),
                                     yaxis = list(title = 'Prevalence'),
                                     zaxis = list(title = 'PPV')))
  return(fig)  
}


#plot function: takes in results and outputs the NPV figure 
makeNPVplot=function(myresults){
  figN <- plot_ly(myresults, x = ~spec, y = ~p, z = ~NPV,
                  marker = list(color = ~sens, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
  figN <- figN %>% add_markers()
  figN <- figN %>% layout(scene = list(xaxis = list(title = 'Specificity'),
                                       yaxis = list(title = 'Prevalence'),
                                       zaxis = list(title = 'NPV')))
  return(figN)               
}

############################################################
###input data for 3D plots in paper########
############################################################

#note: you change these parameters to get a different plot
prev.input=seq(0.01 ,0.3,0.005)
sens.input=seq(0.8,1,by=0.01)
spec.input=seq(0.9,1,by=0.001)

############################################################
###get results and produce plots
############################################################

results1=fullPPVNPVresults(sens.input,spec.input,prev.input)
PPVplot=makePPVplot(results1)
NPVplot=makeNPVplot(results1)

#####print plots################
PPVplot
NPVplot
