#program for sensitivity analysis
#compare prev vs test positivity

############################################################
#install/load necessary packages
############################################################

library(MKmisc)
library(ggplot2)

############################################################
#functions
############################################################

#general function to get PPV and NPV from a vector
#should be input sens, then spec, then prev
getPPVandNPV=function(x){
  #predValues(x$sens.read,x$spec.ready,x$prev.ready)
  predValues(x[1],x[2],x[3])
}


#input vectors of sensitivity, specificity, and prev
#output results (PPV, NPV, sens, spec, prev)
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

#Plot PPV by prevalence and output results to pdf
ggplotPPVbyPrev=function(dataset,myfilename){
  newplot=
    ggplot(data = dataset,mapping = aes(x = p, y = PPV,color=factor(sens),shape=factor(spec))) + 
    geom_point()+
    geom_line()+
    labs(x="Prevalence",color="Sensitivity",shape="Specificity")
  print(newplot)
  ggsave(newplot,file=myfilename)
}



#new function
getprevfromser=function(p.ser,sens,spec){
  p.est=(p.ser-(1-spec))/(sens-(1-spec))
  return(p.est)
}

#do it for all rows
prevallrows=function(xmat){
  getprevfromser(xmat[1],xmat[2],xmat[3])
}

#
#put everything together
estimatepfromser=function(PPVresults){
  parms=cbind(PPVresults$p,PPVresults$sens,PPVresults$spec)
  #print(parms)
  p.better=apply(parms,1,prevallrows)
  #print(p.better)
  #print(dim(parms))
  #print(dim(p.better))
  #print(length(p.better))
  finalp=cbind(parms,p.better)
  #print(finalp)
  colnames(finalp)=cbind("p.ser","sens","spec","p.better")
  #now we have the results of NPV and PPV with other variables contained in the 
  #object called final
  final.pdat=as.data.frame(finalp)
  return(final.pdat)
}


#############Chelsea###########
Chelsea.sens=.8866 
Chelsea.spec=.9063
prev.Chelsea=0.315
Chelsea.PVs=fullPPVNPVresults(Chelsea.sens,Chelsea.spec,prev.Chelsea)
##########these are the differences#############
#calculate updated prevalence estimate
Chelsea.phat=estimatepfromser(Chelsea.PVs)
#new prevalence estimate minus positivity
Chelsea.phat$p.better-Chelsea.phat$p.ser
#new matrix
fornewPV.Chelsea=cbind(Chelsea.phat$sens,Chelsea.phat$spec,Chelsea.phat$p.better)
Chelsea.PVhat=t(apply(fornewPV.Chelsea,1,getPPVandNPV))
#updated PPV
Chelsea.PVhat

#differences
Chelsea.PVhat[,1]-Chelsea.PVs$PPV

###############New York############
#old estimates
prev.NY=c(0.012, 0.019, 0.022, 0.024, 0.03, 0.036,0.06,0.117,0.138,0.139,0.167,0.199,0.212)
spec.NY=c(0.93,0.973,0.988,0.995)
sens.NY=c(0.805,0.88,0.928)
NY.PVs=fullPPVNPVresults(sens.NY,spec.NY,prev.NY)

#calculate updated prevalence estimate
NY.phat=estimatepfromser(NY.PVs)
#new matrix
fornewPV.NY=cbind(NY.phat$sens,NY.phat$spec,NY.phat$p.better)

#how many are prevalenc estimates are positive
prevgood=fornewPV.NY[fornewPV.NY[,3]>0,] #120
negprev=fornewPV.NY[fornewPV.NY[,3]<=0,] #36
#differences in prevalence when phat is positive
prevdiffs=NY.phat$p.better[fornewPV.NY[,3]>0]-NY.phat$p.ser[fornewPV.NY[,3]>0]
hist(prevdiffs)
summary(prevdiffs)

#info with new PPVs
NY.PVhat=t(apply(prevgood,1,getPPVandNPV))
newPPVinfo=cbind(NY.PVhat[,1],prevgood)

#extract only matching PPVs from old estimate
oldPPVinfo=NY.PVs[fornewPV.NY[,3]>0,]

#difference in PPV
NYdiffPPV=newPPVinfo[,1]-oldPPVinfo$PPV
hist(NYdiffPPV)
summary(NYdiffPPV)

#new graph
NYnewforPlot=as.data.frame(newPPVinfo)
names(NYnewforPlot)=c("PPV","sens","spec","p")
ggplotPPVbyPrev(NYnewforPlot,"sensitivityNY.pdf")
