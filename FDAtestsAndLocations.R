#program to generate FDA table, plots, and area studies
#written by Naomi Brownstein, May 2020
#updated Nov 24, 2020

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

##########################################33
#load FDA test information: change this to your directory
#updated Nov 24, 2020
##########################################33

setwd("M:\\lab\\Lab_Brownstein\\Research\\ScreeningTests\\Antibodies\\R1\\code_drafts")
tests.full=read.csv("formattedTestsFull.csv",header=TRUE)
tests.combOnly=read.csv("formattedTestsCombOnly.csv",header=TRUE)


############################################################
###prepare data for FDA tests in paper########
############################################################
#parameters (full: all tests)
FDA.sens=cbind(tests.full$sensLCL,tests.full$sens,tests.full$sensUCL)
FDA.spec=cbind(tests.full$specLCL,tests.full$spec,tests.full$specUCL)
FDA.prev=seq(0.005,0.1,0.005)

#parameters (combined only: remove extra IgG and IgM)
FDA.sens.combOnly=cbind(tests.combOnly$sensLCL,tests.combOnly$sens,tests.combOnly$sensUCL)
FDA.spec.combOnly=cbind(tests.combOnly$specLCL,tests.combOnly$spec,tests.combOnly$specUCL)



############################################################
###generate plots for each serology test########
############################################################
#graphs for all tests 
#no need to rerun for subset
collectallpdfnames=rep("",nrow(tests.full))
for(i in 1:nrow(tests.full)){
  FDA.results=fullPPVNPVresults(FDA.sens[i,],FDA.spec[i,],FDA.prev)
  pdfname=paste(tests.full$Test[i],tests.full$Type[i],".pdf")
  ggplotPPVbyPrev(FDA.results,pdfname)
  collectallpdfnames[i]=pdfname
}




############################################################
###output results of tests to send to Excel########
############################################################
#all
#prevalences for PPV/false positives
prev1pct=rep(0.01,nrow(tests.full))
prev5pct=rep(0.05,nrow(tests.full))
prev10pct=rep(0.1,nrow(tests.full))

#subset
prev1pctCombOnly=rep(0.01,nrow(tests.combOnly))
prev5pctCombOnly=rep(0.05,nrow(tests.combOnly))
prev10pctCombOnly=rep(0.1,nrow(tests.combOnly))


#input parameters for p=0.01
calc1=cbind(tests.full$sens,tests.full$spec,prev1pct)
calc1LCL=cbind(tests.full$sensLCL,tests.full$specLCL,prev1pct)
calc1UCL=cbind(tests.full$sensUCL,tests.full$specUCL,prev1pct)

#subset
calc1.combOnly=cbind(tests.combOnly$sens,tests.combOnly$spec,prev1pctCombOnly)
calc1LCL.combOnly=cbind(tests.combOnly$sensLCL,tests.combOnly$specLCL,prev1pctCombOnly)
calc1UCL.combOnly=cbind(tests.combOnly$sensUCL,tests.combOnly$specUCL,prev1pctCombOnly)

#input parameters for p=0.05
calc5=cbind(tests.full$sens,tests.full$spec,prev5pct)
calc5LCL=cbind(tests.full$sensLCL,tests.full$specLCL,prev5pct)
calc5UCL=cbind(tests.full$sensUCL,tests.full$specUCL,prev5pct)

#subset
calc5.combOnly=cbind(tests.combOnly$sens,tests.combOnly$spec,prev5pctCombOnly)
calc5LCL.combOnly=cbind(tests.combOnly$sensLCL,tests.combOnly$specLCL,prev5pctCombOnly)
calc5UCL.combOnly=cbind(tests.combOnly$sensUCL,tests.combOnly$specUCL,prev5pctCombOnly)


#input parameters for p=0.1
calc10=cbind(tests.full$sens,tests.full$spec,prev10pct)
calc10LCL=cbind(tests.full$sensLCL,tests.full$specLCL,prev10pct)
calc10UCL=cbind(tests.full$sensUCL,tests.full$specUCL,prev10pct)

#subset
calc10.combOnly=cbind(tests.combOnly$sens,tests.combOnly$spec,prev10pctCombOnly)
calc10LCL.combOnly=cbind(tests.combOnly$sensLCL,tests.combOnly$specLCL,prev10pctCombOnly)
calc10UCL.combOnly=cbind(tests.combOnly$sensUCL,tests.combOnly$specUCL,prev10pctCombOnly)


#calculate PPV for p=0.01
FDA1=t(apply(calc1,1,getPPVandNPV))
FDA1LCL=t(apply(calc1LCL,1,getPPVandNPV))
FDA1UCL=t(apply(calc1UCL,1,getPPVandNPV))

#subset
FDA1.combOnly=t(apply(calc1.combOnly,1,getPPVandNPV))
FDA1LCL.combOnly=t(apply(calc1LCL.combOnly,1,getPPVandNPV))
FDA1UCL.combOnly=t(apply(calc1UCL.combOnly,1,getPPVandNPV))


#calculate PPV for p=0.05
FDA5=t(apply(calc5,1,getPPVandNPV))
FDA5LCL=t(apply(calc5LCL,1,getPPVandNPV))
FDA5UCL=t(apply(calc5UCL,1,getPPVandNPV))

#subset
FDA5.combOnly=t(apply(calc5.combOnly,1,getPPVandNPV))
FDA5LCL.combOnly=t(apply(calc5LCL.combOnly,1,getPPVandNPV))
FDA5UCL.combOnly=t(apply(calc5UCL.combOnly,1,getPPVandNPV))

#calculate PPV for p=0.1
FDA10=t(apply(calc10,1,getPPVandNPV))
FDA10LCL=t(apply(calc10LCL,1,getPPVandNPV))
FDA10UCL=t(apply(calc10UCL,1,getPPVandNPV))

#subset
FDA10.combOnly=t(apply(calc10.combOnly,1,getPPVandNPV))
FDA10LCL.combOnly=t(apply(calc10LCL.combOnly,1,getPPVandNPV))
FDA10UCL.combOnly=t(apply(calc10UCL.combOnly,1,getPPVandNPV))


#put each estimate with the upper and lower bounds
#full
FDA1withCL=cbind(FDA1LCL[,1],FDA1[,1],FDA1UCL[,1])
FDA5withCL=cbind(FDA5LCL[,1],FDA5[,1],FDA5UCL[,1])
FDA10withCL=cbind(FDA10LCL[,1],FDA10[,1],FDA10UCL[,1])
FDAPPVests=cbind(FDA1withCL,FDA5withCL,FDA10withCL)

#subset
FDA1withCL.combOnly=cbind(FDA1LCL.combOnly[,1],FDA1.combOnly[,1],FDA1UCL.combOnly[,1])
FDA5withCL.combOnly=cbind(FDA5LCL.combOnly[,1],FDA5.combOnly[,1],FDA5UCL.combOnly[,1])
FDA10withCL.combOnly=cbind(FDA10LCL.combOnly[,1],FDA10.combOnly[,1],FDA10UCL.combOnly[,1])
FDAPPVests.combOnly=
  cbind(FDA1withCL.combOnly,FDA5withCL.combOnly,FDA10withCL.combOnly)

#write.csv(FDAPPVests,"FDAPPV.csv")
#write.csv(FDAPPVests.combOnly,"FDAPPVcombOnly.csv")


#convert to false positives
FDAfalsep1=cbind(1-FDA1UCL[,1],1-FDA1[,1],1-FDA1LCL[,1])
FDAfalsep5=cbind(1-FDA5UCL[,1],1-FDA5[,1],1-FDA5LCL[,1])
FDAfalsep10=cbind(1-FDA10UCL[,1],1-FDA10[,1],1-FDA10LCL[,1])
FDAfps=cbind(FDAfalsep1,FDAfalsep5,FDAfalsep10)
#write.csv(round(100*FDAfps,1),"FDAfalsepos.csv")


#subset
FDAfp1.co=cbind(1-FDA1UCL.combOnly[,1],
                    1-FDA1.combOnly[,1],1-FDA1LCL.combOnly[,1])
FDAfp5.co=cbind(1-FDA5UCL.combOnly[,1],
                    1-FDA5.combOnly[,1],1-FDA5LCL.combOnly[,1])
FDAfp10.co=cbind(1-FDA10UCL.combOnly[,1],
                     1-FDA10.combOnly[,1],1-FDA10LCL.combOnly[,1])
FDAfps.co=cbind(FDAfp1.co,FDAfp5.co,FDAfp10.co)
colnames(FDAfps.co)=c("fp1LCL","fp1PE","fp1UCL",
                   "fp5LCL","fp5PE","fp5UCL",
                   "fp10LCL","fp10PE","fp10UCL")

FDAfalsep1.co=paste(as.character(format(100*(1-FDA1.combOnly[,1]),digits=1,nsmall=1)),
                    "(",
                    as.character(format(100*(1-FDA1UCL.combOnly[,1]),digits=1,nsmall=1)),
                    ",",
                    as.character(format(100*(1-FDA1LCL.combOnly[,1]),digits=1,nsmall=1)),
                    ")"
                      )
FDAfalsep5.co=paste(as.character(format(100*(1-FDA5.combOnly[,1]),digits=1,nsmall=1)),
                    "(",
                    as.character(format(100*(1-FDA5UCL.combOnly[,1]),digits=1,nsmall=1)),
                    ",",
                    as.character(format(100*(1-FDA5LCL.combOnly[,1]),digits=1,nsmall=1)),
                    ")"
)
FDAfalsep10.co=paste(as.character(format(100*(1-FDA10.combOnly[,1]),digits=1,nsmall=1)),
                    "(",
                    as.character(format(100*(1-FDA10UCL.combOnly[,1]),digits=1,nsmall=1)),
                    ",",
                    as.character(format(100*(1-FDA10LCL.combOnly[,1]),digits=1,nsmall=1)),
                    ")"
)
####################################
#new code: automate tables
FDAtable.combOnly=
  cbind(
    as.character(tests.combOnly$Test),
    as.character(tests.combOnly$testType),
    as.character(tests.combOnly$sens_and_CI),
    as.character(tests.combOnly$spec_and_CI),
    as.character(tests.combOnly$N),  
    FDA1
    FDAfalsep5.co,
    FDAfalsep10.co,
      #FPR y% (CI) (for each prev chosen), 1 dec
    as.character('\\')  )
write.table(FDAtable.combOnly,"FDAtable.combOnly.txt",sep="&",row.names=FALSE)

#results only (no latex formatting)
FDAtable.combOnly.num=cbind(tests.combOnly,FDAfps.co)
write.csv(FDAtable.combOnly.num,"FDAtable.combOnly.num.csv")

#output names for easy input into LaTeX
command=cbind(#" ",
  #as.character("hfill"),
  #as.character("subfloat"),sep="\\"),
  #as.character("\\hfill\\subfloat\["),
  as.character(tests.full$Test),
  #as.character("]{\includegraphics[width=0.45\textwidth]{'"),
  as.character(collectallpdfnames)#,
  #as.character("\label{"),
#  as.character(tests.full$Test),
#  as.character("}}")
#sep=""
)
#\subfloat[Cellex]{\includegraphics[width=0.45\textwidth]{Cellex Combined .pdf} \label{Cellex}}  
colnames(command)
write.table(command,"pdfnames.txt",row.names=FALSE,col.names=c("company","pdf"))
##################################

#########################################################
###############################Area studies
#########################################################

#############
#santa clara
#99.5% specificity with a 95% PI (98.8%,99.8)
#81.8% sensitivity (64.2%.91.0%
santaclara.spec=c(.988,.995,0.998)#specificity was 99.5\% (95CI 99.2-99.7\%) and 
santaclara.sens=c(0.642,0.818,0.91)
santaclara.prev=seq(0,0.05,0.001)

santaclara.PVs=fullPPVNPVresults(santaclara.sens,santaclara.spec,santaclara.prev)
ggplotPPVbyPrev(santaclara.PVs,"santaclaraResults.pdf")

############NY###################
#statewide 13.9\%; 
#Regionally, Long Island at 16.7, 
#New York City at 21.2, 
#Westchester, Rockland 11.7 and 
#rest of state, 3.6
#update May 2 included now
#1.2,1.9,2,2: relatively unaffected areas
prev.NY=c(0.012, 0.019, 0.022, 0.024, 0.03, 0.036,0.06,0.117,0.138,0.139,0.167,0.199,0.212)
#Test: Wadsworth IgG
#specificity 93-100 (Ann link)
#FDA Wadsworth 98.8 CI 97.3,99.5
spec.NY=c(0.93,0.973,0.988,0.995)

#sensitivity: not listed 
#FDA sensitivity (wadsworth): .88 (.805,.928)
sens.NY=c(0.805,0.88,0.928)
NY.PVs=fullPPVNPVresults(sens.NY,spec.NY,prev.NY)
ggplotPPVbyPrev(NY.PVs,"NewYorkResults.pdf") #NewYorkPPVplot

########chelsea
#chelsea: Biomedomics
#https://www.biomedomics.com/?fldl=3050
Chelsea.sens=.8866 
Chelsea.spec=.9063
prev.Chelsea=0.315
Chelsea.PVs=fullPPVNPVresults(Chelsea.sens,Chelsea.spec,prev.Chelsea)
Chelsea.PVs
