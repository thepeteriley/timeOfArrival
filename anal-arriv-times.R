# Driver routine to read in the results of the CCMC shock arrival times forecasts, perform some statistics on them, 
# and make the figures found in the paper: 
#
# "Forecasting the Arrival Time of Coronal Mass Ejections: Analysis of the CCMC CME Scoreboard"
#
# to be, or already published in "Space Weather". 
#
# by:
#
# Pete Riley\affil{1}, Leila Mays\affil{2}, Jesse Andries\affil{3}, Tanja Amerstorfer\affil{4},  
# Douglas Biesecker\affil{5}, V\'eronique Delouille\affil{3}, Mateja Dumbovi\'{c}\affil{6,7}, Xueshang Feng\affil{8}, Edmund Henley\affil{9}, Jon A. Linker\affil{1}, Christian M\"{o}stl\affil{4}, Marlon Nu\~{n}ez\affil{10},Vic Pizzo\affil{5}, Manuela Temmer\affil{4}, W.K. Tobiska\affil{11}, C. Verbeke\affil{12}, Matthew J West\affil{3}, and Xinhua Zhao\affil{6}
}
# \affiliation{1}{Predictive Science Inc., San Diego, USA}
# \affiliation{2}{NASA/GSFC, Greenbelt, MD 20771, USA}
# \affiliation{3}{Solar-Terrestrial Center of Excellence, Royal Observatory of Belgium, Ringlaan 3, B-1180 Brussels, Belgium}
# \affiliation{4}{Space Research Institute, Austrian Academy of Sciences, 8042 Graz, Austria}
# \affiliation{5}{Space Weather Prediction Center, NOAA, Boulder, Colorado, USA}
# \affiliation{6}{Institute of Physics, University of Graz, Graz, Austria}
# \affiliation{7}{Hvar Observatory, Faculty of Geodesy, University of Zagreb, Zagreb, Croatia}
# \affiliation{8}{SIGMA Weather Group, State Key Laboratory of Space Weather, National Space Science Center, Chinese Academy of Sciences, Beijing 100190, China}
# \affiliation{9}{Met Office, FitzRoy Road, Exeter, Devon, UK}
# \affiliation{10}{Department of Languages and Computer Sciences, Universidad de M\'{a}laga, Málaga, Spain}
# \affiliation{11}{Space Environment Technologies, Pacific Palisades, CA 90272, USA}
# \affiliation{12}{Centre for Mathematical Plasma-Astrophysics, KU Leuven, Leuven, Belgium}
#
# Written 09/21/17 by PR. Continuously updated through July 18, 2018, at which point it was frozen for publication. 
#
# This is intended to provide support for the analysis presented in the aforementioned paper; however,
# you are free to take/modify this code for their own purpose(s). 
# 
# Disclaimer: 
#
# It includes all the gory details and transient changes so that I, or someone else 
# can recover all/earlier results. It's not supposed to be "production" code. 
#
# Please contact Pete Riley (pete@predsci.com) with any issues or suggested changes - 
# I'm more than happy to help set this up if you want to test YOUR forecasts against 
# those in the CCMC Scoreboard. 
#
# At some point, I'll update this to be a general purpose tool for testing/assessing other 
# forecasts to try to improve on the ones that are currently in the CCMC. 
#
# USAGE: 
#
# Point the "file" variable to the location of the csv file (included with the GitHub repository).
#
# choose one of the plot/analysis options (myPlot)
#
# Dependencies: Uses some R packages that are not installed by default. 
# You'll be prompted to install (via error messages) if you don't have them. 
#
options(digits=3)
require(reshape2)

# read in the CCMC scoreboard data:

#file = "/Users/pete/Dropbox/research2018/arrival-times/ccmc-scoreboard-091517.csv"
file = "/Users/pete/Dropbox/research2018/arrival-times/ccmc-scoreboard-040318.csv"

rawData = read.csv(file)

arrDiff = rawData$Difference..hrs.
arrMethod= rawData$Method
cmeTime= rawData$CME
leadTime = rawData$Lead.Time..hrs.
timeRaw = rawData$Predicted.Shock.Arrival.Time
ddate   = strptime(timeRaw,"%Y-%m-%dT%H:%M",tz="GMT")
year    = as.numeric(format(ddate, "%Y")) 

# identify the methods and how many each one used:
categories <- unique(arrMethod) 
numberOfCategories <- length(categories)
table(arrMethod)

dfDiff = data.frame(arrDiff,ddate,year,arrMethod)

modelsOld = c("WEC (NOAA/SWPC) ", 
           " Average of all Methods ",
           "WEC (GSFC SWRC) ",
           "Ensemble WEC (GSFC SWRC) ", 
           "WEC (Met Office) ", 
           "SPM ")

models = c(" Average of all Methods ", 
           "WEC (GSFC SWRC) ",
           "SIDC ",
           "WEC (NOAA/SWPC) ", 
           "WEC (Met Office) ", 
           "Ensemble WEC (GSFC SWRC) ")

model_colors = c("blue","red","green","orange","purple","brown")

all_colors = 0*year
color_table = rainbow(30)
#color_table = heat.colors(30)

for (i in 1:length(categories)){
  sub_model  = (rawData$Method == categories[i])
  all_colors[sub_model] = i
}

whModel = 1
sub_model  = (rawData$Method == models[whModel])
diff_model = arrDiff[sub_model]
mean_model = mean(diff_model,na.rm=TRUE)
sd_model = sd(diff_model,na.rm=TRUE)
error <- qnorm(0.975)*sd_model/sqrt(length(sub_model[sub_model == T]))
left <- mean_model-error
right <- mean_model+error
print(paste("mean: ",mean_model))
print(paste("SD: ",sd_model))
print(paste("SE: ",sd_model/sqrt(length(sub_model[sub_model == T]))))
print(paste("5%: ",left))
print(paste("95%: ",right))

whichEvents = 2 # 2 (set to 1 for ALL events, set to 2 for best 28)

# set up some plotting panels

# 1 = time series of everything
# 11 for the legend for Figure 1
# 2 = box plot of all forecasts 
# 3 = histograms for top 6 forecasts
# 33 = Lead time for 6 models
# 4 = plot one modeler's results as a function of time
# 5 = Compute statistics for all of the models in the database
# 6 = Latex output listing all models and number of forecasts
# 7 = lead time vs. delta t
# 8 = bubble plot of lead time vs. delta t

myPlot = 11; 
myCex = 1.00
if (myPlot == 1) {
  plot(ddate,arrDiff,xlab="Time (Years)",col=color_table[all_colors],
       ylab=expression(paste(Delta,"t (hrs)",sep='')),cex=1,pch=19,cex.lab=myCex,cex.axis=myCex)
       #lines(seq(as.Date("2010-01-01"), as.Date("2020-01-01"), by="years"),
       #replicate(11,0),lty=2)
       lines(ddate,0.0*arrDiff,lty=2)
#  legend("topleft", legend=categories,col=color_table, cex=0.5, lwd = 1,title="Models")
#    legend("topright", legend=categories,col=color_table, cex=0.68, lwd = 1,title="Models")

}

if (myPlot == 11) {
  plot(1:2, axes=FALSE, frame.plot=F,labels=F,plot=F)
  legend("bottom",legend=categories,ncol=2,cex=0.6,pch=19,col=color_table,pt) # fill=color_table,
}

if (myPlot == 2) {
  boxplot(arrDiff~year, dfDiff,xlab="Time (Years)",ylab=expression(paste(Delta,"t (hrs)",sep='')))
  myboxstats = boxplot(arrDiff~year, dfDiff,plot=F)
  boxstats
#  xtable(boxstats)
}

# set up a vector of a subset of events

if (whichEvents == 1) {
    samplePop = replicate(length(leadTime),T)
    curveMax = 10
    ymax = 15
} else {
  curveMax = 5
  ymax = 6
samplePop = replicate(length(leadTime),T)


uniqueCME = unique(cmeTime)
nUnique = length(uniqueCME)

countMod = replicate(nUnique,0)

for (i in 1:nUnique) {
 samplePop[i] = F
 subUnique = (cmeTime == uniqueCME[i])
 sumMod = 0
 for (j in 1:6) {
 if (any((rawData$Method[subUnique] == models[j]))) {
   sumMod = sumMod + 1
 }
 }
 print(paste("sumMod:",sumMod,sep=""))
 countMod[i] = sumMod
}

sum(countMod == 6) # answer is 28

# these are the events that have all six models in them
best28 = uniqueCME[(countMod == 6)]

# so now create a vector with the events that have 6 forecasting models with true

for (k in 1:length(leadTime)) {
  samplePop[k] = any(cmeTime[k] == best28)
}

}

if (myPlot == 3) {
  
old.par <- par(mfrow=c(2, 3))

for (whModel in 1:6) {
  sub_model  = ((rawData$Method == models[whModel]) & samplePop)
  lead_model = leadTime[sub_model]
  diff_model = arrDiff[sub_model]
  mean_model = mean(diff_model,na.rm=TRUE)
  median_model = median(diff_model,na.rm=TRUE)
  mae_model = mean(abs(diff_model),na.rm=TRUE)
  sd_model   = sd(diff_model,na.rm=TRUE)
  error <- qnorm(0.975)*sd_model/sqrt(length(diff_model[!is.na(diff_model)]))
  left <- mean_model-error
  right <- mean_model+error
  print(models[whModel])
  print(paste("mean: ",mean_model))
  print(paste("MAE: ",mae_model))
  print(paste("SD: ",sd_model))
  print(paste("SE: ",sd_model/sqrt(length(diff_model[!is.na(diff_model)]))))
  print(paste("5%: ",left))
  print(paste("95%: ",right))
  dataForHist = diff_model
  xrange = c(-60,60)
  dataForHist[dataForHist < xrange[1]] = xrange[1] - 2.5
  dataForHist[dataForHist > xrange[2]] = xrange[2] + 2.5
  histOut = hist(dataForHist, 
  #     main=paste(models[whModel],"/SD=",as.character(format(sd_model,digits=4)),sep=''),
       main=paste(models[whModel],sep=''),
       xlab=expression(paste(Delta,"t (hrs)",sep='')), 
       col=model_colors[whModel],
       xlim=xrange,ylim=c(0,ymax),cex.lab=myCex,cex.axis=myCex,
       breaks=seq(-100,100,by=5))
  #hist(lead_model/5.,breaks=histOut$breaks,add=T,col=model_colors[whModel])
  myDen = density(diff_model,na.rm=T)
  lines(myDen$x,myDen$y*curveMax/max(myDen$y),lty=2)
  text(30,ymax,paste("Median =",as.character(format(median_model,digits=3)),sep=''))
  text(30,ymax-ymax/7.5,paste("St.Dev.=",as.character(format(sd_model,digits=4)),sep=''))
  text1 = paste("Model &   Min. & 1st Qu. &  Median &   Mean & MAE & 3rd Qu.  &  Max. & s.d. \\\\")
  message(text1)
    sub_data  = dfDiff$arrDiff[which(dfDiff$arrMethod==categories[i])]
    text2 = paste(models[whModel],' & ',summary(diff_model)[1],' & ',summary(diff_model)[2],' & ',
                  summary(diff_model)[3],' & ',summary(diff_model)[4],' & ',as.character(format(mae_model,digits=3)),' & ',summary(diff_model)[5],' & ',
                  summary(diff_model)[6], ' & ',as.character(format(sd_model,digits=3)),' \\\\ ')
    message(text2)
  }
par(old.par)

}

if (myPlot == 33) { # This is the Lead Time Plot
  
  old.par <- par(mfrow=c(2, 3))
  
  for (whModel in 1:6) {
    sub_model  = ((rawData$Method == models[whModel]) & samplePop)
    lead_model = leadTime[sub_model]
    mean_model = mean(lead_model,na.rm=TRUE)
    median_model = median(lead_model,na.rm=TRUE)
    mae_model = mean(abs(lead_model),na.rm=TRUE)
    sd_model   = sd(lead_model,na.rm=TRUE)
    error <- qnorm(0.975)*sd_model/sqrt(length(lead_model[!is.na(lead_model)]))
    left <- mean_model-error
    right <- mean_model+error
    print(models[whModel])
    print(paste("mean: ",mean_model))
    print(paste("MAE: ",mae_model))
    print(paste("SD: ",sd_model))
    print(paste("SE: ",sd_model/sqrt(length(lead_model[!is.na(lead_model)]))))
    print(paste("5%: ",left))
    print(paste("95%: ",right))
    if (whModel == 1) {
      lead_model = replicate(length(lead_model),0); 
    }
    xrange = c(0,100)
    dataForHist = lead_model
    dataForHist[dataForHist < xrange[1]] = xrange[1] - 2.5
    dataForHist[dataForHist > xrange[2]] = xrange[2] + 2.5
    histOut = hist(dataForHist, 
                   #     main=paste(models[whModel],"/SD=",as.character(format(sd_model,digits=4)),sep=''),
                   main=paste(models[whModel],sep=''),
                   xlab=expression(paste(Delta,"t (hrs)",sep='')), 
                   col=model_colors[whModel],las=2,breaks=seq(-800,800,by=5), 
                   xlim=xrange,ylim=c(0,ymax),cex.lab=myCex,cex.axis=myCex)
    text(75,ymax,paste("Median =",as.character(format(median_model,digits=3)),sep=''))
    text(75,ymax-ymax/8.,paste("St.Dev.=",as.character(format(sd_model,digits=3)),sep=''))
    myDen = density(lead_model)
    lines(myDen$x,myDen$y*ymax/max(myDen$y),lty=2)
    text1 = paste("Model &   Min. & 1st Qu. &  Median &   Mean & MAE & 3rd Qu.  &  Max. & s.d. \\\\")
    message(text1)
    sub_data  = dfDiff$arrDiff[which(dfDiff$arrMethod==categories[i])]
    text2 = paste(models[whModel],' & ',summary(lead_model)[1],' & ',summary(lead_model)[2],' & ',
                  summary(lead_model)[3],' & ',summary(lead_model)[4],' & ',as.character(format(mae_model,digits=3)),' & ',summary(diff_model)[5],' & ',
                  summary(lead_model)[6], ' & ',as.character(format(sd_model,digits=3)),' \\\\ ')
    message(text2)
  }
  par(old.par)
  
}

if (myPlot == 4) {

  #modelToPlot = "WEC (NOAA/SWPC) "
  #modelToPlot = "SIDC "
  modelToPlot = "WEC (Met Office) "
  
  old.par <- par(mfrow=c(2, 3))
  text1 = paste("Year", " &  Min. & 1st Qu.  & Median  &  Mean & MAE & 3rd Qu.   & Max.  &  s.d. & No. Forecasts \\\\")
  message(text1)
  for (i in 2013:2018) {
    sub_data = dfDiff$arrDiff[which((dfDiff$arrMethod==modelToPlot) & (dfDiff$year==i))]
    mae_model = mean(abs(sub_data),na.rm=TRUE)
    sd_model   = sd(sub_data,na.rm=TRUE)
    myHist = hist(sub_data, 
                  #main = "",
                  main=paste(i),
                  #main=paste(modelToPlot,'/',i,sep=''),
                  xlab=expression(paste(Delta,"t (hrs)",sep='')), 
       col="blue",
       xlim=c(-60,60),
       ylim=c(0,8),
       breaks=seq(-60,60,by=10))  
    sum_data = summary(sub_data)
#    print(sum_data)
#    print(paste("MAE: ",mae_model))
    text2 = paste(i,' & ',summary(sub_data,digits=3)[1],' & ',summary(sub_data)[2],' & ',
          summary(sub_data)[3],' & ',summary(sub_data)[4],' & ', 
          as.character(format(mae_model,digits=3)),' & ', summary(sub_data)[5],' & ',
          summary(sub_data)[6], ' & ',as.character(format(sd_model,digits=3)),' & ',
          length(sub_data),' \\\\ ')
    message(text2)

  }    
  par(old.par)
  old.par <- par(mfrow=c(1,1))
}


if (myPlot == 5) {
  text1 = paste("Model", " &  Min. & 1st Qu.  & Median  &  Mean & MAE & 3rd Qu.   & Max.  &  s.d. & L.T. & No. Forecasts  & NA's \\\\")
  message(text1)
  for (i in 1:numberOfCategories) {
    sub_data  = dfDiff$arrDiff[which(dfDiff$arrMethod==categories[i])]
    sub_lead = leadTime[which(dfDiff$arrMethod==categories[i])]
    mae_model = mean(abs(sub_data),na.rm=TRUE)
    sd_model   = sd(sub_data,na.rm=TRUE)
    median_lead = median(sub_lead,na.rm=TRUE)
    text2 = paste(categories[i],' & ',summary(sub_data,digits=3)[1],' & ',summary(sub_data)[2],' & ',
          summary(sub_data)[3],' & ',summary(sub_data)[4],' & ', 
          as.character(format(mae_model,digits=3)),' & ', summary(sub_data)[5],' & ',
          summary(sub_data)[6], ' & ',as.character(format(sd_model,digits=3)),' & ',
          median_lead, ' & ', length(sub_data), ' & ',summary(sub_data)[7],' \\\\ ')
    message(text2)
  }
}

if (myPlot == 6) {
  text1 = paste("Model Name &  No. Forecasts \\\\")
  message(text1)
  for (i in 1:numberOfCategories) {
    sub_data  = dfDiff$arrDiff[which(dfDiff$arrMethod==categories[i])]
    text2     = paste(categories[i], ' & ', length(sub_data),' \\\\ ')
    message(text2)
  }
}

# lead time vs. delta t for all events

if (myPlot == 7) {   
  plot(leadTime,arrDiff,xlab="Lead Time (hrs)",col=color_table[all_colors],xlim=c(0,130),ylim=c(-70,50),
       ylab=expression(paste(Delta,"t (hrs)",sep='')),cex=1,pch=19,cex.lab=myCex,cex.axis=myCex)
   lines(-10:150,replicate(161,0),lty=2)    
}

# bubble plot showing mean error, lead time, with bubbles the size of the number of forecasts

if (myPlot == 8) {   
  plot(leadTime,arrDiff,xlab="Lead Time (hrs)",col=color_table[all_colors],xlim=c(0,70),ylim=c(-30,30),
       ylab=expression(paste(Delta,"t (hrs)",sep='')),cex=1,pch=19,cex.lab=myCex,cex.axis=myCex,
       type="n")
   lines(-10:100,replicate(111,0),lty=2) 
    color_table_trans <- adjustcolor(color_table, alpha.f = 0.3) 
for (i in 1:numberOfCategories) {
    sub_data  = dfDiff$arrDiff[which(dfDiff$arrMethod==categories[i])]
    sub_lead = leadTime[which(dfDiff$arrMethod==categories[i])]
    mae_model = mean(abs(sub_data),na.rm=TRUE)
    sd_model   = sd(sub_data,na.rm=TRUE)
    median_lead = median(sub_lead,na.rm=TRUE)
    points(median_lead,summary(sub_data)[3],pch=19,col=color_table_trans[i],
    cex=1+3.*log10(length(sub_data)))
    text(median_lead,summary(sub_data)[3]+sample(-2:2,1),categories[i],cex=0.4)
  }

}

