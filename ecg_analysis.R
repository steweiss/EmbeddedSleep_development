hrv_analysis<-function(rr,srate,mainDir){
  #create correct data frame and turn on text info mode
  write.table(rr/srate,paste0(mainDir,"/Data/rr.txt"),row.names = F, col.names = F,sep=",",dec=".")
  hrv_df<-CreateHRVData()
  hrv_df<-SetVerbose(hrv_df,Verbose=F)
  
  #load data
  hrv_df = LoadBeatAscii(hrv_df, "rr.txt",RecordPath = paste0(mainDir,"/Data/"))
  
  # calculate non-interpolated RR intervals
  hrv_df = BuildNIHR(hrv_df)
  
  #filter unacceptable data points
  hrv_df=FilterNIHR(hrv_df)
  
  #interpolation neccessary for spectral analysis
  hrv_df.freq = InterpolateNIHR (hrv_df, freqhr = 1)
  hrv_df.freq = CreateFreqAnalysis(hrv_df.freq)
  
  #Calculate and Plot Powerbands7make transparent
  hrv_df.freq = CalculatePowerBand(hrv_df.freq, indexFreqAnalysis= 1,shift=2,size=30,
                                   type = "fourier",
                                   bandtolerance = 0.01, relative = FALSE)
  
  
  #Create and Print Time Analysis
  hrv_df.time = CreateTimeAnalysis(hrv_df, size = 300,interval = 7.8125)
  
  
  # Nonlinear Analysis and Poincare Plot
  hrv_df.nonlin = CreateNonLinearAnalysis(hrv_df)
  
  return(list(hrv_df.freq,hrv_df.time,hrv_df.nonlin))
}


ecg_refs<-function(controls,cases){
  source("C:/Users/Stefan/Nextcloud/Arbeit/MentaLab/Publikation/05_Statistical Analysis/Code/for repo/eXg_filter.R")
  source("C:/Users/Stefan/Nextcloud/Arbeit/MentaLab/Publikation/05_Statistical Analysis/Code/for repo/qrs_detection.R")
  control_dat<-list()
  case_dat<-list()
  rr_df<-list()
  
  for(i in 1:length(controls)){
    
    ###Reading as Numeric Time Series/outsource if ones has been read to read directly processed
    control_dat[i]<-data.frame(as.numeric(as.character(
      read.csv(paste0(getwd(),"/data/",controls[i],".csv"))[,2]
    ))[-1])
    
    ##Naming Datasets for reidentification
    names(control_dat)[i]<-controls[i]
    
    ###Plotting ECG Signal, make toggle plot
    
    times<-ts(control_dat[[i]],frequency=200)
    p<-autoplot(times*10,xlab="Time [sec]", ylab="unfiltered Voltage [mV]")+theme_classic()
    
    # output[[controls[i]]] <- renderPlot({
    ggplotly(p)
    # })
    filts<-eXg_filter(times,
                      samp_freq=200,
                      passband=c(4,40),
                      pb_ripple=5,
                      stopband=c(1,60),
                      sb_ripple=20,
                      "menta")
    
    p<-autoplot(filts[[3]],xlab="Time [sec]", ylab="filtered Voltage [mV]")+theme_classic()
    
    rrs<-qrs_detection(filts)
    for_plot<-data.frame(y=rep(max(filts[[3]])+max(filts[[3]])/10,length(rrs)),x=rrs+1)
    texty<-p+geom_point(data=for_plot,aes(x=x,y=y))
    
    
    rr_df[[i]]<-data.frame(rrs)
    names(rr_df)[i]<-controls[i]
    
    #write.table(rr_df[i],paste0(getwd(),"/data/rr_",cases[i],".txt"),col.names = FALSE,row.names = FALSE)
    colnames(rr_df[[i]])[1]<-"y"
    rr_df[[i]]["x"]<-seq(0,length(rr_df[[i]][,1])-1)
    
    ###Plottin RR distances, make transparency plot
    # output[cases[i]] <- renderPlot({
    ggplot(data=rr_df[[i]],aes(y=x,x=y))+geom_point()+xlab("RR Distance in seconds")+ ylab("")+ggtitle(cases[i])
    # })
    
    
  }
  
  for(j in 1:length(cases)){
    i<-j+length(controls)
    ###Reading as Numeric Time Series/outsource if ones has been read to read directly processed
    case_dat[j]<-data.frame(as.numeric(as.character(
      read.csv(paste0(getwd(),"/data/",cases[j],".csv"))[,2]
    ))[-1])
    
    ##Naming Datasets for reidentification
    names(case_dat)[j]<-cases[j]
    
    times<-ts(case_dat[[j]],frequency=360)
    p<-autoplot(times,xlab="Time [sec]", ylab="V2 Lead position [mV]")+theme_classic()
    
    
    ###Plotting ECG Signal, make toggle plot
    # output[[cases[i]]] <- renderPlot({
    ggplotly(p)
    # })
    
    rr_df[[i]]<-data.frame(qrs_detection(eXg_filter(times,
                                                    samp_freq=360,
                                                    passband=c(2,40),
                                                    pb_ripple=5,
                                                    stopband=c(1,60),
                                                    sb_ripple=20)))
    names(rr_df)[i]<-cases[j]
    #write.table(rr_df[i],paste0(getwd(),"/data/rr_",cases[i],".txt"),col.names = FALSE,row.names = FALSE)
    #rr_df[i]<-data.frame(rr_df[i])
    colnames(rr_df[[i]])[1]<-"y"
    rr_df[[i]]["x"]<-seq(0,length(rr_df[[i]][,1])-1)
    
    ###Plottin RR distances, make transparency plot
    # output[cases[i]] <- renderPlot({
    ggplot(data=rr_df[[i]],aes(y=x,x=y))+geom_point()+xlab("RR Distance in seconds")+ ylab("")+ggtitle(cases[j])
    # })
    
    
  }
  return(rr_df)
}