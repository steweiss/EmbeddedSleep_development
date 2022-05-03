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
