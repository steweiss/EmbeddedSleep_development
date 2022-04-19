butterworth_filter<-
  function(ts,samp_freq=200,
           passband=c(2,40),
           pb_ripple=5,
           stopband=c(1,60),
           sb_ripple=20){
    Fs = samp_freq;                                                     # Sampling Frequency (Guess)
    Fn = Fs/2;                                                          # Nyquist Frequency
    Ts = 1/Fs;                                                          # Sampling Interval
    Wp = passband/Fn;                                                   # Passband
    Ws = stopband/Fn;                                                   # Stopband
    Rp = pb_ripple;                                                     # Passband Ripple
    Rs = sb_ripple;                                                     # Stopband Ripple
    btord <- buttord(Ws,Wp,Rp,Rs)
    b<-butter(btord$n, btord$Wc, "pass")$b
    a<-butter(btord$n, btord$Wc, "pass")$a
    #filter and align ecg according to specs above
    return(py$ECG_Filter(a,b,ts))
  }


# mean_centre<-function(df,L=30*250){
#   mean_df=downsample(df,30*250)
#   sq_space=seq(0,nrow(df),by=L)
#   for(i in 1:ncol(df)){
#     df[,i]=period.apply(df[,i],sq_space,function(x)mean_ing)
#   }
#   
# }
# mean_ing=function(df)df-mean(df)

# ssa.svd<-function(df,L){
#   st_klt_clust<-as.numeric(input$klt_clust)
#   if(st_klt_clust==0){
#     x
#   }else{
#     ## chose based on stationarity - e.g. nutrlan or not not - eigen, additionally töplitz only for stationary 
#     df.ssa<-ssa(df$ts,L=30*250,kind="mssa",svd.method="propack")
#     
#     
#     eig.values<-df.ssa$sigma
#     
#     sig.comps <- grouping.auto.wcor(df.ssa,nclust=2,groups=c(1:50),method="complete")
#     
#     eeg_klt_gr <- grouping.auto.wcor(df.ssa,nclust=as.numeric(20),groups=c(1:50),method="complete")
#     ## REconstructin by Components, e.g. seasonality, sine-cosine frequency components (use period) and noise (by correlation matrix)
#     sig.split <- reconstruct(df.ssa, groups = eeg_klt_gr)
#     plot(sig.split,plot.method = "xyplot")
#     
#     
#     residuals(eeg_klt_re)
#     # klt_dat<-x
#     # for(i in 1:as.numeric(input$klt_clust)){
#     #   klt_dat<-cbind(klt_dat,as.numeric(eeg_klt_re[[as.character(i)]]))
#     #   colnames(klt_dat)[i+1]<-paste("PC-Signal",i)
#     # }
#     # klt_dat<-cbind(klt_dat,"PC-uncovered"=as.numeric(residuals(eeg_klt_re)))
#     # klt_dat<-as.ts(klt_dat,frequency=200)
#     # klt_dat[,-1][,1:(ncol(klt_dat)-(2+input$final_component))]
#   }
# }

ssa.svd<-function(df,L,chans,eigs){
  ## chose based on stationarity - e.g. nutrlan or not not - eigen, additionally töplitz only for stationary
  if(ncol(df)>1){
    df.ssa<-ssa(df[,chans],L=L,neig=eigs,kind="mssa",svd.method="propack")
  }else{
    df.ssa<-ssa(df[,chans],L=L,neig=eigs,kind="1d-ssa",svd.method="propack")
    
  }
  return(df.ssa)
}

grouping.ssa<-function(df,g,nc){
  sig.comps  <- grouping.auto.wcor(df,nclust=nc,groups=g,method="complete")
  return(sig.comps)
  ## Reconstructin by Components,
  ## e.g. seasonality, sine-cosine frequency components (use period) and noise (by correlation matrix)
  #residuals(eeg_klt_re)
  # klt_dat<-x
  # for(i in 1:as.numeric(input$klt_clust)){
  #   klt_dat<-cbind(klt_dat,as.numeric(eeg_klt_re[[as.character(i)]]))
  #   colnames(klt_dat)[i+1]<-paste("PC-Signal",i)
  # }
  # klt_dat<-cbind(klt_dat,"PC-uncovered"=as.numeric(residuals(eeg_klt_re)))
  # klt_dat<-as.ts(klt_dat,frequency=200)
  # klt_dat[,-1][,1:(ncol(klt_dat)-(2+input$final_component))]
}

recons<-function(issa,gclust,stacked=F){
  if(!is.list(gclust)){
    print("Clusters must be list")
    return(issa)
  }
  # residuals first
  tmp=residuals(issa)
  colnames(tmp)[1]<-"Unexplained"
  
  ## test for full recons
  ## remember sum eigenvalues
  gcls=c()
  eig_rank=c()
  for (i in gclust) {
    gcls=c(gcls,i)
    # print(i)
    eig_rank=c(eig_rank,sum(issa$sigma[i]))
    
  }
  
  if(!all(1:max(gcls) %in% gcls)){
    sings=(1:max(gcls))[!(1:max(gcls) %in% gcls)]
    clustsigs=list()
    for(i in 1:length(gclust)){
      clustsigs[[i]]=gclust[[i]]
    }
    singsigs=as.list(sings)
    clustsigs=append(clustsigs,singsigs)
    sig.split <- reconstruct(issa,clustsigs,drop.attributes=F)
    
  }else{
    
    sig.split <- reconstruct(issa,gclust,drop.attributes=F)
    
  }
  
  # if incomplete clustering, reconstruct so
  # if(){
  # sig.split <- reconstruct(df,g,drop.attributes=T)
  # }else{
  #
  # }
  
  
  ## add signals according to eigvalues and return as cumsum timeseries
  if(stacked){
    ##wrong!
    for(i in 1:length(gclust)){
      # tmp2=tmp[,i]
      
      # for(j in 1:i){
        tmp2<- reconstruct(issa,list(unlist(gclust[1:i])),drop.attributes=F)$F1
        # tmp2=tmp2+sig.split[[order(eig_rank,decreasing=T)[j]]]
      # }
      colnames(tmp2)<-paste0("Clust",i)
      tmp=merge(tmp,tmp2)
    }
  }else{
    for(i in 1:length(sig.split)){
      tmp2=sig.split[[order(eig_rank,decreasing=T)[i]]]

      colnames(tmp2)<-paste0("Clust",order(eig_rank,decreasing=T)[i])
      tmp=merge(tmp,tmp2)
    }
    
  }
  
  
  return(tmp)
  
}



kalman_filter<-function(df,autoarma=T){
  if(autoarma){
    buildFun <- function(x) {
      res<-auto.arima(df,allowdrift = F)
      ar_sum<-sum(grepl("ar",names(res$coef)))
      ma_sum<-sum(grepl("ma",names(res$coef)))
      if(ar_sum+ma_sum!=length(res$coef))print("ERROR!!!")
      m <- dlmModARMA(ar=res$coef[1:ar_sum], ma=res$coef[(1+ar_sum):(length(res$coef))])
      return(m)
    }
  }else{
    buildFun <- function(x,p,q) {
      res<-arima(df,order=c(p,0,q))
      ar_sum<-sum(grepl("ar",names(res$coef)))
      ma_sum<-sum(grepl("ma",names(res$coef)))
      if(ar_sum+ma_sum!=length(res$coef))print("ERROR!!!")
      m <- dlmModARMA(ar=res$coef[1:ar_sum], ma=res$coef[(1+ar_sum):(length(res$coef))])
      return(m)
    }
  }
  fit <- dlmMLE(df,parm = rep(0,8), build = buildFun)
  
  dlmExG <- buildFun(fit$par)
  
  filt_dat<-dlmFilter(df, dlmExG)
  
  # p_dat<-dlmSmooth(filt_dat)
  
  return(fil_dat$s[-1])
}

rel_bp<-function(df){
  res_df<-as.data.frame.matrix(matrix(0L,nrow=floor(length(df)/(30*250)),ncol=6))
  colnames(res_df)<-c("Delta (<4 Hz)","Theta (4-8 Hz)","Alpha (8-18 Hz)","Beta (18-30 Hz)","Gamma (30-50 Hz)","Total")
  for(i in 1:(floor(length(df)/(30*250)))){
    ft <- fft(df[(1+(i-1)*(30*250)):((i)*(30*250))])
    fre <-(1:length(ft))*250/length(ft)
    A <- (Re(ft)^2 + Im(ft)^2)^.5
    PSD <- cbind.data.frame(fre,A)
    delta<-sum(PSD[fre<4,"A"])
    theta<-sum(PSD[fre>=4&fre<8,"A"])
    alpha<-sum(PSD[fre>=8&fre<18,"A"])
    beta<-sum(PSD[fre>=18&fre<30,"A"])
    gamma<-sum(PSD[fre>=30&fre<50,"A"])
    total<-delta+theta+alpha+beta+gamma
    res_df[i,]<-c(delta,theta,alpha,beta,gamma,total)
  }
  return(res_df)
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