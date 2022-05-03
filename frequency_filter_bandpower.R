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