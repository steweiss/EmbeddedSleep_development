#################### SHOULD BE UPDATED IN A BASE STARTER SCRIPT ####################
# if(shiny_running()){
init.Date=Sys.time()
init.time=strftime(init.Date, format = "%H:%M")
base.nodes<<-data.frame(
  id = c(1:7),
  label=c("EEG","ECG","EOG","EMG","Notch","High Pass","Low Pass"),
  group=c(rep("biological",4),rep("digital",3)),stringsAsFactors = F)

# init.ExG=import_raw("/srv/shiny-server/UItest/EDF/Stefan W, 32 yo, Apnoe, no Medication_.edf")
# init.srate=init.ExG[["srate"]]
# init.ExG=round(data.frame(init.ExG[["signals"]]),0)
if(demoMode==F){
  ## A small converted data set
  init.ExG<<-round(data.frame(fread(paste0(mainDir,"/Data/DATA021_eeg_sub.csv"), header=T)[,-1])*1000*1000,0)
  init.srate<<-250
  init.channels<<-colnames(init.ExG)
  
  base.nodes<<-data.frame(
    id = c(1:7),
    label=c("EEG","ECG","EOG","EMG","Notch","High Pass","Low Pass"),
    group=c(rep("biological",4),rep("digital",3)),stringsAsFactors = F)
  
  init.nodes<<-update_nodes(base.nodes,init.ExG)
  biosigs<-subset(init.nodes,group=="biological")
  physigs<-subset(init.nodes,group=="physical")
  digsigs<-subset(init.nodes,group=="digital")
  
  init.edges=data.frame(from = physigs$id[5:8], to = biosigs$id[1])
  init.edges<-rbind(init.edges,data.frame(from = physigs$id[1], to = biosigs$id[2]))
  init.edges<-rbind(init.edges,data.frame(from = physigs$id[2:3], to = biosigs$id[3]))
  init.edges<-rbind(init.edges,data.frame(from = physigs$id[4], to = biosigs$id[4]))
  init.edges<-rbind(init.edges,data.frame(from = biosigs$id[1], to = digsigs$id[2:3]))
  init.edges$id<-1:nrow(init.edges)
  init.edges$title=""
  init.edges$value=1
  
  
  init.eeg.fits=findexg(init.nodes,init.edges,1)
  init.ecg.fits=findexg(init.nodes,init.edges,2)
  init.eog.fits=findexg(init.nodes,init.edges,3)
  init.emg.fits=findexg(init.nodes,init.edges,4)
  
  init.freq_filter.temp=freq_filter(init.ExG,init.srate,50,c(2,30),init.eeg.fits)
  colnames(init.freq_filter.temp)<-paste0("Filtered",colnames(init.freq_filter.temp))
  init.freq_filter<-as.xts(init.freq_filter.temp,
                           seq(as.POSIXct("2020-03-03 00:40:00"),
                               by=1/init.srate,
                               length=nrow(init.freq_filter.temp)))
  
  init.bandpower<-rel_bp(init.freq_filter.temp[,5])
  init.freq.bandpower<-as.xts((init.bandpower[,1:5]/init.bandpower[,6])*100,seq(as.POSIXct("2020-03-03 00:40:00"), by=30, length=nrow(init.bandpower)))
  
  init.stages=data.frame(
    "Sleep State (W, S1, S2, S3, REM)"=c("W"),
    "Time (xx:xx or xx:xx:xx)"=as.POSIXct(paste(as.Date(as.POSIXct("2020-03-03 00:40:00"))+1,"00:50:00")),
    stringsAsFactors = F)
  
  # init.stages=data.frame(
  #   "Sleep State (W, S1, S2, S3, REM)"=c("W"),
  #   "Time (xx:xx or xx:xx:xx)"=as.POSIXct("2020-03-03 23:55:00"),
  #                        stringsAsFactors = F)
  
  # init.stages=data.frame("Sleep State (W, S1, S2, S3, REM)"=c("W","S1","S2"),"Time (xx:xx or xx:xx:xx)"=c("00:08","00:30","01:02"),stringsAsFactors = F)
  
  # init.stages=data.frame("Stage / Event" = rep("Stage",3),
  #                        "Sleep State (W, S1, S2, S3, REM)"=c("W","S1","S2"),
  #                        "Time (xx:xx or xx:xx:xx)"=c("00:08","00:30","01:02"),
  #                        stringsAsFactors = F)
  
  ## For Visual Display only, optimization routine for quick display
  init.orig.chan=init.freq_filter[,paste0("Filtered",init.eeg.fits)]
  init.orig.ssa=init.orig.chan
  init.orig.disp=init.orig.chan
  init.disp=downsample(init.freq_filter[,paste0("Filtered",init.eeg.fits)],250*60)
  init.resol=c(first(index(init.freq_filter)),last(index(init.freq_filter)))
  
  
  ## Singular Spectrum Analysis initial call
  # tmp=round(downsample(window(init.freq_filter,start=first(index(init.freq_filter)),end=(first(index(init.freq_filter))+60*5)),10),6)
  tmp=downsample(window(init.freq_filter,start=as.POSIXct("2020-03-03 00:40"),end=as.POSIXct("2020-03-03 00:58")),25)
  init.ssa=ssa.svd(window(tmp[,paste0("Filtered","ch5")],start=first(index(tmp)),end=first(index(tmp))+60*60),30*10,1,50)
  
  ## Update graph display, currently commented
  #init.nodes=update_nodes_components(init.nodes,init.ssa)
  init.ssa_g=grouping.ssa(init.ssa,g=1:(ncol(init.ssa$U)),nc=3)
  
  ## Update graph display, currently commented
  #init.nodes=update_nodes_cluster(init.nodes,init.ssa_g)
  #init.edges=update_edges_cluster(init.nodes,init.edges,init.ssa_g)
  
  ## Reconstruction, only for display, currently soft implemented only
  init.recons=recons(init.ssa,init.ssa_g,stacked=F)
  #}
  
  # init.nodes=update_nodes_recons(init.nodes)
}else{
  init.nodes=base.nodes
  init.edges=data.frame("from","to")
  init.eeg.fits=c()
  init.ecg.fits=c()
  init.emg.fits=c()
  init.eog.fits=c()
  
}


if(shiny_running()){
  
  ## Reactive Environement for graph network
  graph_data <<- reactiveValues(
    nodes = init.nodes,
    edges=init.edges,
    EEG=init.eeg.fits,
    ECG=init.ecg.fits,
    EOG=init.eog.fits,
    EMG=init.emg.fits
  )
  
  ## Main Reactive Object
  if(demoMode==F){
    ExG<<-reactiveValues(
      ts=as.xts(init.ExG,seq(as.POSIXct("2020-03-03 00:40:00"), by=1/init.srate, length=nrow(init.ExG))),
      srate=init.srate,
      channels=init.channels,
      freq=init.freq_filter,
      bandpower=init.freq.bandpower,
      stages=init.stages,
      orig.chan=init.orig.chan,
      orig.ssa=init.orig.ssa,
      orit.events=init.orig.ssa,
      orig.disp=init.orig.disp,
      disp=init.disp,
      resol=init.resol,
      posix=as.POSIXct("2020-03-03 00:40:00"),
      sight="time",
      newfile=F,
      ssa=init.ssa,
      gssa=init.ssa_g,
      rssa=init.recons,
      events=data.frame(),
      df_res=data.frame(),
      df_res2=data.frame(),
      pan=0,
      hold=1,
      a=0,
      recording=F
    )
  }else{
    init.ExG=as.data.frame(rep(0,250*60*60*8))
    colnames(init.ExG)="Empty"
    init.stages=
      init.freq_filter=as.xts(init.ExG,seq(init.Date, by=1/250, length=nrow(init.ExG)))
    init.orig.chan=downsample(as.xts(init.ExG,seq(init.Date, by=1/250, length=nrow(init.ExG))),250*60)
    init.stages=data.frame(
      "Sleep State (W, S1, S2, S3, REM)"=c("W"),
      "Time (xx:xx or xx:xx:xx)"=init.Date,
      stringsAsFactors = F)
    
    ExG<<-reactiveValues(
      ts=as.xts(init.ExG,seq(init.Date, by=1/250, length=nrow(init.ExG))),
      srate=250,
      channels=c("Empty"),
      freq=init.freq_filter,
      bandpower=data.frame(),
      stages=init.stages,
      orig.chan=init.orig.chan,
      orig.ssa=init.orig.chan,
      orit.events=init.orig.chan,
      orig.disp=init.orig.chan,
      disp=init.orig.chan,
      resol=c(first(init.freq_filter),last(init.freq_filter)),
      posix=init.Date,
      sight="bird",
      newfile=F,
      events=data.frame(),
      df_res=data.frame(),
      df_res2=data.frame(),
      pan=0,
      hold=1,
      a=0
    )
    
  }
}
