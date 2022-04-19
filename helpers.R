to_day<-function(base, now){
  res=difftime(
    as.Date(base)+1,
    as.POSIXct(paste0(as.Date(base),
                      now))+60*60,units="min")<=
    difftime(as.Date(base,tz="CET")+1,
             base+60*60,units="min")
  return(res)
}


downsample<-function(df,dsrate){
    res=xts()
    sq_space=seq(0,nrow(df),by=dsrate)
    org=df
    df2=df
    df2[!is.na(df2)] <- 1
    df2[is.na(df2)] <- 0
    df[is.na(df)] <- 0
    for(i in 1:ncol(df)){
       # tmp=period.apply(df[,i],sq_space,mean2)
      # tmp_NA=period.apply(df[,i],seq(0,nrow(df),by=dsrate),sumisna)
      # tmp=tmp/(dsrate-tmp_NA)
      
      # df2[,i]=ifelse(is.na(df[,i]),0,1)
      # df[,i]=ifelse(is.na(df[,i]),0,df[,i])
      tmp=period.sum(df[,i],sq_space)/period.sum(df2[,i],sq_space)
      #tmp=tmp/tmp2
      # tmp[length(tmp)]=tmp[length(tmp)]/dsrate
      tmp[length(tmp)]=mean2(org[(nrow(org)-(nrow(org)%% dsrate)):nrow(org),i])
      res=merge(res,tmp)
    }
    colnames(res)<-colnames(df)
    df=res
  return(df)
}
mean2<-function(x)mean(x,na.rm=T)
# sumisna<-function(x)sum(is.na(x))
# sum2<-function(x)sum(x,na.rm=T)

signals2EDF <- function(signals, srs, patId, recId, snames, startdate, starttime,  unitz, 
                        prefilterings, filename,printmes=F){
  if(file.exists(filename)) {file.create(filename,overwrite=T)
    }else{file.create(filename,overwrite=F)}
  fid<-file(filename,'wb')
  numSignals <- length(srs)
  
  ## Mentalab Rounding
  # signals=round(signals,2)
  
  # Calculate recording length (in seconds)
  recLength <- (-1)
  for (k in 1:numSignals){
    recLength <- max(recLength, length(signals[,k])/srs[k])
  }
  recLength <- ceiling(recLength)# Ceil to an integer number
  
  # Asuming 1s databloack duration. Check possible incompatibilities
  blockSizeBytes <- sum(2*srs)
  
  if(blockSizeBytes > 61440){
    print('Yet to be implemented: Signals cannot fit on a 1s datablock. Check for (other block size possibilities)')
  } else {
    blockSize <- 1
    numBlocks <- recLength
  }
  
  general_header_size <- 256#bytes
  one_signal_header_size <- 256#bytes
  
  # Write edf
  
  # FIXED HEADER
  header.version <- 0
  header.local_patient_identification <- patId
  header.local_recording_identification <- recId
  header.startdate_recording <- startdate
  header.starttime_recording <- starttime
  header.num_signals <- numSignals
  header.num_bytes_header <- general_header_size + one_signal_header_size * numSignals
  header.reserved <- ''
  header.duration_data_record <- blockSize
  header.num_data_records <- numBlocks
  
  trimAndFillWithBlanks<-function(txt, maxLength, justify='left'){
    # if(nargin == 2)justify = 'left'
    if(nchar(txt) > maxLength){result = substr(txt,1,maxLength)
    }else{
      if(justify== 'right'){
        result<-paste0(paste0(rep(" ",(maxLength-length(txt))),collapse=""),txt)}
      else{ 
        result<-paste0(txt,paste0(rep(" ",(maxLength-nchar(txt))),collapse=""))
      }
    }
    return(result)
  }
  
  writeBin(charToRaw(trimAndFillWithBlanks(header.version,8)),fid,endian = "little")# version
  writeBin(charToRaw(trimAndFillWithBlanks(header.local_patient_identification,80)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.local_recording_identification,80)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.startdate_recording,8)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.starttime_recording,8)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.num_bytes_header,8)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.reserved,44)),fid,endian="little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.num_data_records,8)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.duration_data_record,8)),fid,endian = "little")
  writeBin(charToRaw(trimAndFillWithBlanks(header.num_signals,4)),fid,endian = "little")
  
  # SIGNAL DEPENDENT HEADER
  signalOffsets <- rep(0, numSignals)# In bytes
  header.signals_info=data.frame(label=rep("",numSignals),
                                 transducer_type=rep("",numSignals),
                                 physical_dimension=rep("",numSignals),
                                 physical_min=rep("",numSignals),
                                 physical_max=rep("",numSignals),
                                 digital_min=rep("",numSignals),
                                 digital_max=rep("",numSignals),
                                 prefiltering=rep("",numSignals),
                                 num_samples_datarecord=rep("",numSignals),
                                 reserved=rep("",numSignals),
                                 sample_rate=rep("",numSignals),
                                 signalOffset=rep("",numSignals),stringsAsFactors = F)
  
  # for(i in 1:ncol(header.signals_info)){
  #   header.signals_info[colnames(header.signals_info)[i]]<-as.character(header.signals_info[colnames(header.signals_info)[i]])
  # }
  
  #Mentalab values
  Pmins<-floor(unlist(lapply(signals,min)))
  # Pmins<-floor(unlist(lapply(signals,max))[seq(2, ncol(signals)*5, 5)])
  Pmaxs<-ceiling(unlist(lapply(signals,max)))
  # Pmaxs<-ceiling(unlist(lapply(signals,quantile))[seq(4, ncol(signals)*5, 5)])
  Dmins<-rep(-32768,numSignals)
  Dmaxs<-rep(32767,numSignals)
  
  if( (  log(    (max(abs(Dmins))+max(abs(Dmaxs))+1) ,2  ) / (max(abs(Pmins))+max(abs(Pmaxs)))    )  > 16 ){ print("Error, Accucracy to high, not 16 bit, either shrink distance !Dmin to Dmax!, or Pmin to Pmax")}
  
  for (k in 1:numSignals){
    header.signals_info[k,"label"] <- snames[k]
    header.signals_info[k,"transducer_type"] <- ''
    header.signals_info[k,"physical_dimension"] <- unitz[k]
    
    header.signals_info[k,"physical_min"] <- Pmins[k]
    header.signals_info[k,"physical_max"] <- Pmaxs[k]
    header.signals_info[k,"digital_min"] <- Dmins[k]
    header.signals_info[k,"digital_max"] <- Dmaxs[k]
    
    header.signals_info[k,"prefiltering"] <- prefilterings[k]
    header.signals_info[k,"num_samples_datarecord"] <- as.character(as.numeric(srs[k])*as.numeric(blockSize))
    header.signals_info[k,"reserved"] <- ''
    # NOTE: The two following are not specific EDF header fields, but are practical for EDF handling
    header.signals_info[k,"sample_rate"] <- as.character(as.numeric(header.signals_info[k,"num_samples_datarecord"])/ 
      as.numeric(header.duration_data_record))
    if (k > 1){
      signalOffsets[k] <- signalOffsets[k - 1] + 2 * as.numeric(header.signals_info[k - 1,"num_samples_datarecord"])
    }
    header.signals_info[k,"signalOffset"] <- signalOffsets[k]
  }
  
  # Write signal-dependent header to file
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(header.signals_info[k,"label"],16)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(header.signals_info[k,"transducer_type"],80)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(header.signals_info[k,"physical_dimension"],8)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(as.character(header.signals_info[k,"physical_min"]),8)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(as.character(header.signals_info[k,"physical_max"]),8)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(as.character(header.signals_info[k,"digital_min"]),8)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(as.character(header.signals_info[k,"digital_max"]),8)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(header.signals_info[k,"prefiltering"],80)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(as.character(header.signals_info[k,"num_samples_datarecord"]),8)),fid,endian="little")
  }
  for (k in 1:numSignals){
    writeBin(charToRaw(trimAndFillWithBlanks(header.signals_info[k,"reserved"],32)),fid,endian="little")
  }
  
  # Check data starting point
  current_position <-seek(fid)# in bytes
  if (header.num_bytes_header!= current_position){
    print('Something wrong could be happening: unexpected position at the beginning of the first data block')
  }
  bytes_full_data_record <- 2 * sum(as.numeric(unlist(header.signals_info["num_samples_datarecord"])))
  
  # DATA WRITING
  
  
  close(fid)
  
  for (k in 1:numBlocks){
    if((k %% 100)==0 & printmes==T)print(k)
    # Initialize datablock
    data <- rep(0, bytes_full_data_record/2)# Num samples per data record
    
    for (k1 in 1:numSignals){
      offsetSignal <- (k - 1) * as.numeric(header.signals_info[k1,"num_samples_datarecord"]) + 1
      onsetSignal <- min(offsetSignal + as.numeric(header.signals_info[k1,"num_samples_datarecord"]) - 1,
                         length(signals[,k1]))
      offsetDataBlock <- as.numeric(header.signals_info[k1,"signalOffset"])/2 + 1
      onsetDataBlock <- offsetDataBlock +  onsetSignal-offsetSignal
      data[offsetDataBlock:onsetDataBlock] <- as.integer(Dmins[k1] + (Dmaxs[k1] - Dmins[k1]) * 
                                                           ((signals[offsetSignal:onsetSignal,k1] - Pmins[k1])/(Pmaxs[k1] - Pmins[k1])))
    }
    if(max(int_to_unit(data))>65536)print("Error")
    py$writeEDFbin(fn,data)
  }
}

freq_filter<-function(x,srate,n,bw,eegs){
  for(i in eegs){
    x[,i]<-py$Notch_Filter(srate, 3, n, 0.5, 3, 'butter', x[,i]) 
    
    x[,i]<- butterworth_filter(x[,i],samp_freq=srate,
                               passband=bw,
                               pb_ripple=round((bw[2]-bw[1])/10),
                               stopband=round(bw*c(-0.5,0.5)+bw),
                               sb_ripple=round((bw[2]-bw[1])/5))
  } 
  return(x)
}


update_nodes<-function(df,df2){
  upd<-data.frame(id=c((nrow(df)+1):(nrow(df)+ncol(df2))),
                  label=colnames(df2),
                  group=rep("physical",ncol(df2)))
  df<-rbind(df,upd)
  return(df)    
}

update_nodes_components<-function(nodes,dssa){
  nodes= subset(nodes,group!="Cluster")
  nodes=subset(nodes,group!="Component")
  res=rbind(nodes,data.frame(id=(nrow(nodes)+1):(nrow(nodes)+ncol(dssa$U)),
                             label=paste0("Comp",seq(1,ncol(dssa$U))),
                             group=rep("Component",ncol(dssa$U))))
  return(res)
}


update_nodes_recons<-function(nodes,rssa){
  res=rbind(nodes,data.frame(id=(nrow(nodes)+1):(nrow(nodes)+length(rssa)),
                             label=paste0("Recon",seq(1,length(rssa))),
                             group=rep("Reconstruction",length(rssa))))
  return(res)}

update_edges_cluster<-function(nodes,edges,gssa){
  nodes2=subset(nodes,group=="Cluster")
  nodes3=subset(nodes,group=="Component")
  edges=edges[!(edges[,"from"] %in% nodes2$id |
    edges[,"to"] %in% nodes2$id|
    edges[,"from"] %in% nodes3$id|
    edges[,"to"] %in% nodes3$id),]
  
  for (i in 1:length(gssa)){
    # Cluster Name - To
    clustn=paste0("Clust",names(gssa)[i])
    
    # Find cluster label
    clustid=nodes$id[nodes$label==clustn]
    
    # Component Name - from
    compn=paste0("Comp",gssa[[names(gssa)[i]]])
    
    ## Find label of component
    compid=nodes$id[nodes$label %in% compn]
    
    ## Write new edge
    edges=bind_rows(edges,
                    data.frame(id=(nrow(edges)+1):(nrow(edges)+length(compn)),
                               from=compid,to=rep(clustid,length(compn)),stringsAsFactors = F))
  }
  for(i in 1:length(gssa)){
    # Cluster Name - To
    clustn=paste0("Clust",names(gssa)[i])
    
    # Find cluster label
    clustid=nodes$id[nodes$label==clustn]
    
    ## Find label of component
    eegid=nodes$id[nodes$label %in% "EEG"]
    
    edges=bind_rows(edges,
                    data.frame(id=(nrow(edges)+1):(nrow(edges)+length(eegid)),
                               from=eegid,to=rep(clustid,length(eegid)),stringsAsFactors = F))
  }
  
  return(edges)
}

update_nodes_cluster<-function(nodes,gssa){
  nodes=subset(nodes,group!="Cluster")
  res=rbind(nodes,data.frame(id=(nrow(nodes)+1):(nrow(nodes)+length(gssa)),
                            label=paste0("Clust",seq(1,length(gssa))),
                            group=rep("Cluster",length(gssa))))
  return(res)
}


