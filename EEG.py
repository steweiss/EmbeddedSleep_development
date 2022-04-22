def Notch_Filter(fs, band, freq, ripple, order, filter_type, data):
  if(freq==0):
    return data
  else:
    from scipy.signal import iirfilter
    from scipy.signal import lfilter
    nyq  = fs/2.0
    low  = freq - band/2.0
    high = freq + band/2.0
    low  = low/nyq
    high = high/nyq
    b, a = iirfilter(order, [low, high], rp=ripple, btype='bandstop',
                     analog=False, ftype=filter_type)
    filtered_data = lfilter(b, a, data)
    return filtered_data

def import_explorepy(fn):
  from explorepy.tools import bin2csv
  bin2csv(fn, do_overwrite=True)
  
  
def glob_bp(df,srate):
  from yasa import bandpower
  return(bandpower(df,win_sec=30, relative=True,sf=srate,
              bands=[(0.5, 4, 'Delta'), (4, 8, 'Theta'), (8, 12, 'Alpha'),
                     (12, 16, 'Sigma'), (16, 30, 'Beta'), (30, 40, 'Gamma')],
              kwargs_welch=dict(average='median', window='hamming')))
      


def Bandpass_Filter(time, band, freq,rpass,rstop,  order, filter_type, data):
    from scipy.signal import iirfilter
    from scipy.signal import filtfilt
    fs   = 1/time
    nyq  = fs/2.0
    low  = freq - band/2.0
    high = freq + band/2.0
    low  = 0.1/nyq
    high = 0.4/nyq
    b, a = iirfilter(order, [low, high],rp=rpass, rs=rstop, btype='bandpass',
                     analog=False, ftype=filter_type)
    filtered_data = filtfilt(b,a, data,axis=0)
    return filtered_data

# import numpy as np
# a = r.filt_coeff[0]
# b = r.filt_coeff[1]

def ECG_Filter(a,b, data):
    from scipy.signal import iirfilter
    from scipy.signal import filtfilt
    filtered_data = filtfilt(b,a, data,axis=0)
    return filtered_data

def Wiener_Filter(x):
  from scipy.signal import wiener
  return wiener(x)


def Explore_Devices():
  from explorepy.tools import bt_scan
  res = bt_scan()
  return res

def Connect_Explore(x):
  from explorepy import Explore
  from explorepy.explore import connect
  explorer = Explore()
  explorer.connect(device_name=x)  # Put your device Bluetooth name
  dat=explorer.acquire()
  return dat

def Start_Recording(explorer):
  from explorepy import Explore
  from explorepy.explore import acquire
  dat=explorer.acquire()
  return dat

def Sleepbrain():
    from visbrain.gui import Sleep
    Sleep().show()

def Exg_plot(dat,srate):
    #https://plot.ly/python/v3/ipython-notebooks/mne-tutorial/
    #https://plot.ly/javascript/ribbon-plots/
    import matplotlib.pyplot as plt
    import plotly.plotly as py
    from plotly import tools
    from plotly.graph_objs import Layout, YAxis, Scatter, Annotation, Annotations, Data, Figure, Marker, Font
    step = 1. / ncol(dat)
    times=row(dat/srate)
    kwargs = dict(domain=[1 - step, 1], showticklabels=False, zeroline=False, showgrid=False)
    # create objects for layout and traces
    layout = Layout(yaxis=YAxis(kwargs), showlegend=False)
    traces = [Scatter(x=times, y=data.T[:, 0])]
    # loop over the channels
    for ii in range(1, n_channels):
        kwargs.update(domain=[1 - (ii + 1) * step, 1 - ii * step])
        layout.update({'yaxis%d' % (ii + 1): YAxis(kwargs), 'showlegend': False})
        traces.append(Scatter(x=times, y=data.T[:, ii], yaxis='y%d' % (ii + 1)))
    # add channel names using Annotations
    annotations = Annotations([Annotation(x=-0.06, y=0, xref='paper', yref='y%d' % (ii + 1),
                                      text=ch_name, font=Font(size=9), showarrow=False)
                          for ii, ch_name in enumerate(ch_names)])
    layout.update(annotations=annotations)
    # set the size of the figure and plot it
    layout.update(autosize=False, width=1000, height=600)
    fig = Figure(data=Data(traces), layout=layout)
    py.iplot(fig, filename='shared xaxis')
    
    
def sleep_spindles(df,srate):
    from yasa import spindles_detect
    from yasa import spindles_detect_multi
    # if(df.shape[1]!=1):
    #   a=spindles_detect_multi(df, srate, ch_names=ch_n, multi_only=True, remove_outliers=True)
    # else:
    a=spindles_detect(df, srate)
    return(a)
  # # Single-channel full command (shows all the default implicit parameters)
  # yasa.spindles_detect(data, sf, hypno=None, include=(1, 2, 3),
  #                      freq_sp=(12, 15), duration=(0.5, 2), freq_broad=(1, 30),
  #                      min_distance=500, downsample=True,
  #                      thresh={'rel_pow': 0.2, 'corr': 0.65, 'rms': 1.5},
  #                      remove_outliers=False)
  #
  # # Multi-channels detection on NREM sleep only (requires an hypnogram)
  # yasa.spindles_detect_multi(data, sf, ch_names, hypno=hypno)
  #
  # # Multi-channels detection on N2 sleep only with automatic outlier rejection
  # yasa.spindles_detect_multi(data, sf, ch_names, hypno=hypno, include=(2), remove_outliers=True)
  #

def slow_waves(df,srate):
    from yasa import spindles_detect
    return(yasa.sw_detect(df, srate))
  # # Single-channel full command (shows all the default implicit parameters)
  # yasa.sw_detect(data, sf, hypno=hypno, include=(2, 3), freq_sw=(0.3, 3.5),
  #                dur_neg=(0.3, 1.5), dur_pos=(0.1, 1), amp_neg=(40, 300),
  #                amp_pos=(10, 150), amp_ptp=(75, 400), downsample=True,
  #                remove_outliers=False)
  #
  # # Multi-channel slow-waves detection on N2 + N3 sleep only (requires an hypnogram)
  # yasa.sw_detect_multi(data, sf, ch_names, hypno=hypno)
  #

def rapid_eyes(loc,roc,srate):
    from yasa import spindles_detect
    return(yasa.rem_detect(loc, roc, sf))
  #
  # # On REM sleep only + all implicit parameters
  # yasa.rem_detect(loc, roc, sf, hypno=hypno, include=4, amplitude=(50, 325),
  #                 duration=(0.3, 1.5), freq_rem=(0.5, 5), downsample=True,
  #                 remove_outliers=False)
      
def spec_power(df,srate,L):
    from scipy.signal import welch
    from yasa import bandpower_from_psd
    return(yasa.bandpower_from_psd(welch(df, srate, nperseg=int(L * srate)), bands=[(1, 4, 'Delta'), (4, 8, 'Theta'),(8, 18, 'Alpha'),(18, 30, 'Beta'),(30, 50, 'Gamma')]))

def writeEDFbin(fn,df):
  import numpy as np
  newFileBytes= np.uint16(df)
  newFile = open(fn, "ab")
  # write to file
  newFile.write(newFileBytes)


def qrs_detect(unfiltered_ecg,fs):
  from ecgdetectors import Detectors
  detectors = Detectors(fs)
  r_peaks = detectors.pan_tompkins_detector(unfiltered_ecg)
  return(r_peaks)
  
def explore_df():
  from pylsl import StreamInlet, resolve_stream
  import pandas as pd
  import numpy as np

  # first resolve an EEG stream on the lab network
  streams = resolve_stream('type', 'ExG')

  # create a new inlet to read from the stream
  inlet = StreamInlet(streams[0])
  sample, timestamp = inlet.pull_sample()
  stamp=np.array(timestamp)
  exg=np.array(sample)
  sample, timestamp = inlet.pull_sample()
  stamp=np.array(timestamp)
  stamp=np.append(timestamp,stamp)
  dfe=np.array(sample)
  exg=np.append([exg],[dfe],axis=0)
  i=1
  while i<250:
    # get a new sample (you can also omit the timestamp part if you're not
    # interested in it)
    sample, timestamp = inlet.pull_sample()
    dft=np.array(timestamp)
    stamp=np.append(timestamp,stamp)
    dfe=np.array(sample)
    exg=np.append(exg,[dfe],axis=0)
    i=i+1
    
  return(exg)

def connect_explore():
  import explorepy
  explore = explorepy.Explore()
  explore.connect(device_name="Explore_1432") # Put your device Bluetooth name
  explore.push2lsl()
  return()
