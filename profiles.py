import matplotlib.pyplot as plt
import numpy as np
from sklearn import preprocessing
import glob
from obspy import read
from obspy.taup import TauPyModel
from scipy.io.wavfile import write
plt.style.use('seaborn')
## for normalizing traces
min_max_scaler11 = preprocessing.MinMaxScaler(feature_range=(-1, 1))
# tauP earth model setting
model = TauPyModel(model="prem")



for i in range(0,9):

    bigstream=np.zeros(1)

    for component in ['Z','R','T']:


        all_disp=glob.glob(f"data_directory/model_{i}/*{component}.sac")
        ## visualizing
        fig,ax = plt.subplots(1,1,figsize=(10,10))

        p_arrival_times,s_arrival_times,scs_arrival_times,s_arrival_dist = [], [], [], []
        ## loop through all files



        for disp in all_disp:
            # read file
            st =  read(disp)

            bigstream=np.append(bigstream,st[0][0:40000],axis=0)
            
            #dist = kilometers2degrees(st[0].stats.sac['dist']) #convert from km to degree
            dist=st[0].stats.sac['gcarc']
            start=st[0].stats.starttime
            #st=st.slice(start+1000,start+2000)
            #print(dist)
            time_axis = st[0].times()
            data_axis = st[0].data
            evdp = st[0].stats.sac['evdp'] #event depth
            arrivals = model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist, phase_list=["P","S","ScS"]) #get P and S wave arrival time using tauP
            #print(arrivals)
            p_arrival = arrivals[0].time
            s_arrival = arrivals[1].time
            scs_arrival=arrivals[2].time
            #data_axis = min_max_scaler11.fit_transform((data_axis).reshape(-1, 1)) #normalize traces
            plt.ylim([65,85])
            
            if component == 'Z':
                plt.xlim([500,1500])

            if component == 'R':
                plt.xlim([500,1500])

            if component == 'T':
                plt.xlim([1000,2000])
                
            ax.plot(time_axis,1e7*data_axis+dist,color='k',lw=0.5)
            #p_arrival_times.append(5*60)
            #s_arrival_times.append(5*60+(s_arrival-p_arrival))
            #scs_arrival_times.append(5*60+(scs_arrival-s_arrival))
            p_arrival_times.append(p_arrival)
            s_arrival_times.append(s_arrival)
            scs_arrival_times.append(scs_arrival)
            s_arrival_dist.append(2*np.mean(data_axis)+dist)
    
        s_arrival_dist = np.array(s_arrival_dist)
        s_arrival_times = np.array(s_arrival_times)
        p_arrival_times = np.array(p_arrival_times)
        scs_arrival_times=np.array(scs_arrival_times)

        s_arrival_sort_idx = np.argsort(s_arrival_dist)
        s_arrival_dist = s_arrival_dist[s_arrival_sort_idx]
        s_arrival_times = s_arrival_times[s_arrival_sort_idx]
        p_arrival_times = p_arrival_times[s_arrival_sort_idx]
        scs_arrival_times=scs_arrival_times[s_arrival_sort_idx]

        #ax.plot(s_arrival_times,s_arrival_dist,'r--',lw=2) #plot S arrivals
        #ax.plot(p_arrival_times,s_arrival_dist,'r--',lw=2) #plot P arrivals
        #ax.plot(scs_arrival_times,s_arrival_dist,'r--',lw=2) #plot ScS arrivals
        ax.set_title(f"model {i}, {component}-component, up to 5 seconds")
        ax.set_ylabel('Distance in degrees',fontsize=14)
        ax.set_xlabel('Time',fontsize=14)
        fig.savefig(f'./figures/model{i}.{component}.eps')
        plt.close(fig)
        fig.clear()

    #print(bigstream)

    data=bigstream
    scaled=np.int16(data/np.max(np.abs(data))*32767)
    write(f'model_{i}_sonif.wav', 51200,scaled)

    
