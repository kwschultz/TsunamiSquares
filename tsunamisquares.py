#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.animation as manimation
import matplotlib.colorbar as mcolorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
# -------
import quakelib
from os import system
import read_ETOPO1


# --------------------------------------------------------------------------------
def make_animation(sim_data, FPS, DPI, ELEV, AZIM, T_MIN, T_MAX, T_STEP, N_STEP):
    # Get ranges
    lon_min,lon_max = sim_data['lon'].min(),sim_data['lon'].max()
    lat_min,lat_max = sim_data['lat'].min(),sim_data['lat'].max()
    z_min,z_max = sim_data['z'].min(),sim_data['z'].max()
    cmap = plt.get_cmap('Blues_r')
    norm = mcolor.Normalize(vmin=z_min, vmax=z_max)
    interp = 'none'
    
    # Split the data up into arrays for each time step
    split_data = np.split(sim_data, np.unique(sim_data['time']).shape[0])

    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='TsunamiSquares', artist='Matplotlib',
            comment='Bump in the middle, with accelerations.')
    writer = FFMpegWriter(fps=FPS, metadata=metadata)

    # Initialize the frame and axes
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlim(lon_min, lon_max)
    plt.ylim(lat_min, lat_max)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    # Increment the time from T_MIN
    TIME = T_MIN

    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = np.unique(first_step['lon']).shape[0]

    surface = None
    with writer.saving(fig, save_file, DPI):
        for index in range(int(N_STEP)):
            # Get the subset of data corresponding to current time
            this_step = split_data[index]
            time = this_step['time'][0]
            
            print "step: "+str(index)+"  time: "+str(time)+" num_points: "+str(len(this_step))
            assert len(this_step) > 0

            X = this_step['lon'].reshape(-1, Ncols)
            Y = this_step['lat'].reshape(-1, Ncols)
            Z = this_step['z'].reshape(-1, Ncols)

            # Plot the surface for this time step
            if surface is None:
                surface = ax.imshow(Z,cmap=cmap,origin='lower',norm=norm,extent=[lon_min,lon_max,lat_min,lat_max],interpolation=interp)
            else:
                surface.set_data(Z)
                
            # Text box with the time
            plt.figtext(0.02, 0.5, 'Time: {:02d}:{:02d}'.format(int(time)/60, int(time)%60), bbox={'facecolor':'yellow', 'pad':5})    
            
                
            writer.grab_frame()
            
            TIME +=T_STEP



# --------------------------------------------------------------------------------
if __name__ == "__main__":
    # Load TsunamiSquares data
    #sim_file = "accel_middle_bump_renormFractions_LxLy_900_dt20.txt"
    
    MODE = "generate"
    
    if MODE == "generate":
        # ====== PARSE ETOPO1 FILE, SAVE SUBSET, EVALUATE EVENT FIELD AT THE LAT/LON, SAVE =====
        ETOPO1_FILE = "ETOPO1_Bed_g_gmt4.grd"
        SAVE_NAME = "local/Channel_Islands.txt"
        MODEL     = "../Desktop/RUNNING/UCERF2/ALLCAL2_VQmeshed_3km.h5"
        EVENTS    = "../Desktop/RUNNING/events_greensTrimmed_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_4kyr_24June2015.h5"
        EVID      = 1157
        # Full range
        MIN_LAT = 33.503
        MAX_LAT = 34.519
        MIN_LON = -120.518
        MAX_LON = -118.883
        # Inner subset
        #MIN_LAT = 33.874
        #MAX_LAT = 34.137
        #MIN_LON = -119.961
        #MAX_LON = -119.35
        # Larger inner subset for tests
        #MIN_LAT = 33.874
        #MAX_LAT = 34.4
        #MIN_LON = -119.961
        #MAX_LON = -119.35
        # --- write grid ------
        lats,lons,bathy = read_ETOPO1.grab_ETOPO1_subset(ETOPO1_FILE,min_lat=MIN_LAT,max_lat=MAX_LAT,min_lon=MIN_LON,max_lon=MAX_LON)
        read_ETOPO1.write_grid(SAVE_NAME,lats,lons,bathy)
        # ---- compute field and write it ------
        system("python ../vq/pyvq/pyvq/pyvq.py --event_file {} --model_file {} --event_id {} --lld_file {} --field_eval".format(EVENTS, MODEL, EVID, SAVE_NAME))
    
    if MODE == "animate":
        sim_file = "local/TS_Channel_Islands_dispField_event1157_dt10_flatBottom_shorter.txt"
        save_file = sim_file.split(".")[0]+".mp4"
        sim_data = np.genfromtxt(sim_file, dtype=[('time','f8'),('lat','f8'),('lon','f8'), ('z','f8')])
        FPS = 5
        DPI = 100
        ELEV = 20
        AZIM = None
        T_MAX,T_MIN = sim_data['time'].max(),sim_data['time'].min()
        T_STEP = np.unique(sim_data['time'])[1] - np.unique(sim_data['time'])[0]
        assert T_STEP > 0
        N_STEP = float(T_MAX-T_MIN)/T_STEP
        # Do it
        make_animation(sim_data, FPS, DPI, ELEV, AZIM, T_MIN, T_MAX, T_STEP, N_STEP)







