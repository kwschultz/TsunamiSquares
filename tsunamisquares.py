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
from mpl_toolkits.basemap import Basemap
import matplotlib.font_manager as mfont
# -------
import quakelib
from os import system
import read_ETOPO1

"""
Tsunami Squares output files are in the following columns:

TIME (secs) |  LON (dec. degrees)  |  LAT (dec. deg.)  |  WATER HEIGHT (meters)  |  BATHY. DEPTH (meters)
  
In the code below: X = LON, Y = LAT, Z = WATER HEIGHT, ALT = BATHYMETRY DEPTH

Water height and bathymetry depth are given in units of meters relative to global mean sea level.
If a square sits on a beach at 5m above sea level and has 4m of water on it, then Z=9m and ALT=5m.

"""


# --------------------------------------------------------------------------------
def make_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP):
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
    
    # I don't know why, but the y-axis is backwards
    ax.invert_yaxis()
    
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    # Increment the time from T_MIN
    TIME = T_MIN

    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = len(np.unique(first_step['lon']))

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
            ALT = this_step['alt'].reshape(-1, Ncols)

            # Plot the surface for this time step
            if surface is None:
                surface = ax.imshow(Z,cmap=cmap,origin='upper',norm=norm,extent=[lon_min,lon_max,lat_max,lat_min],interpolation=interp)
            else:
                surface.set_data(Z)
                
            # Text box with the time
            plt.figtext(0.02, 0.5, 'Time: {:02d}:{:02d}'.format(int(time)/60, int(time)%60), bbox={'facecolor':'yellow', 'pad':5})    
            
                
            writer.grab_frame()
            
            TIME +=T_STEP


# --------------------------------------------------------------------------------
def make_map_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP, save_file):
    # Get ranges
    lon_min,lon_max = sim_data['lon'].min(),sim_data['lon'].max()
    lat_min,lat_max = sim_data['lat'].min(),sim_data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    z_min,z_max = sim_data['z'].min(),sim_data['z'].max()
    cmap = plt.get_cmap('Blues_r')
    norm = mcolor.Normalize(vmin=z_min/60, vmax=-z_min/60)
    interp = 'none'
    landcolor = '#FFFFCC'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    
    # Split the data up into arrays for each time step
    split_data = np.split(sim_data, np.unique(sim_data['time']).shape[0])
    
    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='TsunamiSquares', artist='Matplotlib',
                    comment='Bump in the middle, with accelerations.')
    writer = FFMpegWriter(fps=FPS, metadata=metadata)
    
    # Initialize the frame and axes
    fig = plt.figure()

    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
            llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)

    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    #m.drawcountries()
    #m.drawstates()
    #m.fillcontinents(color=landcolor)
    #m.shadedrelief()
    
    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.95, 0.7, r'water altitude $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    
    # Increment the time from T_MIN
    TIME = T_MIN
    
    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = len(np.unique(first_step['lon']))
    
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
            ALT = this_step['alt'].reshape(-1, Ncols)
            
            # Masked array via conditional, don't color the land unless it has water on it
            masked_data = np.ma.masked_where(np.logical_and(np.array(Z == 0.0),np.array(ALT >= 0.0)), Z)
            
            # Set masked pixels to the land color
            cmap.set_bad(landcolor, 1.0)  # set alpha=0.0 for transparent
            
            # Plot the surface for this time step
            if surface is None:
                surface = m.ax.imshow(masked_data,cmap=cmap,origin='lower',norm=norm,extent=[lon_min,lon_max,lat_max,lat_min],interpolation=interp)
            else:
                surface.set_data(masked_data)
                
            # Text box with the time
            plt.figtext(0.129, 0.82, 'Time: {:02d}:{:02d}'.format(int(time)/60, int(time)%60), bbox={'facecolor':'yellow', 'pad':5})
                
            writer.grab_frame()
        
            TIME +=T_STEP


# =============================================================
def plot_eq_displacements(LLD_FILE, LEVELS, save_file):
    # Read displacement data
    disp_data = np.genfromtxt(LLD_FILE, dtype=[('lat','f8'),('lon','f8'), ('z','f8')],skip_header=3)

    # Data ranges
    lon_min,lon_max = disp_data['lon'].min(),disp_data['lon'].max()
    lat_min,lat_max = disp_data['lat'].min(),disp_data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    z_min,z_max = disp_data['z'].min(),disp_data['z'].max()
    z_lim = max(np.abs(z_min),np.abs(z_max))
    cmap = plt.get_cmap('seismic')
    norm = mcolor.Normalize(vmin=-z_lim, vmax=z_lim)
    interp = 'cubic'
    landcolor = '#FFFFCC'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)

    # Initialize the frame and axes
    fig = plt.figure()
    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    m.fillcontinents(color=landcolor, zorder=0)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'displacement $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)

    # Reshape into matrices
    Ncols = len(np.unique(disp_data['lon']))
    Nrows = len(np.unique(disp_data['lat']))
    
    X = disp_data['lon'].reshape(Nrows, Ncols)
    Y = disp_data['lat'].reshape(Nrows, Ncols)
    Z = disp_data['z'].reshape(Nrows, Ncols)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELS)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(Z <= LEVELS[zero_above]),np.array(Z >= LEVELS[zero_below])), Z)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELS, cmap=cmap, norm=norm, extend='both', zorder=1)

    plt.savefig(save_file,dpi=100)
    print("Saved to "+save_file)

# =============================================================
def bathy_topo_map(LLD_FILE, save_file):
    # Read bathymetry/topography data
    data = np.genfromtxt(LLD_FILE, dtype=[('lat','f8'),('lon','f8'), ('z','f8')],skip_header=3)

    # Data ranges
    lon_min,lon_max = data['lon'].min(),data['lon'].max()
    lat_min,lat_max = data['lat'].min(),data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    cmap = plt.get_cmap('terrain')
    interp = 'none'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    
    # Reshape into matrices
    Ncols = len(np.unique(data['lon']))
    Nrows = len(np.unique(data['lat']))
    
    X = data['lon'].reshape(Nrows, Ncols)
    Y = data['lat'].reshape(Nrows, Ncols)
    Z = data['z'].reshape(Nrows, Ncols)
        
    # catch any nan values
    masked_data = np.ma.masked_invalid(Z)
    cmap.set_bad('red')
    
    # Color limits
    z_min,z_max = masked_data.min(),masked_data.max()
    z_lim = max(np.abs(z_min),np.abs(z_max))
    norm = mcolor.Normalize(vmin=-z_lim, vmax=z_lim)

    # Initialize the frame and axes
    fig = plt.figure()
    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'elevation $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)    
        
    # Plot the contours
    #m.contourf(X, Y, masked_data, 100, cmap=cmap, norm=norm, extend='both', zorder=1)
    m.ax.imshow(masked_data,cmap=cmap,origin='lower',norm=norm,extent=[lon_min,lon_max,lat_max,lat_min],interpolation=interp)

    plt.savefig(save_file,dpi=100)
    print("Saved to "+save_file)
    

# --------------------------------------------------------------------------------
if __name__ == "__main__":
    
    MODE = "animate"
    
    if MODE == "generate":
        # ====== PARSE ETOPO1 FILE, SAVE SUBSET, EVALUATE EVENT FIELD AT THE LAT/LON, SAVE =====
        ETOPO1_FILE = "ETOPO1_Bed_g_gmt4.grd"
        SAVE_NAME = "local/Channel_Islands_largest_subset.txt"
        MODEL     = "../VQModels/UCERF2/ALLCAL2_VQmeshed_3km.h5"
        EVENTS    = "../Desktop/RUNNING/events_greensTrimmed_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_4kyr_24June2015.h5"
        EVID      = 1157
        # Full range
        #MIN_LAT = 33.503
        #MAX_LAT = 34.519
        #MIN_LON = -120.518
        #MAX_LON = -118.883
        # =================================
        # Larger subset
        MIN_LAT = 33.75
        MAX_LAT = 34.3
        MIN_LON = -120.2
        MAX_LON = -119.2
        # --- write grid ------
        lats,lons,bathy = read_ETOPO1.grab_ETOPO1_subset(ETOPO1_FILE,min_lat=MIN_LAT,max_lat=MAX_LAT,min_lon=MIN_LON,max_lon=MAX_LON)
        read_ETOPO1.write_grid(SAVE_NAME,lats,lons,bathy)
        # ---- compute field and write it ------
        system("python ../vq/PyVQ/pyvq/pyvq.py --field_eval  --event_file {} --model_file {} --event_id {} --lld_file {} ".format(EVENTS, MODEL, EVID, SAVE_NAME))
    
    if MODE == "animate":
        sim_file = "local/Channel_Islands_interp_larger_EQ_sample_flatBottom.txt"
        save_file = sim_file.split(".")[0]+"_colorChop_noLines.mp4"
        sim_data = np.genfromtxt(sim_file, dtype=[('time','f8'),('lat','f8'),('lon','f8'), ('z','f8'), ('alt','f8')])
        FPS = 10
        DPI = 100
        T_MAX,T_MIN = sim_data['time'].max(),sim_data['time'].min()
        T_STEP = np.unique(sim_data['time'])[1] - np.unique(sim_data['time'])[0]
        assert T_STEP > 0
        N_STEP = float(T_MAX-T_MIN)/T_STEP
        # Do it
        #make_map_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP, save_file)
        make_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP)

    if MODE == "eq_field_plot":
        Levels = [-.3, -.2, -.1, -.05, -.008, .008, .05, .1, .2, .3]
        plot_eq_displacements("local/Channel_Islands_dispField_event1157.txt",Levels, "disp_map.png")
        
    if MODE == "bathy":
        #Levels = [-3, -.2, -.1, -.05, -.008, .008, .05, .1, .2, .3]
        #bathy_topo_map("local/Channel_Islands.txt",Levels, "bathy_map.png")
        bathy_topo_map("local/Channel_Islands_interp_larger.txt", "bathy_map_interp_larger_imshow.png")
        
    if MODE == "interp":
        # ====== PARSE ETOPO1 FILE, SAVE SUBSET, EVALUATE EVENT FIELD AT THE LAT/LON, SAVE =====
        ETOPO1_FILE = "ETOPO1_Bed_g_gmt4.grd"
        SAVE_NAME = "local/Channel_Islands_interp_larger_subset.txt"
        MODEL     = "../VQModels/UCERF2/ALLCAL2_VQmeshed_3km.h5"
        EVENTS    = "../Desktop/RUNNING/events_greensTrimmed_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_4kyr_24June2015.h5"
        EVID      = 1157
        # Full range
        #MIN_LAT = 33.4
        #MAX_LAT = 34.6
        #MIN_LON = -120.6
        #MAX_LON = -118.8
        # =================================
        # Larger subset for EQ sampling
        MIN_LAT = 33.75
        MAX_LAT = 34.3
        MIN_LON = -120.2
        MAX_LON = -119.2
        # --- write grid ------
        lats,lons,bathy = read_ETOPO1.grab_ETOPO1_subset_interpolated(ETOPO1_FILE,min_lat=MIN_LAT,max_lat=MAX_LAT,min_lon=MIN_LON,max_lon=MAX_LON)
        read_ETOPO1.write_grid(SAVE_NAME,lats,lons,bathy)
        # ---- compute field and write it ------
        system("python ../vq/PyVQ/pyvq/pyvq.py --field_eval  --event_file {} --model_file {} --event_id {} --lld_file {} ".format(EVENTS, MODEL, EVID, SAVE_NAME))


if MODE == "eq_field_eval":
        LLD_NAME = "local/Channel_Islands_interp_larger_subset.txt"
        MODEL     = "../VQModels/UCERF2/ALLCAL2_VQmeshed_3km.h5"
        EVENTS    = "../Desktop/RUNNING/events_greensTrimmed_ALLCAL2_VQmeshed_3km_EQSim_StressDrops_4kyr_24June2015.h5"
        EVID      = 1157
        # ---- compute field and write it ------
        system("python ../vq/PyVQ/pyvq/pyvq.py --field_eval  --event_file {} --model_file {} --event_id {} --lld_file {} ".format(EVENTS, MODEL, EVID, LLD_NAME))




