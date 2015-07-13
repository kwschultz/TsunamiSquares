#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation
import matplotlib.colorbar as mcolorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --------------------------------------------------------------------------------
def make_animation_2D(sim_data, FPS, DPI, ELEV, AZIM, T_MIN, T_MAX, T_STEP, N_STEP):
    # Get ranges
    lon_min,lon_max = sim_data['lon'].min(),sim_data['lon'].max()
    lat_min,lat_max = sim_data['lat'].min(),sim_data['lat'].max()
    z_min,z_max = sim_data['z'].min(),sim_data['z'].max()
    cmap = plt.get_cmap('Blues_r')
    norm = mcolor.Normalize(vmin=z_min, vmax=z_max)
    #norm = mcolor.Normalize(vmin=-10, vmax=10)
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
            
            print "step: "+str(index)+"  time: "+str(TIME)+" num_points: "+str(len(this_step))
            assert len(this_step) > 0

            X = this_step['lon'].reshape(-1, Ncols)
            Y = this_step['lat'].reshape(-1, Ncols)
            Z = this_step['z'].reshape(-1, Ncols)

            # Plot the surface for this time step
            if surface is None:
                surface = ax.imshow(Z,cmap=cmap,origin='lower',norm=norm,extent=[lon_min,lon_max,lat_min,lat_max],interpolation=interp)
            else:
                surface.set_data(Z)
                
            writer.grab_frame()
            
            TIME +=T_STEP


def make_animation_3D(sim_data, FPS, DPI, ELEV, AZIM, T_MIN, T_MAX, T_STEP, N_STEP):
    # Get ranges
    lon_min,lon_max = sim_data['lon'].min(),sim_data['lon'].max()
    lat_min,lat_max = sim_data['lat'].min(),sim_data['lat'].max()
    z_min,z_max = sim_data['z'].min(),sim_data['z'].max()

    # Split the data up into arrays for each time step
    split_data = np.split(sim_data, np.unique(sim_data['time']).shape[0])

    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='TsunamiSquares', artist='Matplotlib',
            comment='Bump in the middle, with accelerations.')
    writer = FFMpegWriter(fps=FPS, metadata=metadata)

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface([],[],[])
    ax.set_xlim3d(lon_min, lon_max)
    ax.set_ylim3d(lat_min, lat_max)
    ax.set_zlim3d(z_min, z_max)

    # Increment the time from T_MIN
    TIME = T_MIN

    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = np.unique(first_step['lon']).shape[0]

    surface = None
    with writer.saving(fig, save_file, DPI):
        for index in range(int(N_STEP)):
            # Clear the last surface to plot updated surface
            previous_plot = surface
            if previous_plot is not None:
                ax.collections.remove(previous_plot)
            
            # Get the subset of data corresponding to current time
            this_step = split_data[index]
            
            print "step: "+str(index)+"  time: "+str(TIME)+" num_points: "+str(len(this_step))
            assert len(this_step) > 0

            X = this_step['lon'].reshape(-1, Ncols)
            Y = this_step['lat'].reshape(-1, Ncols)
            Z = this_step['z'].reshape(-1, Ncols)

            # Plot the surface for this time step
            ax.view_init(elev=ELEV, azim=AZIM)
            surface = ax.plot_surface(X,Y,Z,color='#6699FF')
            writer.grab_frame()
            
            TIME +=T_STEP




# --------------------------------------------------------------------------------
# Load TsunamiSquares data
sim_file = "accel_middle_bump_renormFractions_LLDasXYZ_initialV.txt"
save_file = "TS_"+sim_file.split(".")[0]+"_fullRange.mp4"
sim_data = np.genfromtxt(sim_file, dtype=[('time','f8'),('lat','f8'),('lon','f8'), ('z','f8')])
FPS = 2
DPI = 100
ELEV = 20
AZIM = None
T_MAX,T_MIN = sim_data['time'].max(),sim_data['time'].min()
T_STEP = np.unique(sim_data['time'])[1] - np.unique(sim_data['time'])[0]
assert T_STEP > 0
N_STEP = float(T_MAX-T_MIN)/T_STEP
# Do it
make_animation_2D(sim_data, FPS, DPI, ELEV, AZIM, T_MIN, T_MAX, T_STEP, N_STEP)







