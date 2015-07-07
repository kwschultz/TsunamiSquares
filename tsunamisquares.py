#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation


# --------------------------------------------------------------------------------
def make_animation(sim_data, FPS, DPI, ELEV, AZIM, T_MIN, T_MAX, T_STEP, N_STEP):
    # Get ranges
    x_min,x_max = sim_data['x'].min(),sim_data['x'].max()
    y_min,y_max = sim_data['y'].min(),sim_data['y'].max()
    z_min,z_max = sim_data['z'].min(),sim_data['z'].max()

    # Split the data up into arrays for each time step
    split_data = np.split(sim_data, np.unique(sim_data['time']).shape[0])


    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='TsunamiSquares', artist='Matplotlib',
            comment='Bump in the middle, with accelerations.')
    writer = FFMpegWriter(fps=FPS, metadata=metadata)

    # Initialize the frame and axes
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface([],[],[])
    ax.set_xlim3d(x_min, x_max)
    ax.set_ylim3d(y_min, y_max)
    ax.set_zlim3d(z_min, z_max)

    # Increment the time from T_MIN
    TIME = T_MIN

    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = np.unique(first_step['x']).shape[0]

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
sim_file = "accel_middle_bump.txt"
save_file = "TS_"+sim_file.split(".")[0]+".mp4"
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

