#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation

# Load TsunamiSquares data
sim_file = "test_out.txt"
sim_data = np.genfromtxt(sim_file, dtype=[('time','f8'),('x','f8'),('y','f8'), ('z','f8')])
FPS = 5
DPI = 100
T_MAX,T_MIN = sim_data['time'].max(),sim_data['time'].min()
T_STEP = np.unique(sim_data['time'])[1] - np.unique(sim_data['time'])[0]
assert T_STEP > 0
N_STEP = float(T_MAX-T_MIN)/T_STEP

# Get ranges
x_min,x_max = sim_data['x'].min(),sim_data['x'].max()
y_min,y_max = sim_data['y'].min(),sim_data['y'].max()
z_min,z_max = sim_data['z'].min(),sim_data['z'].max()


# Initialize movie writing stuff
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='TsunamiSquares', artist='Matplotlib',
        comment='I hope this works!')
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
with writer.saving(fig, "TS_movie_test.mp4", DPI):
    for _ in range(int(N_STEP)):
        # Clear the plot
        oldplot = surface
        if oldplot is not None:
            ax.collections.remove(oldplot)
        
        # Get the subset of data corresponding to current time
        this_step = sim_data[sim_data['time']==TIME]
        
        print "step: "+str(_)+"  time: "+str(TIME)+" num_points: "+str(len(this_step['x']))
        #assert len(this_step['x']) > 0

        X = this_step['x'].reshape(-1, Ncols)
        Y = this_step['y'].reshape(-1, Ncols)
        Z = this_step['z'].reshape(-1, Ncols)

        # Plot the surface for this time step
        surface = ax.plot_surface(X,Y,Z,color='cyan')
        writer.grab_frame()
        
        TIME +=T_STEP







