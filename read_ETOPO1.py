#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata

def get_subset_indices(min_lat, max_lat, min_lon, max_lon, lats, lons):
    # These are the indices of the closest lat/lon values to
    # (min_lat, max_lat, min_lon, max_lon)
    indices = []
    indices.append(int((np.abs(np.array(lats)-min_lat)).argmin()))
    indices.append(int((np.abs(np.array(lats)-max_lat)).argmin()))
    indices.append(int((np.abs(np.array(lons)-min_lon)).argmin()))
    indices.append(int((np.abs(np.array(lons)-max_lon)).argmin())) 
    # return [min_lat_index, max_lat_index, min_lon_index, max_lon_index]
    return indices


def grab_ETOPO1_subset(file_name, min_lat, max_lat, min_lon, max_lon):
    
    ETOPO1 = Dataset(file_name, 'r')
    lons = ETOPO1.variables["x"][:]
    lats = ETOPO1.variables["y"][:]
    
    # Grab indices for max/min lat/lon bounds
    minLat, maxLat, minLon, maxLon = get_subset_indices(min_lat, max_lat, min_lon, max_lon, lats, lons)
    bathy = ETOPO1.variables["z"][minLat:maxLat,minLon:maxLon]
    lons,lats = np.meshgrid(lons[minLon:maxLon],lats[minLat:maxLat])
    print("== Selected {} points ({}x{}) from {}".format(bathy.size,bathy.shape[1],bathy.shape[0],file_name))   
    print("---- Lats: {} to {},   Lons: {} to {}".format(min_lat, max_lat, min_lon, max_lon))
            
    return lats,lons,bathy
    
    
def write_grid(out_file_name, lats, lons, bathy):
    outfile = open(out_file_name, 'w')
    
    outfile.write("# N_lats\n")
    outfile.write("# N_lons\n")
    outfile.write("{} {}\n".format(lons.shape[0],lons.shape[1]))
    outfile.write("##################################\n")
    # Write vertices from top left (Northwest) to bottom right (Southeast)
    for i in list(reversed(range(lons.shape[0]))):
        for j in range(lons.shape[1]):
            outfile.write("{}\t{}\t{}\n".format(lats[i][j],lons[i][j],bathy[i][j]))
    outfile.close()
    print("output written to {}".format(out_file_name))
    
    
    
def grab_ETOPO1_subset_interpolated(file_name, min_lat, max_lat, min_lon, max_lon, factor=3):
    debug = False

    ETOPO1 = Dataset(file_name, 'r')
    lons = ETOPO1.variables["x"][:]
    lats = ETOPO1.variables["y"][:]

    # Extend the bounds to add a buffer, ensures that interpolation has enough data on the boundaries
    min_lat_buff = min_lat*0.98
    max_lat_buff = max_lat*1.02
    min_lon_buff = min_lon*1.02
    max_lon_buff = max_lon*0.98
    
    # Grab indices for larger buffered max/min lat/lon bounds
    minLat_buff, maxLat_buff, minLon_buff, maxLon_buff = get_subset_indices(min_lat_buff, max_lat_buff, min_lon_buff, max_lon_buff, lats, lons)
    lons_buff,lats_buff = np.meshgrid(lons[minLon_buff:maxLon_buff],lats[minLat_buff:maxLat_buff])
    bathy_buff = ETOPO1.variables["z"][minLat_buff:maxLat_buff,minLon_buff:maxLon_buff]

    # Grab indices for requested lat/lon bounds
    minLat, maxLat, minLon, maxLon = get_subset_indices(min_lat, max_lat, min_lon, max_lon, lats, lons)
    lons,lats = np.meshgrid(lons[minLon:maxLon],lats[minLat:maxLat])
    bathy = ETOPO1.variables["z"][minLat:maxLat,minLon:maxLon]
    
    print("== Selected {} points ({}x{}) from {}".format(bathy.size,bathy.shape[1],bathy.shape[0],file_name))   
    print("---- Lats: {} to {},   Lons: {} to {}".format(min_lat, max_lat, min_lon, max_lon))
    
    # Flatten everything into 1D arrays
    bathy_flat = bathy_buff.flatten()
    lons_buff_flat = lons_buff.flatten()
    lats_buff_flat = lats_buff.flatten()
    points = zip(lons_buff_flat,lats_buff_flat)
    
    # Create the grid for interpolation
    grid_lon_vals = np.linspace(min_lon, max_lon, num=int(len(lons)*factor))
    grid_lat_vals = np.linspace(min_lat, max_lat, num=int(len(lats)*factor))
    grid_lons, grid_lats = np.meshgrid(grid_lon_vals, grid_lat_vals)

    if debug:
        print("buffered")
        print(minLat_buff, maxLat_buff, minLon_buff, maxLon_buff)
        print(min_lat_buff, max_lat_buff, min_lon_buff, max_lon_buff)
        print(lons_buff)
        print(" ")
        print("original")
        print(minLat, maxLat, minLon, maxLon)
        print(min_lat, max_lat, min_lon, max_lon)

    # Interpolate
    bathy_interp = griddata(points, bathy_flat, (grid_lons, grid_lats), method='cubic')        
            
    print("** Interpolated {} points ({}x{})".format(bathy_interp.size,bathy_interp.shape[1],bathy_interp.shape[0]))
    return grid_lats, grid_lons, bathy_interp

    
    
    
    

    
    
    
    
