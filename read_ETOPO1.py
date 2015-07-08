#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset

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
    
    
    
if __name__ == "__main__":

    ETOPO1_FILE = "ETOPO1_Bed_g_gmt4.grd"    
    SAVE_NAME = "test.txt"
    MIN_LAT = 36.0
    MAX_LAT = 36.07
    MIN_LON = -125.0
    MAX_LON = -124.92
    lats,lons,bathy=grab_ETOPO1_subset(ETOPO1_FILE,min_lat=MIN_LAT,max_lat=MAX_LAT,min_lon=MIN_LON,max_lon=MAX_LON)
    write_grid(SAVE_NAME,lats,lons,bathy)
    
    
    
    
