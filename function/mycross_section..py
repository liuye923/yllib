import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib
import netCDF4 as nc
import xarray as xr
import pandas as pd
#import wrf
from shapely.geometry.polygon import Polygon
#from cartopy import crs as ccrs
#from cartopy.feature import NaturalEarthFeature, OCEAN, LAND, LAKES
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cmaps
import datetime
from copy import deepcopy
import glob
import os
from scipy.interpolate import griddata

#from .commonfunc import CommonFunc

#========================================================================
# LOAD DATA
#========================================================================
#--- MODEL level 1 ---
#class CurrentFunc(CommonFunc):
class CurrentFunc():
    def __init__(self):
        pass
    
    def extract_data_along_line(self, data, xlong, xlat, lat_start, lon_start, lat_end, lon_end, num_points):
        # Generate the points along the line
        latitudes_line = np.linspace(lat_start, lat_end, num_points)
        longitudes_line = np.linspace(lon_start, lon_end, num_points)
        points = np.vstack((latitudes_line, longitudes_line)).T
        
        # Flatten the coordinate grids and the data for interpolation
        points_to_interpolate = np.column_stack((xlat.ravel(), xlong.ravel()))
        data_values = data.ravel()
        
        # Perform the interpolation
        interpolated_values = griddata(points_to_interpolate, data_values, points, method='linear')
        
        return interpolated_values

    def cross_section_averaged_over_rotated_rectangle(self, data, xlat, xlong, vertices, direction, num_lines, num_points):
        if direction==0:
            lat_step_start = (vertices[2][1] - vertices[1][1]) / num_lines
            lat_step_end   = (vertices[3][1] - vertices[0][1]) / num_lines
            lon_step_start = (vertices[2][0] - vertices[1][0]) / num_lines
            lon_step_end   = (vertices[3][0] - vertices[0][0]) / num_lines
            start_point_lat = vertices[1][1]
            start_point_lon = vertices[1][0]
            end_point_lat   = vertices[0][1]
            end_point_lon   = vertices[0][0]        
        averaged_data = []
        lines = []
        for i in range(num_lines):
            # Calculate new start and end points for each line
            new_start_lat = start_point_lat + i * lat_step_start
            new_start_lon = start_point_lon + i * lon_step_start
            new_end_lat   = end_point_lat + i * lat_step_end
            new_end_lon   = end_point_lon + i * lon_step_end
    
            lines.append([[new_start_lon, new_end_lon], [new_start_lat, new_end_lat]])
            
            # Extract data along this line
            line_data = self.extract_data_along_line(
#                data.data, data['XLONG'].data, data['XLAT'].data, 
                data, xlong, xlat, 
                new_start_lat, new_start_lon, new_end_lat, new_end_lon,
                num_points
            )
            
            # Calculate the average of the extracted data
            if len(line_data) > 0:
                averaged_data.append(np.mean(line_data))
            else:
                averaged_data.append(None)
        
#        return averaged_data, lines
        print(averaged_data)
        return np.array(averaged_data)

    def cross_section_averaged_over_rotated_rectangle_3D(self, data, dim, xlat, xlong, vertices, direction, num_lines, num_points):
        from functools import partial
        partial_func = partial(self.cross_section_averaged_over_rotated_rectangle, 
            xlat=xlat, xlong=xlong, 
            vertices=vertices, direction=direction, num_lines=num_lines, num_points=num_points)
            
        return xr.apply_ufunc(partial_func, data,
#            input_core_dims=[*dim],
            input_core_dims=[['south_north','west_east']],
            vectorize=True,# !Important!
            output_core_dims=[['line']],
            )

mycross_section = CurrentFunc()

if __name__ == '__main__':
#    vertices = [
#        [-95.09765264,  29.374     ],
#        [-95.50234736,  29.374     ],
#        [-95.50234736,  30.853     ],
#        [-95.09765264,  30.853     ]
#    ]
    vertices = (np.array([-95.09765264,  29.374     ]),
         np.array([-95.50234736,  29.374     ]),
         np.array([-95.50234736,  30.853     ]),
         np.array([-95.09765264,  30.853     ]))
    da = xr.open_dataarray('../../test_data.nc')
    print(da)
#    averages = mycross_section.cross_section_averaged_over_rotated_rectangle(da.isel(level=1)*100, vertices, 0, 20, 10)
#    averages = mycross_section.cross_section_averaged_over_rotated_rectangle_3D(da.isel(level=1)*100, ['south_north','west_east'], da['XLAT'].data, da['XLONG'].data, vertices, 0, 20, 10)
    averages = mycross_section.cross_section_averaged_over_rotated_rectangle_3D(da*100, ['south_north','west_east'], da['XLAT'].data, da['XLONG'].data, vertices, 0, 20, 10)
    print(averages)

