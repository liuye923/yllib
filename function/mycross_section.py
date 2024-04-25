import numpy as np
import netCDF4 as nc
import xarray as xr
import pandas as pd
from shapely.geometry.polygon import Polygon
import datetime
from copy import deepcopy
import glob
import os
from scipy.interpolate import griddata

#========================================================================
# LOAD DATA
#========================================================================
#--- MODEL level 1 ---
class CurrentFunc():
    def __init__(self):
        pass

    def calculate_geographic_rectangle(self, midpoint1, midpoint2, width_km):
        # Constants
        EARTH_RADIUS_KM = 6371.0  # Approximate radius of the earth in kilometers
        DEG_TO_KM = np.pi * EARTH_RADIUS_KM / 180  # Factor to convert degrees to kilometers
        
        # Convert midpoints to numpy arrays
        mp1 = np.array(midpoint1)
        mp2 = np.array(midpoint2)
        
        # Vector from the first to the second midpoint
        vector = mp2 - mp1
        
        # Calculate half-width in degrees
        half_width_deg = (width_km / 2) / DEG_TO_KM
        
        # Normalize the vector
        vector_norm = np.linalg.norm(vector)
        unit_vector = vector / vector_norm
        
        # Rotate 90 degrees to find perpendicular
        perp_vector = np.array([-unit_vector[1], unit_vector[0]]) * half_width_deg
        
        # Calculate vertices assuming square (midpoints are parallel to sides)
        vertex1 = mp1 + perp_vector
        vertex2 = mp1 - perp_vector
        vertex3 = mp2 - perp_vector
        vertex4 = mp2 + perp_vector
        
        # Calculate angle from East
        angle_from_east = np.degrees(np.arctan2(vector[1], vector[0]))
        
        return (vertex1[::-1], vertex2[::-1], vertex3[::-1], vertex4[::-1]), angle_from_east
    
    def extract_data_along_line(self, data, 
        xlong, xlat, 
        lat_start, lon_start, lat_end, lon_end,
        num_points, method='linear'):
        # Generate the points along the line
        latitudes_line = np.linspace(lat_start, lat_end, num_points)
        longitudes_line = np.linspace(lon_start, lon_end, num_points)
        points = np.vstack((latitudes_line, longitudes_line)).T
        
        # Flatten the coordinate grids and the data for interpolation
        points_to_interpolate = np.column_stack((xlat.ravel(), xlong.ravel()))
        data_values = data.ravel()
        
        # Perform the interpolation
        interpolated_values = griddata(points_to_interpolate, data_values, points, method=method)
        
        return interpolated_values

    def calculate_lines(self, vertices, num_lines, direction):
        if direction==0:
            lat_step_start = (vertices[2][1] - vertices[1][1]) / num_lines
            lat_step_end   = (vertices[3][1] - vertices[0][1]) / num_lines
            lon_step_start = (vertices[2][0] - vertices[1][0]) / num_lines
            lon_step_end   = (vertices[3][0] - vertices[0][0]) / num_lines
            start_point_lat = vertices[1][1]
            start_point_lon = vertices[1][0]
            end_point_lat   = vertices[0][1]
            end_point_lon   = vertices[0][0]        
        lines = []
        for i in range(num_lines):
            # Calculate new start and end points for each line
            new_start_lat = start_point_lat + i * lat_step_start
            new_start_lon = start_point_lon + i * lon_step_start
            new_end_lat   = end_point_lat + i * lat_step_end
            new_end_lon   = end_point_lon + i * lon_step_end
            lines.append([[new_start_lon, new_end_lon], [new_start_lat, new_end_lat]])
        return lines

    def averaged_over_lines(self, data, xlat, xlong, lines, num_points):
        averaged_data = []
        for line in lines:
            new_start_lon, new_end_lon = line[0]
            new_start_lat, new_end_lat = line[1]
            
            # Extract data along this line
            line_data = self.extract_data_along_line(
                data, xlong, xlat, 
                new_start_lat, new_start_lon, new_end_lat, new_end_lon,
                num_points
            )
            
            # Calculate the average of the extracted data
            if len(line_data) > 0:
                averaged_data.append(np.mean(line_data))
            else:
                averaged_data.append(None)
        
        return np.array(averaged_data)

    def averaged_over_lines_3D(self, data, dim, xlat, xlong, lines, num_points):
        from functools import partial
        partial_func = partial(self.averaged_over_lines, 
            xlat=xlat, xlong=xlong, 
            lines=lines, num_points=num_points)
            
        return xr.apply_ufunc(partial_func, data,
            input_core_dims=[dim],
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
#    vertices = (np.array([-95.09765264,  29.374     ]),
#         np.array([-95.50234736,  29.374     ]),
#         np.array([-95.50234736,  30.853     ]),
#         np.array([-95.09765264,  30.853     ]))
    midpoint1 = (29.374, -95.30)  
    midpoint2 = (30.853, -95.30)  
    width_km = 45  # Example width in kilometers
    
    vertices, angle = mycross_section.calculate_geographic_rectangle(midpoint1, midpoint2, width_km)
    print("Vertices of the rectangle are:", vertices)
    print('angle from east:', angle)
    da = xr.open_dataarray('../../test_data.nc')
    lines = mycross_section.calculate_lines(vertices, 20, 0)
    averages = mycross_section.averaged_over_lines_3D(da*100, ['south_north','west_east'], 
        da['XLAT'].data, da['XLONG'].data, lines, 10)
    print(averages)

