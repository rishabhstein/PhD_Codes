#!/usr/bin/env python3
"""
Created on Thu April 17 12:50:27 2020

@author: rishabhstein
"""

import numpy as np

def seaborn_plot_normalise(image_obj, norm_xy = 'Y', index_subplot=0):    
    line = image_obj.get_lines()[index_subplot]
    xd = line.get_xdata()
    yd = line.get_ydata()

    xd2 = (xd - xd.min(0)) / xd.ptp(0)
    yd2 = (yd - yd.min(0)) / yd.ptp(0)
    
    
    # #normalize points
    if(norm_xy == 'Y'):
        return(xd, yd2)
    elif(norm_xy == 'X'):
        return(xd2, yd)
    elif(norm_xy == 'XY'):
        return(xd, yd2)
    