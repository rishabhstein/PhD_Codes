#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:50:27 2019

@author: rishabhstein
"""

import numpy as np

def makeF(data, name='F', num=0,):
    
    N = np.size(data)
    data = np.reshape(data,(N,),order='F')
    
    fr = open('0.org/'+name+'org')
    fw = open('0.org/'+name+str(num), 'w')
    
    for k in range(19):
        line = fr.readline()
        if k == 13:
            line = '    object      F' + str(num)+';'+'\n'
        fw.write(line)
    
    data = np.array(data,float)
    fw.write('internalField   nonuniform List<scalar>\n')  
    fw.write(str(N) + '\n(\n')
    for k in range(N):
        if data[k] == 1.0:
            fw.write(str(0.99) + '\n')
        else:
            fw.write(str(data[k]) + '\n')
    fw.write(')\n;\n\n')
    
    write = False
    while 1:
        line=fr.readline()
        if not line: break
        if line == 'boundaryField\n':
            write = True
        if write:
            fw.write(line)   

    fr.close()    
    fw.close()


def calculate_porosity(img, img_grid_info, start):
    """This function returns the average porosity 
        value of the input grids"""
    dims = len(img_grid_info);
    if dims == 3:
        tmp_x, tmp_y, tmp_z = img_grid_info;
        tmp_int_x, tmp_int_y, tmp_int_z = int(tmp_x), int(tmp_y), int(tmp_z);
        p = 0.; 
        start_x, start_y, start_z = start;

        for i in range(start_x * tmp_int_x, (start_x +1) * tmp_int_x):
            for j in range(start_y * tmp_int_y, (start_y +1) * tmp_int_y):
                for k in range(start_z * tmp_int_z, (start_z +1) * tmp_int_z):
                    p += img[i,j,k]
        return p/(tmp_x*tmp_y*tmp_z);
    elif dims == 2:
        tmp_x, tmp_y = img_grid_info;
        tmp_int_x, tmp_int_y = int(tmp_x), int(tmp_y);
        p = 0.; 
        start_x, start_y = start;
        for i in range(start_x * tmp_int_x, (start_x +1) * tmp_int_x):
            for j in range(start_y * tmp_int_y, (start_y +1) * tmp_int_y):
                p += img[i,j]
        return (p/(tmp_x*tmp_y));


def get_porosity(img, Nx = 1, Ny = 1, Nz = 1):
    """This function takes an image/stack_image as input and
    return an array of scaled porosity field for whole geometry
    for given mesh numbers in x,y,z direction respectively"""
    dims = np.shape(img);
    if len(dims) == 3:
        size_x, size_y, size_z = np.shape(img);

        axial_fraction_x = size_x/Nx; #No. of pixels in x direction of grids
        axial_fraction_y = size_y/Ny; #No. of pixels in y direction of grids
        axial_fraction_z = size_z/Nz; #No. of pixels in y direction of grids

        if ((axial_fraction_x / int(axial_fraction_x) != 1.0) or 
           (axial_fraction_y / int(axial_fraction_y) != 1.0) or 
           (axial_fraction_z / int(axial_fraction_z) != 1.0)):
            print("Please check the ratio of grids and pixels\n");
            print("Ratio of grids and pixels is: %f, %f, %f \n" %(axial_fraction_x, axial_fraction_y, axial_fraction_z));
            return 0;
        
        img_grid_info = [axial_fraction_x, axial_fraction_y, axial_fraction_z];
        a = np.zeros([Nx,Ny,Nz]);
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    loc = [i, j, k];   
                    a[i,j,k] = calculate_porosity(img, img_grid_info, loc)       
        return (a, np.sum(a)/np.size(a));
    
    elif len(dims) == 2:
        size_x, size_y = np.shape(img);        
        axial_fraction_x = size_x/Nx; #No. of pixels in x direction of grids
        axial_fraction_y = size_y/Ny; #No. of pixels in y direction of grids
        
        if ((axial_fraction_x / int(axial_fraction_x) != 1.0) or 
           (axial_fraction_y / int(axial_fraction_y) != 1.0)):
            print("Please check the ratio of grids and pixels\n");
            print("Ratio of grids and pixels is: %f, %f\n" %(axial_fraction_x, axial_fraction_y));
            return 0;
        
        img_grid_info = [axial_fraction_x, axial_fraction_y];
        a = np.zeros([Nx,Ny]);

        for i in range(Nx):
            for j in range(Ny):
                loc = [i, j];   
                a[i,j] = calculate_porosity(img, img_grid_info, loc)       
        return (a, np.sum(a)/np.size(a));



    
    
def crop_image(filename, size):
    """This function takes image as input and crop it 
    according to given size.
    Also, if the image is in 16-bit or 64 bit it will scale it to
    8-bit image"""
    pix_x, pix_y, pix_x2, pix_y2 = size;
    import cv2 as cv;
    img = cv.imread(filename, cv.IMREAD_GRAYSCALE)
    crop_img = img[pix_x: pix_x2, pix_y:pix_y2]
    if (np.amax(crop_img)) > 255:
        crop_img = crop_img/(255**2);
    elif (np.amax(crop_img)) > 1:
        crop_img = crop_img/255;
    return crop_img;
