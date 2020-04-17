#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 12:50:27 2019

@author: rishabhstein
"""

import numpy as np
import os
import cv2
import math

def makeF(data, case = None, name='F', num=0,):
    '''
    This function creates the F0 file for porousFoam
    Input:
        data: nd-array as porosity field
        case: '2D' or '3D'
    '''
    
    
    N = np.size(data)
    data = np.reshape(data,(N,),order='F')
    
    if case == '2D':
        os.system('cd 0.org && cp Forg2Dcyclic Forg')
    elif case == '3D':
        os.system('cd 0.org && cp Forg3DzeroGrad Forg')    
        
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
    """
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    If Input Image is binary and pores are black then
    porosity = 1 - sum of intensity of white pixels
    If pores are white then
    porosity = sum of intensity of white pixels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    This function takes an image/stack_image as input and
    return an array of scaled porosity field for whole geometry
    for given mesh numbers in x,y,z direction respectively
    """
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
        pore_matrix = np.zeros([Nx,Ny,Nz]);
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    loc = [i, j, k];   
                    pore_matrix[i,j,k] = calculate_porosity(img, img_grid_info, loc)       
        return (pore_matrix, np.sum(pore_matrix)/np.size(pore_matrix));
    
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
        pore_matrix = np.zeros([Nx,Ny]);

        for i in range(Nx):
            for j in range(Ny):
                loc = [i, j];   
                pore_matrix[i,j] = calculate_porosity(img, img_grid_info, loc)       
        return (pore_matrix, np.sum(pore_matrix)/np.size(pore_matrix));



    
    
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




def Aligning_images(img, ref_point, shift = None):
    """
    =============Still Untested=============
    
    
    Input:
        image: Ndarray
            It could be a single image or a 3D array of images.
        
        ref_point: Ndarray
            It is an array of points (rows = 3, cols = 2) for Affine transformation.
            example = [[x0,y0], [x1,y1], [x2,y2]]
        
        shift:
            list [x_shift, y_shift] of shift in x and y axis from reference image.
            If slices are more than one than shift will be divided into slices.
    Return:
        It returns the aligned image or image stack   
    """
    
    img = img.copy();
    dims = img.shape;
    
    if len(dims)==3:
            no_of_slices, rows, cols = img.shape;
    elif len(dims)==2:
            rows, cols = img.shape;
            no_of_slices = 1;
            
    if (np.size(ref_point) != 6):
        print("Reference point array's dimensions are not correct \n");
        print('It should be in the form of [[x0,y0], [x1,y1], [x2,y2]]');
        return 0;
    
    
    diffx, diffy = shift
    
    thetax = math.asin(diffx/no_of_slices)
    thetay = math.asin(diffy/no_of_slices)

    transform_matrix = np.ndarray((no_of_slices,2))
    
    #reference points of first image
    pt1 = np.float32(ref_point)    
    
    align_img = [];
 
    for i in range(no_of_slices):
        transform_matrix[i][0] = int(i * math.sin(thetax)) #shift in X
        transform_matrix[i][1] = int(i * math.sin(thetay)) #shift in Y   
        
        if len(dims) == 2:
            img_2_transform = img;
        else:
            img_2_transform = img[i,:,:];
            
        #Points for affine transformation
        pt2 = np.float32([[pt1[0,0]+transform_matrix[i][0],pt1[0,1]+transform_matrix[i][1]],
                          [pt1[1,0]+transform_matrix[i][0],pt1[1,1]+transform_matrix[i][1]],
                          [pt1[2,0]+transform_matrix[i][0],pt1[2,1]+transform_matrix[i][1]]])

        matrix = cv2.getAffineTransform(pt1, pt2)
        align_img.append(cv2.warpAffine(img_2_transform, matrix, (cols, rows)));
    
    align_array = np.reshape(align_img, dims);

    return align_array
