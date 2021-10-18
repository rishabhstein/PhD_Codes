#!/usr/bin/env python3
"""
Created on Thu May 16 12:50:27 2019

@author: rishabhstein
"""

import numpy as np

def roughness(dims, f):
    '''Create a random mineral field
        dims - tuple with system dimensions
        f - list of mineral porosity and perturbation'''
    
    data  = np.ones(dims)*f[0]
    data += (2*np.random.random(dims)-1)*np.sqrt(3.)*f[1]
    
    return data

def create_porosity_roughness(dims, porosity = 1, noise = 1):
    '''Create a porosity field with dims, dimensions and
    noise'''
    data  = np.ones(dims)*porosity
    data += (2*np.random.random(dims)-1)*noise
    return data


def porosity_profile(img, axis = None, sample_type = None, void_fraction = 1):
    """
    This function calculates porosity for each slice
    and return a list of porosity values for all slices
    
    Input:
        img: Nd-array as image data
        
        axis: Along which axis the porosity profile needs to be calculated
        
        sample_type: This has 3 options
            a. 3_Phase :  Thresholding in Pores, Micropores and Grains.
            
            b. Linda_et_al :  Thresholding in Pores, Micropores and Grains but Micropores are linearly correlated to Pores and Grains.
            
            c. Core_1_Phase : Core scans(with core and extra part of image) are used with Thresholding in 1 and 0.
            
            d. Crop_1_Phase : Cropped part of scans(cropped part of sample with no extra sides) are used with thresholding 1 and 0.
    
    
    """
    
    
    img = img.copy();
    phi = [];
    
    
    if (sample_type == '3_Phase'):
        for i in img:
            n = i[i > 0].size;
            phi.append(100*np.sum(i[i>0])/n/255);
            
    elif (sample_type == 'Linda_et_al'):
        #This part only calculate for one slice
        i = img;
        n = i[i > 0].size
        phi_m = np.sum(i[i==255]);
        phi_g = np.sum(i[i==1]);
        tmp = i[i < 255];
        phi_micro =  np.sum(tmp[tmp>1])
        
        porosity = (100*(phi_m + phi_g + void_fraction*phi_micro)/n/255)
        phi.append(porosity)
            
            
    elif (sample_type == 'Core_2_Phase'):
        print("This part is still need to implement");
        
    elif (sample_type == 'Crop_1_Phase'):
        for i in img:
            n = i.size;
            phi.append(100*np.sum(i)/n/255);
            
    return phi;



def n_weighted_moment(values, weights, n):
    values = np.array(values)
    if weights == 1:
        weights = np.ones(values.shape)
        
    assert n>0 & (values.shape == weights.shape)
    w_avg = np.average(values, weights = weights)
    w_var = np.sum(weights * (values - w_avg)**2)/np.sum(weights)

    if n==1:
        return w_avg
    elif n==2:
        return w_var
    else:
        w_std = np.sqrt(w_var)
        return np.sum(weights * ((values - w_avg)/w_std)**n)/np.sum(weights)



def Segmentation_Linda_et_al(img, segmentation_bound):
    
    """This method of 3-phase thresholding is taken from Linda et al's article "Experimental Characterization of Porosity Structure 
    and Transport Property Changes in Limestone Undergoing Different Dissolution Regimes"
    
    Porosity_pore = 1;
    Porosity_grain = 0;
    Porosity_micro_pore =  phi =  Gm - Gx/(Gm - Gv) 
    
    Input: 
        img: 2D array of an image
        
        segmentation_bound: array of [1x2] for lower and upper limit of threshold for macropores and grains respectively
        
    Return:
        return the 3phase segmented image
    
    """
    
    dims = img.shape;  
    Th_V, Th_G = segmentation_bound
    
    
    #Convert image to (pixel_I <= 1)
#     if (np.max(img) > 0 & np.max(img < 255)): #Test for 8-bit image
#         img = img/255;

    if (np.max(img) > 255 & np.max(img) < 255**2): #Convert 16 bit to 8 bit image
        img = (img/255).astype(np.uint8);
    elif (np.max(img) > 255**2):                   #Convert 32 bit to 8 bit image
        img = (img/255**2).astype(np.uint8);      
    
    
    """Needs to be implemented using histographic data"""
    Gm = Th_G;   #Threshold_grain
    Gv = Th_V;   #Threshold_pores

    slice2 = img.copy();
    phi = [];
    
    slice2[img >= Th_G] = 1; #Grain
    slice2[(img <= Th_V) & (img != 0)] = 255; #Pore
    
    slice2[(img == 0)] = 0; #Background
    
    tmp = img[(img < Th_G) & (img > Th_V)];
    slice2[(img < Th_G) & (img > Th_V)] = 255*((Gm - tmp)/(Gm - Gv)) #Micro_Pores
    
    #     Old but Slow method    
    #        for i in range(dims[0]):
    #            for j in range(dims[1]):
    #                if slice[i,j] <= Th_V :
    #                    phi.append(1)
    #                elif (slice[i,j] > Th_V) & (slice[i,j] < Th_G):
    #                    phi.append((Gm - slice[i,j])/(Gm-Gv))
    #                elif slice[i,j] >= Th_G & slice[i,j] != 255 :
    #                    phi.append(0.01)
    #                elif slice[i,j] == 255:
    #                    phi.append(0)

    #        phi = np.reshape(phi*255,[dims[0],dims[1]])



    return slice2.astype(np.uint8);