from skimage import measure, morphology, segmentation
from skimage.morphology import reconstruction, skeletonize_3d
import scipy.ndimage.morphology as morph
from scipy import ndimage
from pylab import *

from tifffile import imsave as tifsave
from skimage import util 

import sys
sys.path.append('../')
import lib.Analyze_Skeleton as AS

#Reload modules
import importlib
importlib.reload(AS)

class Wormhole_from_tomography:
    
    
    def __init__(self, wormhole_from_differencing_stack, suffix, voxel_size, downscaling_factor):
        self.wormhole_from_differencing_stack = wormhole_from_differencing_stack;
        self.filled_wormhole = [];
        self.Volume_of_wormhole = 0;
        self.Surface_area_of_wormhole = 0;
        
        self.info = {"File_name": suffix,
                    "voxel_size": voxel_size,
                    "downscaling_factor": downscaling_factor,
                    "Volume_in_voxel":self.Volume_of_wormhole,
                    "Volume in cm3": 0,
                    "Surface_area_in_voxel": self.Surface_area_of_wormhole,
                    "Surface_area_in_cm2": 0
                    }
                    

    def Connected_components_extraction(self, th_value, point):
                                      
        z, x, y = point;
        
        #Segmentation
        mask = self.wormhole_from_differencing_stack >= th_value
        self.filled_wormhole = segmentation.flood(mask, (z, y, x), connectivity = 1)

        #Unfilled wormhole's volume
        self.Volume_of_wormhole = self.filled_wormhole.sum()
        
        return self.Volume_of_wormhole, self.filled_wormhole
    
    
    def Connected_component_extraction_of_binary_wormholes(self, mask_img, point):
        z, x, y = point;
        
        self.filled_wormhole = segmentation.flood(mask_img, (z, y, x), connectivity = 1)

        #Unfilled wormhole's volume
        self.Volume_of_wormhole = self.filled_wormhole.sum()
        
        return self.Volume_of_wormhole, self.filled_wormhole
    

    def Fill_holes_in_wormhole(self, if_update_volume = True):
        unfilled_wormhole_copy = self.filled_wormhole
        self.filled_wormhole = morph.binary_fill_holes(unfilled_wormhole_copy)
                                        
        if if_update_volume:
            self.Volume_of_wormhole = self.filled_wormhole.sum()
            self.info['Volume_in_voxel'] = self.Volume_of_wormhole
        
        return self.Volume_of_wormhole, self.filled_wormhole
        
        
    def Fill_holes_in_wormhole_external_image(self, sample_image, if_update_class = False):
        filled_wormhole = morph.binary_fill_holes(sample_image)
        Volume_of_wormhole = filled_wormhole.sum()
        
        if if_update_class:
            self.filled_wormhole = filled_wormhole
            self.Volume_of_wormhole = Volume_of_wormhole
            self.info['Volume_in_voxel'] = self.Volume_of_wormhole
            
        return Volume_of_wormhole, filled_wormhole

                                        
                                        
    def Surface_area_Legland_et_al(self):
        """
        This function calculates the surface area of 3D object in binary images.
        """
        
        bimg = self.filled_wormhole;
        
        assert len(bimg.shape) == 3, "bimg should be a volume"
        vol = bimg.sum()

        # weights for each direction
        w = zeros(13)
        w[0:3] = 0.04577789120476 * 2
        w[3:9] = 0.03698062787608 * 2
        w[9:] =  0.03519563978232 * 2

        # number of connected components in different directions
        nc = zeros(13)

        # primary perpendicular planes
        nc[0] = vol - logical_and(bimg[:,:-1,:], bimg[:,1:,:]).sum()

        nc[1] = vol - logical_and(bimg[:-1,:,:], bimg[1:,:,:]).sum()
        nc[2] = vol - logical_and(bimg[:,:,:-1], bimg[:,:,1:]).sum()

        d1 = 1; d2 = 1; d3 = 1;

        # plane diagonals
        nc[3] = vol - logical_and(bimg[1:,:-1,:], bimg[:-1,1:,:]).sum()
        nc[4] = vol - logical_and(bimg[:-1,:-1,:], bimg[1:,1:,:]).sum()
        nc[5] = vol - logical_and(bimg[:,1:,:-1], bimg[:,:-1,1:]).sum()
        nc[6] = vol - logical_and(bimg[:,:-1,:-1], bimg[:,1:,1:]).sum()
        nc[7] = vol - logical_and(bimg[1:,:,:-1], bimg[:-1,:,1:]).sum()
        nc[8] = vol - logical_and(bimg[:-1,:,:-1], bimg[1:,:,1:]).sum()

        nc[9] = vol - logical_and(bimg[:-1,:-1,:-1], bimg[1:,1:,1:]).sum()
        nc[10] = vol - logical_and(bimg[1:,:-1,:-1], bimg[:-1,1:,1:]).sum()
        nc[11] = vol - logical_and(bimg[:-1,1:,:-1], bimg[1:,:-1,1:]).sum()
        nc[12] = vol - logical_and(bimg[1:,1:,:-1], bimg[:-1,:-1,1:]).sum()

        d12 = hypot(d1, d2); d13 = hypot(d1, d3); d23 = hypot(d2,d3);
        d123 = sqrt(d1**2 + d2**2 + d3**2)

        s = 4*d1*d2*d3 * ( nc[0]*w[0]/d1 + nc[1]*w[1]/d2 + nc[2]*w[2]/d3 + \
                                          (nc[3] + nc[4])*w[3]/d12 + (nc[5] + nc[6])*w[5]/d13 + (nc[7] + nc[8])*w[7]/d23 + \
                                          (nc[9] + nc[10] + nc[11] + nc[12])*w[9]/d123 )

        

        self.Surface_area_of_wormhole = s
        self.info['Surface_area_in_voxel'] = s;
        
        return self.Surface_area_of_wormhole
    
    
    
    def Volume_and_Surface_area_Cm2(self):  
        
        # Calculating Volume and surface area in cm scale
        
        vol_tmp = self.Volume_of_wormhole*(self.info['downscaling_factor']*self.info['voxel_size'])**3/1000**3/10**3
        sa_tmp = self.Surface_area_of_wormhole*(self.info['downscaling_factor']*self.info['voxel_size'])**2/1000**2/100**2;
 
        self.info["Volume in cm3"] = vol_tmp;
        self.info["Surface_area_in_cm2"]  = sa_tmp
        return vol_tmp, sa_tmp
    
    
    def Save_extracted_wormhole(self, save_as_bool = True):

        bool_img = self.filled_wormhole
        
        if save_as_bool:
            tifsave(self.info['File_name'], bool_img, dtype = uint8)
        else:
            img_uint8 =  np.uint8(bool_img * 255)
            tifsave(self.info['File_name'], img_uint8, dtype = uint8)
    

    def Apply_median_filter(self, filter_size = 2):
        self.filled_wormhole = ndimage.median_filter(self.filled_wormhole, size = filter_size)
        
    def Skeletonise_the_wormhole(self, outlet_node = 0, 
                                 inlet_node = 'NULL',
                                 if_apply_median_filter = True,
                                 if_save_skeleton_tif = False):
        
        if if_apply_median_filter:
            filterd_image = ndimage.median_filter(self.filled_wormhole, size = 2)
            skeleton = skeletonize_3d(filterd_image)
        else:
            skeleton = skeletonize_3d(self.filled_wormhole)
            
        if if_save_skeleton_tif:
            skeleton_file_name = self.info['File_name'].strip('.tif') + '_skeleton.tif'
            tifsave(skeleton_file_name, skeleton, dtype = uint8)
            
        
        load_skeleton = AS.Analyze_Skeleton(skeleton, suffix = self.info['File_name'])
        load_skeleton.Skeletonized_image_to_NetworkX_Graph()
        isolated_node = load_skeleton.Remove_isolated_nodes()
        load_skeleton.Nx_to_vtk();
        
        length_of_wormhole, wormhole_network_graph = load_skeleton.Find_Dominant_Wormhole(outlet_node = outlet_node, 
                                                                                          given_inlet_node = inlet_node)
        
        self.info['Length_of_wormhole'] = length_of_wormhole 
        self.info['Length_of_wormhole_cm'] = length_of_wormhole * self.info['voxel_size'] * 1e-4
        self.info['Tortuosity'] = load_skeleton.Tortuosity()
        self.info['Wastefullness'] = load_skeleton.Wastefulness()
        

        return wormhole_network_graph
    
    
    

######## Not a part of CLASS###########
    
    
    
def Find_linear_region_of_curve(data, window_of_linearity = 10,                                   
                                closeness_to_zero_derv_one = 2, 
                                closeness_to_zero_derv_two = 2,
                                step = 10):
    
    derV1 = [-(data[i] - data[i+1])/step for i in range(len(data)-1)] #First derivative
    der2V = [(derV1[i] - derV1[i+1])/step for i in range(len(derV1)-1)] #Second derivative
    
    roi = [derV1.index(k) for k in derV1 if abs(k) < closeness_to_zero_derv_one]
    print("Region of interest starts from " + str(roi[0]))
    
    Min_Threshold = 0;
    
    for i in roi: #range(len(der2V)-window_of_linearity):
        if der2V[i] is not nan:
            if (abs(der2V[i]) < closeness_to_zero_derv_two):
                base_value = der2V[i]
                count = 0;
                for j in range(1, window_of_linearity):
                    if ((abs(der2V[i+j]) < closeness_to_zero_derv_two) 
                        and (abs(der2V[i+j] - base_value) < closeness_to_zero_derv_two)):
                        
                        count += 1;
                        if count == window_of_linearity - 1:
                            print("Linear section found")
                            Min_Threshold = i
                            return Min_Threshold;
                    
        else:
            print("Nan found,  check data for error")
                    
    return 0;


def Find_threshold_with_10_percent_volume_decrement(data, Min_Threshold_index):
    list_volume = []
    list_volume.append(Min_Threshold_index)
    def count_iteration(count):
        if count == 0:
            return True;
        else: 
            return False;
    
    count_10 = 0; count_20 = 0; count_30 = 0;    
    for i in range(Min_Threshold_index, len(data)):
        percent_change_volume = 100*(1 - data[i]/data[Min_Threshold_index])
       # print(percent_change_volume)

        if percent_change_volume//10 == 1 and count_10 == 0 :
            count_10 = 1;
            list_volume.append(i)
        elif percent_change_volume//20 == 1 and count_20 == 0 :
            count_20 = 1;
            list_volume.append(i)
        elif percent_change_volume//30 == 1 and count_30 == 0 :
            count_30 = 1;
            list_volume.append(i)

            
    return list_volume


def Find_tip_position_from_binary_tiff_stack_by_maximum_z(img):
    """
    INPUT:
        Binary image of extracted wormhole

    RETURN:
        The tip position


        1.    This function gives a rough idea of extracting tip in the slice by finiding the closest slice to OUTPUT
        or Max(z). 
        2. Assumption is that there is only one tip at that slice because it takes average of white pixel in that slice

    """
    pz, px, py = np.shape(img)
    if np.max(img) == 255:
        img_max = 255
    else:
        img_max = 1

    for i in range(pz-1, 0, -1):
        if img[i].any() == img_max:
            tip_z = i
            tip_y, tip_x = np.where(img[i] == img_max)
            tip_x = int(np.mean(tip_x))
            tip_y = int(np.mean(tip_y))

            # px - tip_x is used to match the origin for paraview
            return tip_x, py-tip_y, tip_z 
    return 0
    