#!/usr/bin/env python3
"""
Created on Thu May 17 12:50:27 2020

@author: rishabhstein
"""

import numpy as np
from pathlib import Path
import vtk

def load_ARM_data(path, prefix = '', time = 'First', file_type = '.out'):
    """ This function reads the output data of network model
    Currently following files can be provided
     For '.out' format : d_prefix.out, l_prefix.out, q_prefix.out 
     For '.csv' format : data_prefix.csv which contains diameters and flow rate can be provided.
     
    INPUT: 
           path  : Path of the dir containing the output files
           prefix: Prefix for files if any as string
           
           time  : string 'First' of 'Last'. 
                   Still need to be implemented for custom time steps
           
           file_type: '.out' for .out files and '.csv' for csv file fro PARAVIEW 
             
    OUTPUT:
            Dictionary of Network with the properties from file.
    
    """

    
    network = {}
    if (file_type == '.out'):
        # Reading diameter files 
        file = Path(path.resolve(), 'd' + prefix +'.out')
        dia = [];  count_steps = [];
        with open(file,'r') as f:
            data = f.readlines()
            count_steps.extend([i+1 for i in range(len(data)) if data[i][0]=='#'])
            step_len = count_steps[1] - count_steps[0] - 1 

            if(time == 'First'):
                for i in range(count_steps[0], step_len + count_steps[0]):
                    dia.extend(data[i].strip('\t').split())
            elif(time == 'Second'):
                for i in range(count_steps[1], step_len + count_steps[1]):
                    dia.extend(data[i].strip('\t').split()) 
            elif(time == 'Last'):
                for i in range(count_steps[-1], step_len + count_steps[-1]-1):
                    dia.extend(data[i].strip('\t').split())
           
        network['diameter'] = np.vstack([float(i) for i in dia])

#         file= Path(path.resolve(), 'l' + prefix +'.out')
#         lth = []; count_steps = [];
#         with open(file,'r') as f:
#             data = f.readlines()
#             count_steps.extend([i+1 for i in range(len(data)) if data[i][0]=='#'])
#             step_len = count_steps[1] - count_steps[0] - 1 

#             if(time == 'First'):
#                 for i in range(count_steps[0], step_len + count_steps[0]):
#                     lth.extend(data[i].strip('\t').split())
#             elif(time == 'Last'):
#                 for i in range(count_steps[-1], step_len + count_steps[-1]-1):
#                     lth.extend(data[i].strip('\t').split())

#        network['length'] = np.vstack([float(i) for i in lth])

        file= Path(path.resolve(), 'q' + prefix +'.out')
        q = []; count_steps = [];
        with open(file,'r') as f:
            data = f.readlines()
            count_steps.extend([i+1 for i in range(len(data)) if data[i][0]=='#'])
            step_len = count_steps[1] - count_steps[0] - 1 

            if(time == 'First'):
                for i in range(count_steps[0], step_len + count_steps[0]):
                    q.extend(data[i].strip('\t').split())
            if(time == 'Second'):
                for i in range(count_steps[1], step_len + count_steps[1]):
                    q.extend(data[i].strip('\t').split())            
            elif(time == 'Last'):
                for i in range(count_steps[-1], step_len + count_steps[-1]-1):
                    q.extend(data[i].strip('\t').split())

        network['flow_rate'] = np.vstack([float(i) for i in q])
        
        #""" Trimming the zero diameter pores and related flow_rates if any """ 

        trim_pores = np.where(np.any(network['diameter'] == 0, axis = 1))[0]
        dia_after_trimming = np.delete(network['diameter'], trim_pores, axis=0)
        flow_after_trimming = np.delete(network['flow_rate'], trim_pores, axis=0)

        network['diameter'] = np.vstack(dia_after_trimming)
        network['flow_rate'] =np.vstack(flow_after_trimming)


    elif (file_type == '.csv'):
        import pandas as pd
        file = Path(path.resolve(), 'data'+ prefix +'.csv')
        net_data_csv = pd.read_csv(filepath_or_buffer=file,
                                    usecols=['Diameter', 'Flow_Rate'])
        
        network['diameter'] = np.vstack(net_data_csv['Diameter'])
        network['flow_rate'] =np.vstack(net_data_csv['Flow_Rate'])
        
        # trimming the zero diameter pores and related flow_rates if any
        trim_pores = np.where(np.any(network['diameter'] == 0, axis = 1))[0]
        dia_after_trimming = np.delete(network['diameter'], trim_pores, axis=0)
        flow_after_trimming = np.delete(network['flow_rate'], trim_pores, axis=0)
        
        network['diameter'] = np.vstack(dia_after_trimming)
        network['flow_rate'] =np.vstack(flow_after_trimming)

    return network


def load_ARM_data_fromVTK(file_name):
    
    """This function reads the Network model output VTK files
        ==========================
        "Without Pore-Merging"
        ==========================
        and returns a Dictionary containing points and scalars in the file
        
        INPUT:
            VTK file saved during simulations
           
        OUTPUT:
            A Dictionary containing cordinates of network points and all scalar fields
               
    """
    #Dictionary where the data of VTK file
    scalars = {}
    
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(file_name)
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    points = np.array(data.GetPoints().GetData())
    pointdata = data.GetPointData()
    celldata = data.GetCellData()

    #Extracting Cell bounds or pore nodes
    CellArray = data.GetCells()
    pores = CellArray.GetData()
    no_of_pores = CellArray.GetNumberOfCells()
    
    pores_array = np.array(pores)
    list_of_pore_connection = []

    for i in range(no_of_pores):
        list_of_pore_connection.append([pores_array[j] for j in range(i*3, i*3+3)])
    scalars['Cell_Nodes'] = list_of_pore_connection;
    
    #No of properties saved in VTK files
    no_of_fields = reader.GetNumberOfScalarsInFile() 

    scalars['Points'] = points;
    for i in range(no_of_fields): 
        name_of_field = reader.GetScalarsNameInFile(i)
        scalars[name_of_field] = np.array(celldata.GetArray(name_of_field))
           
                
    return scalars



def Diameter_based_tip_position(LIST_OF_VTK_FILES, beta_d = 2):
    """
    This function takes input of network.vtk files and extract 
    tip position of wormhole based on evolved diameters
    
    d_evolved >= beta_d * d_initial
    
    Input:
        1. List of VTK files including path
        2. beta_d ratio of d_evolved/d_initial
        
    Output:
        List of tip_position
        
    """
    tip_position_d = [];
    data = load_ARM_data_fromVTK(LIST_OF_VTK_FILES[0]); #First file 'network_0.vtk'
    
    for i in LIST_OF_VTK_FILES[1:]:
        data2 = load_ARM_data_fromVTK(i) #iterating over all other vtk files
        
        dissolved_media_d = (data2['Diameter'] > beta_d*data['Diameter'])
        tmp_d = [i for i, x in enumerate(dissolved_media_d) if x]
        
        dissolved_pore_connections = [data2['Cell_Nodes'][j] for j in tmp_d]

        list_centroid_pores = [];
        for k in range(len(dissolved_pore_connections)):
            node_1 = dissolved_pore_connections[k][1] #Node-1 connected to Pore
            node_2 = dissolved_pore_connections[k][2] #Node-2 connected to Pore
            list_centroid_pores.append((data2['Points'][node_1][1] + data2['Points'][node_2][1])/2.0)

        if len(list_centroid_pores) != 0 :
            tip_position_d.append(np.max(list_centroid_pores))
    return np.sort(tip_position_d)
        
    
    
def Concentration_based_tip_position(LIST_OF_VTK_FILES, beta_c = 0.1):
    """
    This function takes input of network_%d.vtk files and extract 
    tip position of wormhole based on evolved Concentration
    
    c_evolved < beta_c
    
    Input:
        1. List of VTK files including path
        2. beta_d ratio of d_evolved/d_initial
        
    Output:
        List of tip_positions
        
    """
    
    tip_position_c = [];
    data = load_ARM_data_fromVTK(LIST_OF_VTK_FILES[0]);
    
    for i in LIST_OF_VTK_FILES[1:]:
        data2 = load_ARM_data_fromVTK(i)
        dissolved_media_c = (data2['Concentration'] > beta_c)
        tmp_c = [i for i, x in enumerate(dissolved_media_c) if x]
        
        dissolved_pore_connections = [data2['Cell_Nodes'][j] for j in tmp_c]

        list_centroid_pores = [];
        for k in range(len(dissolved_pore_connections)):
            node_1 = dissolved_pore_connections[k][1] #Node-1 connected to Pore
            node_2 = dissolved_pore_connections[k][2] #Node-2 connected to Pore
            list_centroid_pores.append((data2['Points'][node_1][1] + data2['Points'][node_2][1])/2.0)

        if len(list_centroid_pores) != 0 :
            tip_position_c.append(np.max(list_centroid_pores))
    return np.sort(tip_position_c)


def ARM_VTKtoCSV(path_to_fileIn, path_to_fileOut, fname):
    import vtk
    import csv
    """This function converts the VTK Polydata file to 
    CSV files. 
    
    Currently save only Diameter and Flow_Rate
    
    
    Input: 
        path to file.vtk
        path to file.csv
        file name
    Output:
        Save file.csv to path_to_fileOut/file.csv
    
    """
    fileIn  = Path(path_to_fileIn.resolve(), fname)
    fileOut  = Path(path_to_fileOut.resolve(), fname.replace('.vtk', '.csv'))
    reader = vtk.vtkDataSetReader()
    reader.SetFileName(str(fileIn))
    reader.ReadAllScalarsOn()
    reader.Update()
    
    cell_obj = reader.GetOutput()

    points = cell_obj.GetPoints()

    dia = cell_obj.GetCellData().GetArray('Diameter')
    flow = cell_obj.GetCellData().GetArray('Flow_Rate')

    table = vtk.vtkDataObjectToTable()
    table.SetInputData(cell_obj)
    table.Update()
    table.GetOutput().AddColumn(dia)
    table.GetOutput().AddColumn(flow)
    table.Update()

    writer = vtk.vtkDelimitedTextWriter()
    writer.SetInputConnection(table.GetOutputPort())
    writer.SetFileName(str(fileOut))
    writer.Update()
    writer.Write()    