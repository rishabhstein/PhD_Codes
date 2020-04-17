#!/usr/bin/env python3
"""
Created on Thu May 17 12:50:27 2020

@author: rishabhstein
"""

import numpy as np
from pathlib import Path

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

            if(time == 'first'):
                for i in range(count_steps[0], step_len + count_steps[0]):
                    dia.extend(data[i].strip('\t').split())
            elif(time == 'last'):
                for i in range(count_steps[-1], step_len + count_steps[-1]-1):
                    dia.extend(data[i].strip('\t').split())

        network['diameter'] = np.vstack([float(i) for i in dia])

        file= Path(path.resolve(), 'l' + prefix +'.out')
        lth = []; count_steps = [];
        with open(file,'r') as f:
            data = f.readlines()
            count_steps.extend([i+1 for i in range(len(data)) if data[i][0]=='#'])
            step_len = count_steps[1] - count_steps[0] - 1 

            if(time == 'first'):
                for i in range(count_steps[0], step_len + count_steps[0]):
                    lth.extend(data[i].strip('\t').split())
            elif(time == 'last'):
                for i in range(count_steps[-1], step_len + count_steps[-1]-1):
                    lth.extend(data[i].strip('\t').split())

        network['length'] = np.vstack([float(i) for i in lth])

        file= Path(path.resolve(), 'q' + prefix +'.out')
        q = []; count_steps = [];
        with open(file,'r') as f:
            data = f.readlines()
            count_steps.extend([i+1 for i in range(len(data)) if data[i][0]=='#'])
            step_len = count_steps[1] - count_steps[0] - 1 

            if(time == 'first'):
                for i in range(count_steps[0], step_len + count_steps[0]):
                    q.extend(data[i].strip('\t').split())
            elif(time == 'last'):
                for i in range(count_steps[-1], step_len + count_steps[-1]-1):
                    q.extend(data[i].strip('\t').split())

        network['flow_rate'] = np.vstack([float(i) for i in q])
    
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
