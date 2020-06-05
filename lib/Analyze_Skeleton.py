from skan import skeleton_to_csgraph, summarize
import networkx as nx
import numpy as np
import vtk 
from vtk.util.colors import banana, plum 
import copy
import itertools



class Analyze_Skeleton:
    """
    This class Analyze the skeletonized image using "Skan" and "NetworkX".
    Following properties of Wormhole skeleton are calculated:
    
    1. Wormhole Length
    2. Tortuosity
    3. Wastefulness
    4. Horton number for branching (NOT IMPLEMENTED YET)
    5. Fractal dimension           (NOT IMPLEMENTED YET)
    
    """
    def __init__(self, skeletonized_image, suffix = 'Sample'):
        self.skeletonized_image = skeletonized_image
        self.node_connection = []
        self.node_cordinates = []
        self.G = nx.Graph();
        self.inlet_node = 0;
        self.outlet_node = 0;
        self.pos = {}
        self.length_of_wormhole = 0;
        self.suffix = suffix
        self.dominant_wormhole = []
        
        
    def Skeletonized_image_to_NetworkX_Graph(self):
        
        #This function convert skeletonized image to NetworkX graph
        
        
        dims = np.shape(self.skeletonized_image)
        self.node_connection, self.node_cordinates,  deg = skeleton_to_csgraph(self.skeletonized_image)
        self.G = nx.from_scipy_sparse_matrix(self.node_connection)

        #converting 2D graph to 3D
        if len(dims) == 2: 
            node_cordinates_copy = np.empty(shape = (len(self.node_cordinates), 3))
            for i in range(len(self.node_cordinates)):
                node_cordinates_copy[i] = [self.node_cordinates[i,0], self.node_cordinates[i,1], 0]
      
            self.node_cordinates = node_cordinates_copy
        
        # Setting the position of nodes in the NetworkX Graph for 2D or 3D    
        self.pos = {i: (self.node_cordinates[i,0], self.node_cordinates[i,1], self.node_cordinates[i,2]) for i in range(len(self.node_cordinates))}
           
        for n, p in self.pos.items():
            self.G.nodes[n]['pos'] = p

            
    def Remove_isolated_nodes(self):
        G_copy = self.G
        isolated_nodes = list(nx.isolates(G_copy))
        G_copy.remove_nodes_from(isolated_nodes)
        self.G = G_copy            
        return isolated_nodes
            
    def Nx_to_vtk(self, isolated_nodes = [], addtional_suffix = ''): 
    # set node positions 
        np={} 
       
        node_pos = self.pos
        for n in self.G.nodes(): 
            try: 
                np[n]=node_pos[n] 
            except KeyError: 
                raise networkx.NetworkXError

          # Generate the polyline for the spline. 
        points = vtk.vtkPoints() 
        edgeData = vtk.vtkPolyData() 

          # Edges 

        lines = vtk.vtkCellArray()          
            
            
        for n in self.G.nodes():
            (x,y,z) = node_pos[n]
            points.InsertPoint(n,x,y,z)
     
        #Filling zeros at deleted nodes
        try:
            for i in isolated_nodes:
                points.InsertPoint(int(i),0,0,0)
        except:
            pass
        
        tmp_u = []; tmp_v = [];
        for e in self.G.edges(): 
            u=e[0] 
            v=e[1] 
            lines.InsertNextCell(2)  
            lines.InsertCellPoint(u) 
            lines.InsertCellPoint(v)



        edgeData.SetPoints(points) 
        edgeData.SetLines(lines)

        writer = vtk.vtkXMLPolyDataWriter();
        writer.SetFileName(str(self.suffix)+ addtional_suffix +".vtp");

        writer.SetInputData(edgeData)
        writer.Write()

    def Delete_Loops_from_Skeleton(self):
        #It's better to find the dominant wormhole first
        self.Find_Dominant_Wormhole()
        try:
            if_main_wormhole_path = False
            
            cycles = nx.find_cycle(self.G)
            junction_nodes = [];
            loop_source_sync = [] #make a copy
            for c in cycles:
                u = c[0]; v = c[1];
                if len(list(self.G.neighbors(u))) > 2 :
                    junction_nodes.append(u)
                if u in self.dominant_wormhole:
                    if_main_wormhole_path = True

            shortes_path = []

            if(if_main_wormhole_path == False):
                if (len(junction_nodes) > 2):
                    for i,j in itertools.combinations(junction_nodes, 2):
                        shortes_path = nx.shortest_path(self.G,  i, j)
                        if len([w for w in junction_nodes if w in shortes_path]) > 2:
                            loop_source_sync.extend([i,j])
                else:
                    loop_source_sync = junction_nodes

            elif (if_main_wormhole_path == True):
                    print("This cycle is in main wormhole path")
                    if (len(junction_nodes) > 2):
                        for i,j in itertools.combinations(junction_nodes, 2):
                            if i in self.dominant_wormhole and j in self.dominant_wormhole:
                                shortes_path = nx.shortest_path(self.G,  i, j)
                                if len([w for w in junction_nodes if w in shortes_path]) > 2:
                                    loop_source_sync.extend([i,j])
                                    break;
                                else:
                                    loop_source_sync = [jn for jn in junction_nodes if jn in self.dominant_wormhole]

                    else:
                        loop_source_sync = junction_nodes



      #      print(loop_source_sync)
            loop_paths = list(nx.all_simple_paths(self.G, loop_source_sync[0], loop_source_sync[1]))   

            path_length = np.sort(list(len(i) for i in loop_paths))
            tmp_e2 = [];
            for i in loop_paths:
                if len(i) > path_length[0]:
                    e1 = loop_source_sync[1]
                    e2 = i[-2]
                    if e2 not in tmp_e2 and e2 not in nx.shortest_path(self.G, loop_source_sync[0], loop_source_sync[1]):
                        self.G.remove_edge(e1,e2)
                        print("Removed edges are: ", e1, e2)
                        print("from the cycle: ", i)
                        tmp_e2.append(e2)

        except:
            print('No cycle found.')

                    
                
                
    def Find_Dominant_Wormhole(self):
        
        
        #This function is works best for 3D images, because of uncertainity in choosing
        # Inlet and Outlet node in 2D
        #No. of slices should be in Z-Axis (i.e First axis here)
        
        
        inletz, _, _ = np.argmin(self.node_cordinates, axis=0)
        outletz, _, _ = np.argmax(self.node_cordinates, axis=0)
        
        #NetworkX starts number of nodes from 1 so
        inletz = inletz+1

        self.inlet_node = inletz
        self.outlet_node = outletz
        
        self.dominant_wormhole = nx.shortest_path(self.G, inletz,  outletz)
        self.length_of_wormhole = nx.shortest_path_length(self.G, inletz,  outletz, weight='weight')
        return self.length_of_wormhole, self.dominant_wormhole
        
    def Tortuosity(self):
        return self.length_of_wormhole/self.node_cordinates[self.outlet_node][0]        
        
    def Wastefulness(self):
        #Using sum of sparce matrix but every connection is doubled so
        #divide by half
        sum_of_all_branches_including_wormhole = self.node_connection.sum()/2
        return sum_of_all_branches_including_wormhole/self.node_cordinates[self.outlet_node][0]
        
    def Horton_number(self):
        print("Horton needs to be implimented")
        
    def Fractal_dimension(self):
        print("Thumbprint of God is yet to be implemented")
    