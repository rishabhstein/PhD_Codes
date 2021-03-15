from skan import skeleton_to_csgraph, summarize, Skeleton
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
        self.isolated_nodes = []
        
        
    def _Longest_branch_in_loop(self, junction_nodes):
        
        #Protected functions which return the largest branch in the loop

        tmp_ln = 0;
        source = 0;
        sync = 0;
        for i,j in itertools.combinations(junction_nodes, 2):
            shortest_path_ln = nx.shortest_path_length(self.G,  i, j)
                #Specific case to breaking the Longest branch
            if tmp_ln < shortest_path_ln:
                tmp_ln = shortest_path_ln
                source = i; sync = j;
        return [source,sync]



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

        if len(self.node_cordinates) != np.shape(self.node_connection)[0]:
            print("WARNING:: Different shapes of node_cordinates and node_connection ")
            print("Number of nodes are: " + str(len(self.node_cordinates)) + 
                " and shape of node_connection matrix is: ", str(np.shape(self.node_connection)))
            self.pos = {i: (self.node_cordinates[i,0], self.node_cordinates[i,1], self.node_cordinates[i,2]) 
                                    for i in range(np.shape(self.node_connection)[0])}
        else:
            self.pos = {i: (self.node_cordinates[i,0], self.node_cordinates[i,1], self.node_cordinates[i,2]) 
                        for i in range(len(self.node_cordinates))}

        for n, p in self.pos.items():
            self.G.nodes[n]['pos'] = p

    
    def Simplified_skeleton_for_Mariusz(self, tmp_G = 0, tmp_img = 0, update_G = False):

        """
        This function removes the intermediate nodes(with only two connectivity/degree)
        between two junction nodes and replace it with one branch according to the 
        requirement for Mariusz

        """
        
        print("WARNING: This function convert cordinates to 3D which will throw error for 2D plotting using nx_draw")

        
        if tmp_G == 0:
            G_simplified = self.G.copy();
            img = self.skeletonized_image;
        elif type(tmp_G) != 'int':
            G_simplified = tmp_G.copy();
            img = tmp_img;
            
        for i in self.G.degree:

            if i[1] == 2:

                node = i[0];
                edges = self.G.edges(node) if type(tmp_G) != 'int' else tmp_G.edges(node)
                edges = list(edges.__iter__())

                a0, b0 = edges[0];
                a1, b1 = edges[1];

                e0 = a0 if a0!=node else b0
                e1 = a1 if a1!=node else b1

                G_simplified.remove_node(node)


        #Using skan to extract the junction to junction/end end points

        data = summarize(Skeleton(img))

        list_1 = data['node-id-src']
        list_2 = data['node-id-dst']

        for i,j in zip(list_1,list_2):
            G_simplified.add_edge(i,j)

           
        if update_G is True:
                self.G = G_simplified
                
        return G_simplified
    
    
    
    def Remove_isolated_nodes(self):
        
        #This function delete the isolated nodes and return it

        G_copy = self.G.copy()
        self.isolated_nodes = list(nx.isolates(G_copy))
        G_copy.remove_nodes_from(self.isolated_nodes)
        self.G = G_copy            
    
    
    def Nx_to_Mayavi(self, tmp_G = 0):
        
        from mayavi import mlab
        nop={} 
       
        node_pos = self.pos
        for n in self.G.nodes(): 
            try: 
                nop[n]=node_pos[n] 
            except KeyError: 
                raise networkx.NetworkXError

    # numpy array of x,y,z positions in sorted node order

        xyz = np.array(list([[node_pos[v][0], node_pos[v][1], node_pos[v][2]] for v in self.G.nodes()]))
        
        scalars = np.ones((self.G.number_of_nodes()));
        
#         if self.outlet_node != 0:
#             for i in range(len(scalars)):
#                 if (node_pos[i] == self.inlet_node or node_pos[i] == self.outlet_node) :
#                     scalars[i] = 5;
#                 elif node_pos[i] in self.isolated_nodes:
#                     scalars[i] = 0

#         #Filling zeros at deleted nodes
#         try:
#             for i in self.isolated_nodes:

#                 xyz[i][0] = 0
#                 xyz[i][1] = 0
#                 xyz[i][2] = 0


#         except:
#             pass
        
        pts = mlab.points3d(
        xyz[:, 0],
        xyz[:, 1],
        xyz[:, 2],

        scale_factor=1,

        )

        pts.mlab_source.dataset.lines = np.array(list(self.G.edges()))
        tube = mlab.pipeline.tube(pts, tube_radius=0.1)
        mlab.pipeline.surface(tube, color=(1, 0., 0.))
        mlab.show()
    
    
    def Nx_to_vtk(self, addtional_suffix = ''): 
        """
        This function is written using VTK module
        Input:
            1. Network graph G
            2. isolated_nodes from "Remove_isolated_nodes" method to avoid the data error by fillinf zeros there
            3. Suffix (string) to save name of files

        """

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
            for i in self.isolated_nodes:
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
        writer.SetFileName(str(self.suffix).strip('.tif')+ addtional_suffix +".vtp");

        writer.SetInputData(edgeData)
        writer.Write()

    def Delete_Loops_from_Skeleton(self):
        """
            This function delete loops from the main skeleton in following ways:
            1. If loop is a part of the main Wormhole then the largest branch between 
            source and sync node will be cut at sync node

            2. If loop is not a part of wormhole and contains only 2 junction nodes 
            then the longest branch will be cut at Sync node
            
            3. If loop is not a part of wormhole and contains 3 or more junction nodes 
            then the longest branch will be cut at Sync node 


        """
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

            shortest_path = []

            if(if_main_wormhole_path == False):
                if (len(junction_nodes) > 2):
                    for i,j in itertools.combinations(junction_nodes, 2):
                        shortest_path = nx.shortest_path(self.G,  i, j)
                        if len([w for w in junction_nodes if w in shortest_path]) > 2:
                            loop_source_sync.extend([i,j])
                            break;
                        else:
                            #Specific case to breaking the Longest branch
                            loop_source_sync.extend(self._Longest_branch_in_loop(junction_nodes)) 
                            break;
                else:
                    loop_source_sync = junction_nodes

            elif (if_main_wormhole_path == True):
                    print("This cycle is in main wormhole path")
                    if (len(junction_nodes) > 2):
                        for i,j in itertools.combinations(junction_nodes, 2):
                            if i in self.dominant_wormhole and j in self.dominant_wormhole:
                                shortest_path = nx.shortest_path(self.G,  i, j)
                                if len([w for w in junction_nodes if w in shortest_path]) > 2:
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

                    
        
        
    def Find_Dominant_Wormhole(self, outlet_node = 0, given_inlet_node = 'NULL'):
        #This function works best for 3D images, because of uncertainity in choosing
        # Inlet and Outlet node in 2D
        #No. of slices should be in Z-Axis (i.e First axis here)
        
        print("Finding dominant wormhole")

        inlet_in_isolated_node = False
        
        
        non = len(self.node_cordinates) #non of total nodes
        
        print("%=====================Looking for Inlet node=========================%\n")
        inletz, _, _ = np.argmin(self.node_cordinates, axis=0)
        
        if given_inlet_node == 'NULL':
            if inletz in self.isolated_nodes: 
                print("Node-" + str(inletz) + "is not connected to main wormhole")
                inlet_in_isolated_node = True
                tmp_io = 0;

            while (inlet_in_isolated_node):

                tmp_io += 1 ;
                inletz = tmp_io + np.argmin(self.node_cordinates[tmp_io:non, 0]) #tmp_ioo is added because location of minimum z is always zero if element is excluded

                if inletz in self.isolated_nodes: 
                    print("Node-" + str(inletz) + "is not connected to main wormhole")
                    inlet_in_isolated_node = True
                else:
                    print("Node-"+ str(inletz) + " is connected and chosen as source node")
                    print("Cordinates of the inlet node is " + str(self.node_cordinates[inletz]))
                    inlet_in_isolated_node = False
                    break;
        else:
            inletz = given_inlet_node
            print("\n Node-"+ str(inletz) + ", provided by user is chosen as inlet node")
            print("Cordinates of the inlet node is " + str(self.node_cordinates[inletz]))
    
        print("%=====================Looking for Outlet node=========================%\n")     
        
        if outlet_node == 0:
            outlet_node_found = False
        else:
            outlet_node_found = True
            outletz = outlet_node
            print("\n Node-"+ str(outletz) + ", provided by user is connected and chosen outlet node")
            print("Cordinates of the outlet node is" + str(self.node_cordinates[outletz]))
            
                
        i = 0;
        while (outlet_node_found == False and i < 30):
                
            try:
                outletz = np.argmax(self.node_cordinates[0:non, 0])
                short_path = nx.shortest_path(self.G, inletz,  outletz)
                outlet_node_found = True
              
                print("\n Node-"+ str(outletz) + " is connected and chosen as outlet node")
                print("Cordinates of the outlet node is" + str(self.node_cordinates[outletz]))
                
            except Exception:
                print("Node-" + str(non) + " is not connected to Node-" + str(inletz))
                non -= 1;  i += 1;
                
        self.inlet_node = inletz;
        self.outlet_node = outletz;       
                
        self.dominant_wormhole = nx.shortest_path(self.G, inletz,  outletz)
        self.length_of_wormhole = nx.shortest_path_length(self.G, inletz,  outletz, weight='weight')
        return self.length_of_wormhole, self.dominant_wormhole
        
    def Tortuosity(self):
        """This function calculates the Tortuosity of the dominant wormhole
        
                        Actual length (through branching)
        Tortuosity =    -----------------------------------
                        Minimum length between inlet and outlet


        """
        return self.length_of_wormhole/self.node_cordinates[self.outlet_node][0]        
        
    def Wastefulness(self):
        #Using sum of sparce matrix but every connection is doubled so divided by half
        
        sum_of_all_branches_including_wormhole = self.node_connection.sum()/2
        
        return sum_of_all_branches_including_wormhole/self.node_cordinates[self.outlet_node][0]
        
    def Horton_number(self):
        print("Horton needs to be implimented")
        
    def Fractal_dimension(self):
        print("Thumbprint of God is yet to be implemented")
    