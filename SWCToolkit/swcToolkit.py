import numpy as np
import math
import pandas as pd
import os
import pickle
import copy
import sys
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

class Swctoolkit():
    """
    Originally cloned from: https://github.com/ClaytonBingham/SWCToolkit
    
        
    
    ```
    from roots.swcToolkit import swcToolkit
    swctool = swcToolkit()
    morph1 = swctool.load(morph1_filename)
    morph2 = swctool.load(morph2_filename)
    morph3 = swctool.combine_morphologies(morph1,morph2)
    swctool.to_swc(morph3,'morph3.swc')
    ```
    
    """
    def __init__(self):
        pass
    
    def shift_point(self,a,vec):
        """
        
        Shift a point around an origin according to user specified vector
        
        a - point to be shifted
        vector - x,y,z values to shift point by
        
        returns the shifted point
        
        """
        
        return([a[0]-vec[0],a[1]-vec[1],a[2]-vec[2]])
    
    def move_morphology(self,arbor,vector): #unit - microns
        """
        
        Shift a morphology around an origin according to user specified vector
        
        arbor - dictionary with branch number as keys and sections as values (each section is a list of points)
        vector - x,y,z values to shift all points by
        
        returns the shifted morphology
        
        """
        if any(isinstance(i, list) for i in arbor[0][0]):
            for branch in arbor.keys():
                for section in arbor[branch]:
                    arbor[branch][arbor[branch].index(section)] = [self.shift_point(point[:3],vector)+[point[-1]] for point in section]
        
        else:
            for branch in arbor.keys():
                arbor[branch]= [self.shift_point(point[:3],vector)+[point[-1]] for point in arbor[branch]]
        
        return(arbor)
    
    def rotate_points(self,origin, points, ax=0.0,ay=0.0,az=0.0):
        """
        
        Rotate a point counterclockwise by a given angle around a given origin.
        The angle should be given in radians.
        
        """
        r = Rotation.from_rotvec(np.array([ax,ay,az]))
        rotated_points = np.array(points) - np.array(origin)
        rotated_points = r.apply(rotated_points)
        return(rotated_points+np.array(origin))
    
    def rotate_morphology(self,arbor,origin,about_x=0.0,about_y=0.0,about_z=0.0): 
        """
        
        Rotate a morphology around an origin according to user specified elevation and azimuth
        
        arbor - dictionary with branch number as keys and sections as values (each section is a list of points)
        origin - point about which to rotate the morphology
        elevation - x,y plane in degrees
        azimuth - x,z plane in degrees
        
        returns the rotated morphology
        
        """
        newarbor = {}
        pts = []
        if any(isinstance(i, list) for i in arbor[0][0]):
            for branch in arbor.keys():
                for section in arbor[branch]:
                    for pt in section:
                        pts.append(pt[:3])
            
            rotated_points = self.rotate_points(origin,np.array(pts),ax=np.pi*about_x/180.0,ay=np.pi*about_y/180.0,az=np.pi*about_z/180.0)
            pt_index = 0
            for branch in arbor.keys():
                newarbor[branch] = []
                for section in arbor[branch]:
                    newarbor[branch].append([])
                    for pt in section:
                        newarbor[branch][-1].append(list(rotated_points[pt_index])+[pt[3]])
                        pt_index+=1
        
        else:
            for branch in arbor.keys():
                for pt in arbor[branch]:
                    pts.append(pt[:3])
            
            rotated_points = self.rotate_points(origin,np.array(pts),ax=np.pi*about_x/180.0,ay=np.pi*about_y/180.0,az=np.pi*about_z/180.0)
            pt_index = 0
            for branch in arbor.keys():
                newarbor[branch] = []
                for pt in arbor[branch]:
                    newarbor[branch].append(list(rotated_points[pt_index])+[pt[3]])
                    pt_index+=1
        
        
        return(newarbor)
    
    
    def split_branch_at_bifurcation(self,branch,b_point):
            return(branch[:branch.index(b_point)],branch[branch.index(b_point):])
    
    def reorder_arbor(self,arbor,a,b,c):
        newarbor = {}
        newbranch_count = 0
        for branch in arbor.keys():
            if arbor[branch] == a:
                newarbor[newbranch_count] = b
                newbranch_count+=1
                newarbor[newbranch_count] = c
                newbranch_count+=1
            else:
                newarbor[newbranch_count] = arbor[branch]
                newbranch_count+=1
        
        return(newarbor)
    
    def find_branch_w_bifurcation_point(self,arbor,b_point):
        for branch in arbor.keys():
            if b_point in arbor[branch][1:-1]:
                print(branch)
                return(branch)
        
        print(None)
        return(None)
    
    def split_branches_at_bifurcations(self,arbor):
        b_points = [arbor[b][0] for b in arbor.keys()]
#        new_arbor = copy.deepcopy(arbor)
        for point in b_points:
            a = self.find_branch_w_bifurcation_point(arbor,point)
            if a is not None:
                a = arbor[a]
            
            else:
                continue
            
            b,c = self.split_branch_at_bifurcation(a,point)
            arbor = self.reorder_arbor(arbor,a,b,c)
        return(arbor)
    
    def prepend_headers(self,filename):
        with open(filename,'r',encoding='utf-8') as f:
            dat = f.readlines()
            if dat[0] == 'n T x y z R P\n':
                return()
        with open(filename,'w',encoding='utf-8') as f:
            f.write('n T x y z R P\n')
            for line in dat:
                f.write(line)
                
    
    def cut_headers(self,filename):
        with open(filename,'r',encoding='utf-8') as f:
            dat = f.readlines()
        if dat[0] == 'n T x y z R P\n':            
            with open(filename,'w',encoding='utf-8') as f:
                for line in dat[1:]:
                    f.write(line)
            return()
        else:
            return()
    
    def sort_index(self,li):
        newli = sorted(li)
        return([li.index(l) for l in newli])
    
    def find_children(self,bif,np_dict):
        children = []
        for n in np_dict.keys():
            if np_dict[n] == bif:
                children.append(n)
        return(children)
    

    
    
    def load_morph_and_section_types(self,fname):
        label_dict=dict(zip([1,2,3,4,5,6],['soma','axon','basal','apical','tuft','custom']))
        morph = self.load(fname)
        dat = self.load_swc_as_table(fname)
        labels = {}
        for key in morph.keys():
            labels[key] = []
            for section in morph[key][1:]:
                labels[key].append(label_dict[int(dat.loc[(dat['x']==section[0])&(dat['y']==section[1])&(dat['z']==section[2])]['T'])])
        
        return(labels,morph)
        
    def load_swc_as_table(self,fname):
        self.prepend_headers(fname)
        dat = pd.read_table(fname,delim_whitespace=True)
        self.cut_headers(fname)
        for col in ['n','T','P']:
            dat[col] = dat[col].astype('int')
        for col in ['x','y','z','R']:
            dat[col] = dat[col].astype('float')
        
        return(dat)
    
    def load(self,fname,asTree=True,interbifurcated=True):
        dat = self.load_swc_as_table(fname)
        if asTree:
            if not interbifurcated:
                arbor = {}
                n_point = dict(zip(dat['n'],zip(dat['x'],dat['y'],dat['z'],dat['R'])))
                np_dict = dict(zip(dat['n'],dat['P']))
                branchnumber = 0
                arbor[0] = [n_point[1]]
                for n in np_dict.keys():
                    if np_dict[n] == -1:
                        continue
                    if np_dict[n] != n-1:
                        branchnumber+=1
                        arbor[branchnumber] = [n_point[np_dict[n]]]
                    
                    arbor[branchnumber].append(n_point[n])
                
                
                return(arbor)
            
            else:
                arbor = {}
                n_point = dict(zip(dat['n'],zip(dat['x'],dat['y'],dat['z'],dat['R'])))
                np_dict = dict(zip(dat['n'],dat['P']))
                branchnumber = 0
                arbor[0] = [n_point[1]]
                for n in np_dict.keys():
                    if np_dict[n] == -1:
                        continue
                    if np_dict[n] != n-1:
                        branchnumber+=1
                        arbor[branchnumber] = [np_dict[n]]
                    
                    arbor[branchnumber].append(n)
                
                terminals = [arbor[b][-1] for b in arbor.keys()]
                
                arbor = {}
                pn_dict = dict(zip(dat['P'],dat['n']))
                p_point = dict(zip(dat['P'],zip(dat['x'],dat['y'],dat['z'],dat['R'])))
                branchnumber = 0
                bifrows = [1]
                for n in np_dict.keys():
                    if n == 1:
                        continue
                    if np_dict[n] != n-1:
                        bifrows.append(np_dict[n])
                
                bifrows.sort()
                children = {}
                for b in bifrows:
                    children[b] = self.find_children(b,np_dict)
                
                for childs in children.values():
                    for child in childs:
                        arbor[branchnumber] = [np_dict[child],child]
                        branchnumber+=1
                
                for n,row in enumerate(dat['n']):
                    if row-1 not in terminals and row-1 not in bifrows:
                        for branch in arbor.keys():
                            if np_dict[row] == arbor[branch][-1]:
                                arbor[branch].append(row)
                                break
                
                for branch in arbor.keys():
                    arbor[branch] = [n_point[n] for n in arbor[branch]]
                
                
                return(arbor)
        
        else:
            x = {0:[dat['x'][0]]}
            y = {0:[dat['y'][0]]}
            z = {0:[dat['z'][0]]} 
            r = {0:[dat['R'][0]]}
            
            
            for p,point in enumerate(dat['x'][1:]):
                inc = dat['P'][p+1]
                if inc != dat['n'][p]:
                    bnum = len(x.keys())
                    x[bnum] = [dat.query('n == '+str(dat['P'][p+1]))['x'].values[0],dat['x'][p+1]]
                    y[bnum] = [dat.query('n == '+str(dat['P'][p+1]))['y'].values[0],dat['y'][p+1]]
                    z[bnum] = [dat.query('n == '+str(dat['P'][p+1]))['z'].values[0],dat['z'][p+1]]
                    r[bnum] = [dat.query('n == '+str(dat['P'][p+1]))['R'].values[0],dat['R'][p+1]]
                    continue
                else:
                    bnum = len(x.keys())-1
                    x[bnum].append(dat['x'][p+1])
                    y[bnum].append(dat['y'][p+1])
                    z[bnum].append(dat['z'][p+1])
                    r[bnum].append(dat['R'][p+1])
            return(x,y,z,r)
    
    
    def to_swc(self,arbor=[],target='swcTest.swc'):
        
        """
        n T x y z R P
        
        n is an integer label that identifies the current point and increments by one from one line to the next.
        
        T is an integer representing the type of neuronal segment, such as soma, axon, apical dendrite, etc. The standard accepted integer values are given below.
        
            0 = undefined
            1 = soma
            2 = axon
            3 = basal dendrite
            4 = apical dendrite
            5 = fork point
            6 = end point
            7 = custom
        
        x, y, z gives the cartesian coordinates of each node.
        
        R is the radius at that node.
        P indicates the parent (the integer label) of the current point or -1 to indicate an origin (soma). 
        """
        
        try:
            with open(target,'wb') as f:
                points = [arbor[0][0]]
                n = [1]
                T = 2
                R = 0.7
                P = [-1]
                nind = 2
                Pind = 1
                f.write(str(n[-1])+' '+str(T)+' '+str(points[0][0])+' '+str(points[0][1])+' '+str(points[0][2])+' '+str(R)+' '+str(P[-1]))
                f.write('\n')
                for key in arbor.keys():
                    for p,point in enumerate(arbor[key][1:]):
                        points.append(point)
                        n.append(nind)
                        if p == 0 and key != 0:  # JMB 2025.04.02: replaced 'is not' with '!='
                            P.append(points.index(arbor[key][0])+1)
                        else:
                            P.append(Pind)
                        nind+=1
                        Pind+=1
                        f.write(str(n[-1])+' '+str(T)+' '+str(point[0])+' '+str(point[1])+' '+str(point[2])+' '+str(R)+' '+str(P[-1]))
                        f.write('\n')
        
        except:
            with open(target,'w',encoding='utf-8') as f:
                points = [arbor[0][0]]
                n = [1]
                T = 2
                R = 0.1
                P = [-1]
                nind = 2
                Pind = 1
                f.write(str(n[-1])+' '+str(T)+' '+str(points[0][0])+' '+str(points[0][1])+' '+str(points[0][2])+' '+str(R)+' '+str(P[-1]))
                f.write('\n')
                for key in arbor.keys():
                    for p,point in enumerate(arbor[key][1:]):
                        points.append(point)
                        n.append(nind)
                        if p == 0 and key != 0:  # JMB 2025.04.02: replaced 'is not' with '!='
                            P.append(points.index(arbor[key][0])+1)
                        else:
                            P.append(Pind)
                        nind+=1
                        Pind+=1
                        f.write(str(n[-1])+' '+str(T)+' '+str(point[0])+' '+str(point[1])+' '+str(point[2])+' '+str(R)+' '+str(P[-1]))
                        f.write('\n')
    
    def find_best_merge_point(self,a,b):
        points_a = []
        for branch in a.keys():
            points_a+=a[branch]
        
        points_a = np.array([np.array(item[:3]) for item in points_a]).reshape(len(points_a,3))
        
        points_b = []
        for branch in b.keys():
            points_b+=b[branch]
        
        points_b = np.array([np.array(item[:3]) for item in points_b]).reshape(len(points_b,3))
        
        distances = cdist(points_a,points_b)
        return(points_a[distances.argmin()])
    
    
    
    def combine_morphologies(self,a,b,merge_point=None):
        if merge_point == None:
            b[0] = self.find_best_merge_point(a,b)+b[0]
        else:
            b[0] = merge_point+b[0]
        
        for branchnumber in b.keys():
            a[len(a.keys())+branchnumber] = b[branchnumber]
        
        return(a)
    

    
    
    
    def convert_to_obj_lines(self,input_filename,output_filename):
        """
        Added 2025.04.02: converts a given SWC file into OBJ lines
        Usage: 
            from swcToolkit import swcToolkit
            swcTool = swcToolkit()
            swcTool.convert_to_obj_lines("./test_retinal_cells.swc" , "./test_converted_output.obj")            
        Status: Tested! 2025.04.02

        Parameters
        ----------
        input_filename : TYPE
            DESCRIPTION. input file containing the SWC file
        output_filename : TYPE
            DESCRIPTION. output file that will contain the OBJ file

        Returns
        -------
        None.

        """
        data_raw = np.loadtxt(input_filename)

        points = data_raw[:, 2:5]  # columns 2, 3, 4 are X, Y, Z
        connection = data_raw[:, 6].astype(int)  # parent IDs, converted to int
    
        #create a mapping from swc ID to vertex index
        id_to_index = {int(data_raw[i, 0]): i + 1 for i in range(len(data_raw))}  # +1 for 1-based index in OBJ
    
        with open(output_filename, 'w') as f:
            for point in points:
                txt_str = 'v ' + str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2]) + '\n'
                f.write(txt_str)
                
            for index, parent_id in enumerate(connection):
                if parent_id == -1 or parent_id == 0:
                    continue
                if parent_id in id_to_index:
                    parent_index = id_to_index[parent_id]
                    txt_str = f'l {index + 1} {parent_index}\n'
                    f.write(txt_str)
        asdf = 1
        
    
    def read_swc_file(self, file_path):
        """
        Added 2025.04.02: Reads the content of a given SWC file 
        Status: TESTED and FUNCTIONAL 
        Parameters
        ----------
        file_path : TYPE
        DESCRIPTION.

        Returns
        -------
        swc_data : TYPE Pandas Dataframe
        DESCRIPTION. swc data in Pandas Dataframe format
        """
        
        # Load SWC file using pandas
        swc_data = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None)
        
        # Replace 'NA' and 'nan' values with -1
        swc_data = swc_data.replace(['NA', np.nan], -1)
        
        # Convert the pandas DataFrame to a NumPy array
        swc_data = swc_data.to_numpy()
        
        # Create a structured array with the specified header
        swc_data = np.core.records.fromarrays(swc_data.T,
                                              names='id, type, x, y, z, radius, parent',
                                              formats='i4, i4, f4, f4, f4, f4, i4')
        return swc_data


    def create_mesh_from_swc(self, swc_data, minRadius=0.005):
        """
        Creates a nesh from the Pandas dataFrame of a SWC file
        Status: TESTED and FUNCTIONAL
        
        Parameters
        ----------
        swc_data : TYPE Pandas dataFrame
            DESCRIPTION. SWC data in a Pandas dataFrame format
        minRadius : TYPE, optional
            DESCRIPTION. The default is 0.005. Lower values will result in a larger file containing more triangles of smaller sizes.

        Returns
        -------
        combined_mesh : TYPE
            DESCRIPTION. A concatenated Trimesh.mesh object 

        """
        
        import trimesh

        # Create an empty list to store all mesh objects
        meshes = []
    
        # Process nodes
        for i, node in enumerate(swc_data):
            # Create a sphere for each node
            sphere = trimesh.creation.icosphere(subdivisions=2, radius=max(node['radius'], minRadius))
            sphere.apply_translation([node['x'], node['y'], node['z']])
            meshes.append(sphere)
    
            # Create a cylinder for each edge
            if node['parent'] != -1:
                parent_node = swc_data[np.where(swc_data['id'] == node['parent'])[0][0]]
    
                # Calculate cylinder properties
                start = np.array([node['x'], node['y'], node['z']])
                end = np.array([parent_node['x'], parent_node['y'], parent_node['z']])
                length = np.linalg.norm(end - start)
    
                if length > 0:  # Check if the length is greater than zero before proceeding
                    direction = (end - start) / length
                    radius = (max(node['radius'], minRadius) + max(parent_node['radius'], minRadius)) / 2
    
                    # Create the cylinder
                    cylinder = trimesh.creation.cylinder(radius=max(radius, minRadius), height=length, sections=16)
    
                    try:
                        cylinder.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], direction))
                    except np.linalg.LinAlgError:
                        # Alternative method to align the vectors
                        axis = np.cross([0, 0, 1], direction)
                        angle = np.arccos(np.dot([0, 0, 1], direction))
                        cylinder.apply_transform(trimesh.transformations.rotation_matrix(angle, axis))
    
                    cylinder.apply_translation((start + end) / 2)
    
                    # Add the cylinder to the list of meshes
                    meshes.append(cylinder)
    
        # Combine all the meshes into a single mesh object
        combined_mesh = trimesh.util.concatenate(meshes)
    
        return combined_mesh
    
    def convert_swc_to_obj(self, swc_file, obj_file):
        """
        # # Example usage
        swc_file = 'SWC_input/output9.swc'
        obj_file = 'OBJ_output/output9.obj'
        
        convert_swc_to_obj(swc_file, obj_file)

        Parameters
        ----------
        swc_file : TYPE
            DESCRIPTION. Name and path of input SWC file
        obj_file : TYPE
            DESCRIPTION. Name and path of output OBJ file

        Returns
        -------
        None.

        """

        swc_data = self.read_swc_file(swc_file)
        mesh = self.create_mesh_from_swc(swc_data)
        mesh.export(obj_file)
    

    def split_swc(self, input_file, output_path): 
        """
        Function separates neuron morphologies contained in one SWC file.
        :param input_file: Path to the input SWC file.
        :param output_path: Path (directory) to save the output file(s).
        """
        #JMB import sys
        #JMB from swcToolkit import swcToolkit  # Import the SWC toolkit
        # JMB import os  # Import os for file operations

        # Extract the file name without extension
        file_name = os.path.splitext(input_file)[0]
        input_path, file_name = os.path.split(file_name)
        orphan_file = 'orphan'
        #swc = swcToolkit() # load swc tool kit for use throughout
        swc_data = self.read_swc_file(input_file)  # Read SWC data
        recursion_limit = len(swc_data) + 1 # changes the recursion limit so that recursive function can be called more than 5 times.
        
        # If the number of points in the SWC file is greater than 20
        # increase the recursion limit to prevent stack overflow.
        if recursion_limit > 20 :  
            sys.setrecursionlimit(recursion_limit)

        # Find soma indices in the SWC data
        # structure_identifier = 1 # column # in swc format for structure identifier
        # parent_sample = 6 # column # in swc format for parent node
        # sample_num = 0 # column # in swc format for sample number
        soma_indices = self.soma_id_and_orphan(swc_data) # array of soma indices from swc data
        dictionary_data = self.data_into_dictionary(swc_data) # 
        num_soma = 0
        returned_nodes = 0
        sample_col = 0
        num_orphans = 0
        r_orphan_nodes = 0

        print (' %d soma node(s) and %d orphan node(s) identified in this file.' % (len(soma_indices[0]), len(soma_indices[1])))

        # Separate neurons based on soma indices and export
        for idx in range(len(soma_indices[0])): 
            num_nodes = 0
            iteration = 1
            swc_data_index = soma_indices[0][idx]
            # print(swc_data_index, 'is the current soma index in the swc_data.\n')
            tree = {1: [swc_data[swc_data_index][sample_col]]}
            new_tree = self.soma_tree(tree, iteration, dictionary_data, num_nodes)
            unfiltered_nrn_data = self.data_idx_to_data(swc_data, new_tree[0])
            returned_nodes += new_tree[1]
            nrn_data = self.redef_idx(unfiltered_nrn_data)
            num_soma += 1 
            # ORIG:         file(nrn_data, output_path, num_soma, file_name, soma_file = True)  # Write se parated neurons to files
            #self.write_to_swc(nrn_data, output_path, num_soma, file_name, soma_file = True)  # Write se parated neurons to files
            self.write_to_swc(nrn_data = nrn_data, output_path = output_path, file_num =  num_soma , input_name = file_name, soma_file = True) 
            #self.write_to_swc(nrn_data, output_path, num_soma, file_name)  # Write se parated neurons to files
            print(returned_nodes, ' nodes were identified to be connected to a different soma.\n')

        for idx in range(len(soma_indices[1])): 
            num_nodes = 0
            iteration = 1
            swc_data_index = soma_indices[1][idx]
            # print(swc_data_index, 'is the current orphan start node index in the swc_data.\n')
            tree = {1: [swc_data[swc_data_index][sample_col]]}
            new_tree = self.soma_tree(tree, iteration, dictionary_data, num_nodes)
            unfiltered_nrn_data = self.data_idx_to_data(swc_data, new_tree[0])
            r_orphan_nodes += new_tree[1]
            nrn_data = self.redef_idx(unfiltered_nrn_data)
            num_orphans += 1 
            #ORIG:         file(nrn_data, output_path, num_orphans, file_name, soma_file = False)  # Write separated neurons to files
            #self.write_to_swc(nrn_data, output_path, num_orphans, file_name, soma_file = False)  # Write separated neurons to files
            self.write_to_swc(nrn_data = nrn_data, output_path = output_path, file_num =  num_orphans , input_name = file_name, soma_file = False)  # Write separated neurons to files
            #self.write_to_swc(nrn_data, output_path, num_orphans, orphan_file)  # Write separated neurons to files
            print(r_orphan_nodes, ' nodes connected to an orphan node.\n')

        # If all nodes have been processed successfully, indicate success
        # Check if all nodes were separated
        if (returned_nodes + r_orphan_nodes == len(swc_data)): 
            print('All nodes separated successfully.\n')
            print('Total number of nodes in swc data:', len(swc_data))
            print('Total number of nodes separated:', returned_nodes + r_orphan_nodes)
            print('Total number of soma nodes: ', len(soma_indices[0]))
            print('Total number of orphan nodes: ', len(soma_indices[1]))
            print('Total number of nodes connected to soma:', returned_nodes)
            print('Total number of nodes connected to orphan nodes:', r_orphan_nodes)
        else: 
            print('Not all nodes were separated. \nPlease check the SWC data for errors, or orphan node connections.')
            print('Total number of nodes in swc data:', len(swc_data))
            print('Total number of nodes separated:', returned_nodes + r_orphan_nodes)
            print('Total number of soma nodes: ', len(soma_indices[0]))
            print('Total number of orphan nodes: ', len(soma_indices[1]))
            print('Total number of nodes connected to soma:', returned_nodes)
            print('Total number of nodes connected to orphan nodes:', r_orphan_nodes)



    def soma_id_and_orphan(self, data): # Seems to work fine
        """
        Find soma indices in the SWC data.
        :param data: SWC data as a pandas DataFrame.
        :return: List of soma indices.
        """
        structure_identifier = 1 # column # in swc format for structure identifier
        parent_sample = 6 # column # in swc format for parent node
        # sample_num = 0 # column # in swc format for sample number
        soma = []
        orphans = []
        for index in range(len(data)):
            if data[index][structure_identifier] == 1 and data[index][structure_identifier] == 1:  # Check if the node is a soma
                soma.append(index)  # Append the soma node to the list of starting points
            elif data[index][parent_sample] == -1 and data[index][structure_identifier] != 1:
                orphans.append(index)
        return soma, orphans

    def find_indices(self, arr, targets):
        """Finds all indices in an array that have the values listed in the targets array
        :param arr: An array of data to be searched through
        :param targets: An array of values to search for
        :return output: An array with the indices of every value that matches a value described in the target array"""
        target_set = set(targets)
        output = [i + 1 for i, v in enumerate(arr) if v in target_set]
        if len(output) == 0: 
            return(None)
        return output
    # Example Usage
    # print(find_indices(p, [12]))

    def add_dict_entries(self, dictionary): 
        """Returns one array as a combination of many smaller arrays each with individual dictionary keys
        :param dictionary: The dictionary full of arrays to combine into one array
        :return output: One array consisting of all array values in the dictionary from every key"""
        output = []
        for i in range(len(dictionary) - 1): 
            output += dictionary[i + 1]
        return output

    def soma_tree(self, dictionary, iteration, data, num_nodes): 
        """Separates individual neurons from network of neurons using swc data and dictionary keys. It finds all nodes connected to the previously located nodes. 
        :param dictionary: A dictionary with each key as an array of nodes connected to the nodes described in the previous key
        :param iteration: The number of iterations the recursive function has completed .
        :param data: The data of parent nodes in dictionary format to be updated each iteration.
        :param num_nodes: The total number of nodes that have been filtered through and isolated into an individual neuron. 
        :return output_arr: A single array of all indices in the swc_data that compose a single neuron.
        :return num_nodes: The total number of nodes that were isolated after isolating one individual neuron."""
        if dictionary.get(iteration) == None: # Checks and returns all indices of raw swc data of a specific soma starting point
            output_arr = self.add_dict_entries(dictionary)
            num_nodes += len(output_arr)
            return output_arr, num_nodes
        else: 
            parents = dictionary.get(iteration) # Recursive function that finds the indices of nodes that make up individual neurons with soma as starting points
            new_parents = self.find_indices(data['p'], parents)
            dictionary.update({(iteration + 1): new_parents})
            iteration += 1
            return self.soma_tree(dictionary, iteration, data, num_nodes)

    def data_into_dictionary(self, data): 
        """Reorganizes raw data into a dictionary format for use in other functions.
        :param data: swc data as a pandas dataFrame.
        :return output_dictionary: dictionary with sample numbers and parent node information of swc data."""
        sample_num = 0 # define swc convention constants
        parent_sample = 6 
        indices_arr = [] # initialize arrays for use in function. 
        parents_arr = []
        for i in range(len(data)): 
            indices_arr.append(data[i][sample_num])
            parents_arr.append(data[i][parent_sample])
        output_dictionary = {'i': indices_arr, 'p': parents_arr}
        return(output_dictionary)
        
    def data_idx_to_data(self, swc_data, idx_arr):
        """Gathers all raw swc data from pandas dataframe from array that holds index values of wanted data. 
        :param swc_data: Raw swc_data in pandas datafram format
        :param idx_array: An array of sample number indices of individual neuron cells.
        :return: Pandas Dataframe of swc raw data with only the indices of wanted data."""
        filtered_data = []
        # print(idx_arr)
        for index in idx_arr: 
            filtered_data.append(swc_data[index - 1])
        return filtered_data

    def redef_idx(self, data): 
        """Redefine indices in the SWC data to match swc format convensions after node isolation.
        :param data: swc data as a pandas dataFrame.
        :return: DataFrame with redefined indices."""
        from pprint import pprint
        import numpy as np
        new_sample_num = 0 # Initialize first sample number for swc format
        data = self.sort_by_id(data)
        # pprint(data)
        for i in range(len(data)): 
            new_sample_num += 1
            current_id = data[i]['id']
            # print("Checking Parents equal to:", current_id) 
            for j in range(len(data)): 
                current_parent = data[j]['parent']
                # print(current_parent)
                if (current_parent == current_id):
                    data[j]['parent'] = new_sample_num
            data[i]['id'] = new_sample_num
        return data      

    def sort_by_id(self, data): 
        """Sorts SWC data by ID value so it starts at the lowest and ends at the highest value to counter overwriting.
        :param data: swc data as a pandas dataFrame.
        :return: DataFrame sorted by ID."""
        id_list = [] 
        new_data = [] 
        for i in range(len(data)): 
            id_list.append(data[i]['id'])
        id_list.sort()
        # print(len(id_list))
        # print(id_list, '\n')
        for i in range(len(id_list)): 
            find_value = id_list[i]
            for j in range(len(data)):
                if data[j]['id'] == find_value:
                    new_data.append(data[j])
                    break
        return new_data

    def write_to_swc_old(self, nrn_data, output_path, file_num, input_name): 
        """
        Write separated neurons to files.
        :param nrn: Separated neuron. 
        :param output_path: Path to save the output files.
        """
        import os
        if not os.path.exists(output_path): # Create export path
            os.makedirs(output_path)
        
        file_name = os.path.join(output_path, f'{input_name}_{file_num}.swc') # creates file name
        # f'
        with open(file_name, 'w') as f: # exports new file to export path
            f.write('# ID, type, x, y, z, radius, parent\n')
            for node in nrn_data:
                f.write(' '.join(map(str, node)) + '\n')

    def write_to_swc(self, nrn_data, output_path, file_num, input_name, soma_file): 
        """
        Write separated neurons to files.
        :param nrn: Separated neuron. 
        :param output_path: Path to save the output files.
        """
        #import os
        print('output path:', output_path, ' input name:', input_name, ' file num:', file_num, ' soma file:', soma_file)
        if not os.path.exists(output_path): # Create export path
            os.makedirs(output_path)
        if soma_file == True: # If the file is a soma file, add soma to the file name
            file_name = os.path.join(output_path, f'{input_name}_{file_num}.swc') # creates file name
        elif soma_file == False: 
            file_name = os.path.join(output_path, f'{input_name}_orphan_{file_num}.swc') # creates file name for orphan nodes
        with open(file_name, 'w') as f: # exports new file to export path
            print('Exporting to:', file_name)
            f.write('# ID, type, x, y, z, radius, parent\n')
            for node in nrn_data:
                f.write(' '.join(map(str, node)) + '\n')