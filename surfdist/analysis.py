import gdist
import numpy as np
import nibabel as nib
from surfdist.utils import surf_keep_cortex, translate_src, recort,AnatomyInputParser,LabelInputParser, create_networkx_graph
import surfdist as sd
from surfdist.load import load_cifti_labels,load_freesurfer_label,get_freesurfer_label, load_gifti_labels,load_FS_annot
#### multiprocessing is used in cifti and gifti distance matrix calculations 
from multiprocessing.pool import Pool as ProcessPool
import multiprocessing


def dist_calc(surf, cortex, source_nodes,recortex=True,maxDist=False):

    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "dist_type" specifies whether to calculate "min", "mean", "median", or "max" distance values
    from a region-of-interest. If running only on single node, defaults to "min".
    """
    surf=AnatomyInputParser(surf)
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    translated_source_nodes = translate_src(source_nodes, cortex)
    ### check to see if we output ditance at all points or within a given radius
    if maxDist==False:
        dist = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)
        dist[~np.isfinite(dist)]=0 ### remove any nan's or infinities if they exist
    else:
        dist = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes,max_distance=maxDist)
        dist[dist==np.max(dist)]=0

    
    if recortex==True:
        dist = recort(dist, surf, cortex)

    return dist

def zone_calc(surf, cortex, source_nodes):
    """
    Calculate closest nodes to each source node using exact geodesic distance along the cortical surface.
    """
    surf=AnatomyInputParser(surf)

    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    dist_vals = np.zeros((len(source_nodes), len(cortex_vertices)))
    for x in range(len(source_nodes)):
        translated_source_nodes = translate_src(source_nodes[x], cortex)
        dist_vals[x, :] = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)

    data = np.argsort(dist_vals, axis=0)[0, :] + 1

    zone = recort(data, surf, cortex)

    del data

    return zone

def dist_calc_matrix(AnatSurf,LabelInput,hemi,exceptions=[],n_cpus=1,fsCort=None,hires=False,maxDist=False):
    """
    Calculate a distance matrix between a set of ROIs defined by a set of labels. 
    
    Inputs are as follows:
    
    AnatSurf: An anatomical surface file which such as freesurfer output or a gifti file. 
    Alternatively a tuple of len(2) containing the vertex indices and vertex faces can be passed. 

    LabelInput: Can be a gifti label file, cifti label file, or freesurfer annotation file. 
    Alternatively list containing three lists:
    The first being the keyu values, the second the source nodes of the label, and the third a list of cortical vertices

    USING FREESURFER ANNOTATIONS: You must specicy the medial wall label of the annotation you're using as the first entry in the exceptions list
    Gifti and Cifti default to a medial wall key value of '???'
    
    hemi: the string 'L' or 'R' specifying which hemisphere the labels are being extracted from

    exceptions: A list of areas to exclude from the distance matrix calculation. 
    There is no need to specify the medial wall. The function excludes that by default

    n_cpus: The number of cpus to use when calculating the distance matrix. 
    
    returns:
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """
    print(f'using {n_cpus} cpus')
    surf=AnatomyInputParser(AnatSurf)
    print(exceptions)
    labels,medialWall=LabelInputParser(LabelInput,hemi,exceptions)
    

    ctx = np.array(range(len(surf[0])))
    cortex = np.delete(ctx, medialWall)

    #### let's edit to allow us to specify a cortex mask here. 

    if hires==True:
        nodes=list(cortex)
        labels=dict(zip(nodes,nodes))

        
    else:

        if len(exceptions)>0:
            print(exceptions)
            for i in exceptions.copy():
                if i=='???':
                    pass
                else:
                    del labels[i]
        
        for l in labels.copy():
            if labels[l].squeeze().shape[0]==0:
                print(f'{l} is an empty label')
                del labels[l]
    
        nodes= list(labels.values())

    if maxDist==False:
        params=[[surf,cortex,nodes[i],'recort=False'] for i in range(len(nodes))]
    else:
        print('using MAXDIST')
        # limit=f'maxDist={maxDist}'
        params=[[surf,cortex,nodes[i],'recort=False','maxDist=',int(maxDist)] for i in range(len(nodes))]
    
    with ProcessPool(processes=n_cpus) as pool:
        dist_roi=pool.starmap(dist_calc,params)
    dist_roi=np.array(dist_roi)
     ##Calculate min distance per region:
    dist_mat = []
    for roi in labels:
        source_nodes=labels[roi]
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    if hires==False:
        return dist_mat,list(labels.keys())
    else:
        return dist_mat,nodes

def GeoDistHeuristic(verts,faces,source,target):
    print(source,target)
    d=gdist.compute_gdist(verts,faces,np.array([source],dtype=np.int32),np.array([target],dtype=np.int32))
    return d[0]

import networkx as nx
def shortest_path(surf,cortex,source,target,method='dijkstra'):
    """ Calculate the shortest path on the cortex using NetworkX and dijkstra's algorithm
    This function returns the vertices in the shortest path between two vertices on the surface.
    It does not make use of the gdist package except for in loading the anatomical files. 
    A cortex mask is required to ensure the mask does not run over the medial wall
     Methods permit shortest path calculation with Dijkstra, bellman ford, or A*. 
      A* may take longer due to the heuristic caluclation of geodesic iteraton between iterations. """
    surf=AnatomyInputParser(surf)
    
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    translated_source = translate_src(source, cortex)[0]
    translated_target = translate_src(target, cortex)[0]
    
    Graph =create_networkx_graph(cortex_vertices, cortex_triangles)

    if method=='dijkstra':
        path=nx.shortest_path(Graph,translated_source,translated_target)
    elif method=='bmf':
        path=nx.shortest_path(Graph,translated_source,translated_target,method='bellman-ford')
    elif method=='Astar':
        path = nx.astar_path(Graph, translated_source, translated_target, 
                         heuristic=lambda node1, node2: GeoDistHeuristic(cortex_vertices, cortex_triangles, node1, node2))
    cortPath=[]
    for p in path:
        #### change to be the value of cort indexed by p
        precort=np.zeros(len(cortex_vertices))
        precort[p]=1
        precort=recort(precort,surf,cortex)
        cortPath.append(np.where(precort==1)[0])
    cortPath=np.hstack(cortPath)

    return cortPath
    