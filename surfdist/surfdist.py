def dist_calc(surf, cortex, source_nodes, dist_type = "min"):

    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "dist_type" specifies whether to calculate "min", "mean", "median", or "max" distance values 
    from a region-of-interest. If running only on single node, defaults to "min".
    """
    
    import gdist
    from utils import surf_keep_cortex, translate_src, recort
    import numpy as np

    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    translated_source_nodes = translate_src(source_nodes, cortex)

    if len(src_new) == 1:
        
        dist_type = "min"
        print "calculating min for single node ROI"
        
    if dist_type == "min":
    
        data = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)
    
    else:
        
        data_nodes = []
        
        for i in translated_source_nodes:
            
            data_nodes.append(gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = np.array(i, ndmin=1)))
            
        data_nodes = np.array(data_nodes)
        
        if dist_type == "mean":
        
            data = np.mean(data_nodes, axis = 0)
        
        if dist_type == "median":
            
            data = np.median(data_nodes, axis = 0)
            
        if dist_type == "max":
        
            data = np.max(data_nodes, axis = 0)
            
    dist = recort(data, surf, cortex)

    del data

    return dist


def zone_calc(surf, cortex, src):
    """
    Calculate closest nodes to each source node using exact geodesic distance along the cortical surface.
    """
    import gdist
    from utils import surf_keep_cortex, translate_src, recort
    import numpy as np

    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    dist_vals = np.zeros((len(source_nodes), len(cortex_vertices)))

    for x in range(len(source_nodes)):
        
        translated_source_nodes = translate_src(source_nodes[x], cortex)
        dist_vals[x, :] = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)

    data = np.argsort(dist_vals, axis=0)[0, :] + 1

    zone = recort(data, surf, cortex)

    del data

    return zone

def dist_calc_matrix(surf, cortex, labels, exceptions = ['Unknown', 'Medial_Wall'], dist_type = "min"):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "labels" specifies the freesurfer label file to use. All values will be used other than those
    specified in "exceptions" (default: 'Unknown' and 'Medial_Wall').
    "dist_type" specifies whether to calculate "min", "mean", "median", or "max" distance values 
    from a region-of-interest. If running only on single node, defaults to "min".
    Returns: dist_matrix and rois names
    """
    
    import gdist
    from utils import surf_keep_cortex, translate_src, recort
    import numpy as np

    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex) 

    label_list = sd.load.get_freesurfer_label(labels, verbose = False)
    rs = np.where([a not in exceptions for a in label_list])[0]
    rois = [label_list[r] for r in rs]

    translated_source_nodes = []
    
    for region in rois:
    
        source_nodes = sd.load.load_freesurfer_label(labels, region, cortex)
        translated_source_nodes.append(translate_src(source_nodes, cortex))
    
    data_matrix = np.zeros((len(rois), len(rois)))
    
    if dist_type == "min":
    
        for i in translated_source_nodes:
        
            data = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = i)
            
            for n, j in enumerate(translated_source_nodes):
                
                data_matrix[n,:] = np.min(data[j])
    
    else:
                
        for m, i in enumerate(translated_source_nodes):
            
            data_nodes = []
            
            for j in i:
            
                data_nodes.append(gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = np.array(j, ndmin=1)))
            
            data_nodes = np.array(data_nodes)            
        
            if dist_type == "mean":

                data = np.mean(data_nodes, axis = 0)

            if dist_type == "median":

                data = np.median(data_nodes, axis = 0)

            if dist_type == "max":

                data = np.max(data_nodes, axis = 0)
                
            del data_nodes
                
            for n, k in enumerate(translated_source_nodes):
                
                if dist_type == "mean":

                    data_matrix[m,n] = np.mean(data[k])

                if dist_type == "median":

                    data_matrix[m,n] = np.median(data[k])

                if dist_type == "max":

                    data_matrix[m,n] = np.max(data[k])

            del data
    
    # take mean of each value to render matrix symmetrical:
    data_matrix = (data_matrix + data_matrix.T) / 2.
            
    return dist_matrix, rois

