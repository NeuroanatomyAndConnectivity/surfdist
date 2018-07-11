import gdist
from surfdist.utils import surf_keep_cortex, translate_src, recort
import numpy as np
    
def dist_calc(surf, cortex, source_nodes):

    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "dist_type" specifies whether to calculate "min", "mean", "median", or "max" distance values 
    from a region-of-interest. If running only on single node, defaults to "min".
    """
    
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    translated_source_nodes = translate_src(source_nodes, cortex)
    data = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)        
    dist = recort(data, surf, cortex)
    del data

    return dist

def zone_calc(surf, cortex, src):
    """
    Calculate closest nodes to each source node using exact geodesic distance along the cortical surface.
    """
    
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    dist_vals = np.zeros((len(source_nodes), len(cortex_vertices)))

    for x in range(len(source_nodes)):
        
        translated_source_nodes = translate_src(source_nodes[x], cortex)
        dist_vals[x, :] = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)

    data = np.argsort(dist_vals, axis=0)[0, :] + 1

    zone = recort(data, surf, cortex)

    del data

    return zone

def dist_calc_matrix(surf, cortex, labels, exceptions = ['Unknown', 'Medial_wall'], verbose = True):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "labels" specifies the freesurfer label file to use. All values will be used other than those
    specified in "exceptions" (default: 'Unknown' and 'Medial_Wall').
    
    returns: 
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """
    
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    
    # remove exceptions from label list:
    label_list = sd.load.get_freesurfer_label(labels, verbose = False)
    rs = np.where([a not in exceptions for a in label_list])[0]
    rois = [label_list[r] for r in rs]
    if verbose:
        print "# of regions: " + str(len(rois))

    # calculate distance from each region to all nodes:
    dist_roi = []
    for roi in rois:
        source_nodes = sd.load.load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_roi.append(gdist.compute_gdist(cortex_vertices, cortex_triangles, 
                                                source_indices = translated_source_nodes))
        if verbose:
            print roi            
    dist_roi = np.array(dist_roi)

    # Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes = sd.load.load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    return dist_mat, rois
