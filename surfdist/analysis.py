import gdist
import numpy as np
import nibabel as nib
from surfdist.utils import surf_keep_cortex, translate_src, recort
import surfdist as sd
from surfdist.load import load_cifti_labels

def dist_calc(surf, cortex, source_nodes,gifti=False):

    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "dist_type" specifies whether to calculate "min", "mean", "median", or "max" distance values
    from a region-of-interest. If running only on single node, defaults to "min".
    """

    if gifti !=False:
        gii=nib.load(surf)
        surf=(gii.darrays[0].data,gii.darrays[1].data)

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
    tmp=[i.decode('utf-8') for i in label_list]
    label_list=tmp
    del tmp
    rois=list(set(label_list)-set(exceptions))
    print(len(label_list))
    print(len(exceptions))
    print(len(rois))

    if verbose:
        print("# of regions: " + str(len(rois)))

    ###calculate distance from each region to all nodes:
    dist_roi = []
    for roi in rois:
        source_nodes = sd.load.load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_roi.append(gdist.compute_gdist(cortex_vertices, cortex_triangles,
                                                source_indices = translated_source_nodes))
        if verbose:
            print(roi)
    dist_roi = np.array(dist_roi)

    ###Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes = sd.load.load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    return dist_mat,rois

def dist_calc_matrixCifti(surfGii,cifti,hemi,verbose=True):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes defined labels in a cifti file.
    "labels" specifies the freesurfer label file to use. All labels in the cifti file except the medial wall will be used.
    returns:
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """
    
    gii=nib.load(surfGii)
    surf=(gii.darrays[0].data,gii.darrays[1].data)
    
    Llabels,Rlabels=load_cifti_labels(cifti)
    ctx=np.array(range(len(surf[0])))

    ### get cortex vertices and remove medial wall from label list 
    if hemi=='L':
        medialWall=Llabels['???']
        del Llabels['???']
        labels=Llabels
    elif hemi=='R':
        medialWall=Rlabels['???']
        del Rlabels['???']
        labels=Rlabels
    cortex=np.delete(ctx,medialWall)
    
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

#     # remove exceptions from label list:
    rois=labels
    if verbose:
        print("# of regions: " + str(len(rois)))

    ###calculate distance from each region to all nodes:
    dist_roi = []
    for roi in rois:
        source_nodes =rois[roi]
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_roi.append(gdist.compute_gdist(cortex_vertices, cortex_triangles,
                                                source_indices = translated_source_nodes))
        if verbose:
            print(roi)
    dist_roi = np.array(dist_roi)

    ###Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes =rois[roi]
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    return dist_mat,rois