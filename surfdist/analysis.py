import gdist
import numpy as np
import nibabel as nib
from surfdist.utils import surf_keep_cortex, translate_src, recort, roiDistance,AnatomyInputParser
import surfdist as sd
from surfdist.load import load_cifti_labels,load_freesurfer_label,get_freesurfer_label, load_gifti_labels
#### multiprocessing is used in cifti and gifti distance matrix calculations 
from multiprocessing.pool import Pool as ProcessPool
import multiprocessing


def dist_calc(surf, cortex, source_nodes, recortex=True, gifti=False):

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
    dist = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)
    if recortex==True:
        dist = recort(dist, surf, cortex)

    return dist

def zone_calc(surf, cortex, source_nodes,gifti=False):
    """
    Calculate closest nodes to each source node using exact geodesic distance along the cortical surface.
    """
    if gifti!=False:
        gii=nib.load(surf)
        surf=(gii.darrays[0].data,gii.darrays[1].data)
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    dist_vals = np.zeros((len(source_nodes), len(cortex_vertices)))
    for x in range(len(source_nodes)):
        translated_source_nodes = translate_src(source_nodes[x], cortex)
        dist_vals[x, :] = gdist.compute_gdist(cortex_vertices, cortex_triangles, source_indices = translated_source_nodes)

    data = np.argsort(dist_vals, axis=0)[0, :] + 1

    zone = recort(data, surf, cortex)

    del data

    return zone


def dist_calc_matrixFS(surf, cortex, labels, exceptions = ['Unknown', 'Medial_wall'], verbose = True):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    This function is speciffic to FreeSurfer inputs.
    "labels" specifies the freesurfer label file to use. All values will be used other than those
    specified in "exceptions" (default: 'Unknown' and 'Medial_Wall').
    example: dist_calc_matrixFS(lh.pial,lh.cortex.label,lh.aparc.annot)
    returns:
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    # remove exceptions from label list:
    label_list = get_freesurfer_label(labels, verbose = False)
    tmp=[i.decode('utf-8') for i in label_list]
    label_list=tmp
    del tmp
    label_list=list(set(label_list)-set(exceptions))


    if verbose:
        print("# of regions: " + str(len(label_list)))

    
    rois=[load_freesurfer_label(labels, roi) for roi in label_list]
    rois=dict(zip(label_list,rois))
   

    nodes= list(rois.values())
    params=[[nodes[i],cortex,cortex_vertices,cortex_triangles] for i in range(len(nodes))]
    
    ###calculate distance from each region to all nodes:
    with ProcessPool(processes=n_cpus-1) as pool:
        dist_roi=pool.starmap(roiDistance,params)
    dist_roi=np.array(dist_roi)

    ###Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes = load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    return dist_mat,list(rois.keys())

def dist_calc_matrixCifti(surfGii, cifti, hemi):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes defined labels in a cifti file.
    "labels" specifies the freesurfer label file to use. All labels in the cifti file except the medial wall will be used.
    example dist_calc_matrixCifti('anat.surf.gii','cifti.dlabel.nii','L')
    returns:
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """

    gii = nib.load(surfGii)
    surf = (gii.darrays[0].data, gii.darrays[1].data)

    Llabels, Rlabels = load_cifti_labels(cifti)
    ctx = np.array(range(len(surf[0])))

    ### get cortex vertices and remove medial wall from label list 
    if hemi == 'L':
        medialWall = Llabels['???']
        del Llabels['???']
        labels = Llabels
    elif hemi == 'R':
        medialWall = Rlabels['???']
        del Rlabels['???']
        labels = Rlabels
    cortex = np.delete(ctx, medialWall)

    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    # remove exceptions from label list:
    rois=labels
    nodes= list(rois.values())
    params=[[nodes[i],cortex,cortex_vertices,cortex_triangles] for i in range(len(nodes))]
    
    with ProcessPool(processes=n_cpus-1) as pool:
        dist_roi=pool.starmap(roiDistance,params)
    dist_roi=np.array(dist_roi)
    
     ###Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes=rois[roi]
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    return dist_mat,list(rois.keys())

def dist_calc_matrixGifti(AnatSurf, giftiLabel,exceptions,n_cpus=1):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes defined labels in a cifti file.
    "labels" specifies the freesurfer label file to use. All labels in the cifti file except the medial wall will be used.
    
    YOU MUST PROVIDE A LIST OF EXCEPTIONS FOR AREAS NOT INCLUDED IN THE CALCUALATION BASED ON YOUR GIFTI FILE 
    example: dist_calc_matrixGifti('anat.surf.gii','labels.label.gii',['???','L_Medial_wall'])
    
    returns:
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """
    print(f'using {n_cpus} cpus')
    surf=AnatomyInputParser(AnatSurf)

    labels= load_gifti_labels(giftiLabel)
    medialWall = labels['???']
  
    ctx = np.array(range(len(surf[0])))
    cortex = np.delete(ctx, medialWall)
    
    for i in exceptions:
        del labels[i]
    
    # remove exceptions from label list:
    rois=labels
    nodes= list(rois.values())
    params=[[surf,cortex,nodes[i],'recort=False'] for i in range(len(nodes))]
    
    with ProcessPool(processes=n_cpus) as pool:
        dist_roi=pool.starmap(dist_calc,params)
    dist_roi=np.array(dist_roi)
    
     ###Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes=rois[roi]
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
    dist_mat = np.array(dist_mat)

    return dist_mat,list(rois.keys())