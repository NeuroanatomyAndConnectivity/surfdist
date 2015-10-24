# -*- coding: utf-8 -*-

import numpy as np

def surf_keep_cortex(surf, cortex):
        
    vertices = np.array(surf[0][cortex], dtype=np.float64)

    keep = np.zeros(len(surf[1]))
    for i in [0,1,2]:
        keep += np.array([item in cortex for item in surf[1][:,i]])
    ind = np.where(keep == 3)
    triangles = np.array(surf[1][ind], dtype=np.int32)
    for c,i in enumerate(cortex):
        triangles[np.where(triangles == i)] = c
    
    return vertices, triangles

def recort(input_data, surf, cortex):
    data = np.zeros(len(surf[0]))
    data[cortex] = input_data
    return data
    
def viz(surf, data):
    vertices = np.array(surf[0], dtype=np.float64)
    triangles = np.array(surf[1], dtype=np.int32)
    from mayavi import mlab
    mlab.triangular_mesh(vertices[:,0], vertices[:,1], vertices[:,2], triangles, scalars=data)    

def distcalc(surf, cortex, src):
    
    import gdist

    vertices, triangles = surf_keep_cortex(surf, cortex)
    dist = gdist.compute_gdist(vertices, triangles, source_indices = src)

    return dist
    