import nibabel as nib
import numpy as np
import os
import surfdist as sd

# calculate and display distance from central sulcus at each node:
base_dir = '/Applications/freesurfer/subjects/bert/'
surf = nib.freesurfer.read_geometry(os.path.join(base_dir, 'surf/lh.pial'))
cort = np.sort(nib.freesurfer.read_label(os.path.join(base_dir, 'label/lh.cortex.label')))
src  = sd.load.load_freesurfer_label(os.path.join(base_dir, 'label/lh.aparc.a2009s.annot'), 'S_central')

dist = sd.surfdist.dist_calc(surf, cort, src)
sd.viz.viz(surf, dist)

# Calculate distances on native surface and display on fsaverage
fsa4 = nib.freesurfer.read_geometry('/Applications/freesurfer/subjects/fsaverage4/surf/lh.sphere.reg')[0]
native = nib.freesurfer.read_geometry('/Applications/freesurfer/subjects/bert/surf/lh.sphere.reg')[0]
idx_fsa4_to_native = sd.utils.find_node_match(fsa4, native)[0]

surf_fsa4 = nib.freesurfer.read_geometry('/Applications/freesurfer/subjects/fsaverage4/surf/lh.pial')
sd.viz.viz(surf_fsa4, dist[idx_fsa4_to_native])


src = [10001, 50]
dist_vals = np.zeros((len(src),len(vertices)))
for x in range(len(src)):
    dist_vals[x,:] = gdist.compute_gdist(vertices, triangles, 
                        source_indices=np.array([src[x]], dtype=np.int32))

zone = np.argsort(dist_vals, axis=0)[-1,:]