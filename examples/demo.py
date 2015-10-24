import nibabel as nib
import numpy as np
import os
import surfdist as sd
from surfdist import viz, load, utils, surfdist

# calculate and display distance from central sulcus at each node:
base_dir = '/Applications/freesurfer/subjects/bert/'
surf = nib.freesurfer.read_geometry(os.path.join(base_dir, 'surf/lh.pial'))
cort = np.sort(nib.freesurfer.read_label(os.path.join(base_dir, 'label/lh.cortex.label')))

# load central sulcus nodes
src  = sd.load.load_freesurfer_label(os.path.join(base_dir, 'label/lh.aparc.a2009s.annot'), 'S_central')

# calculate distance
dist = sd.surfdist.dist_calc(surf, cort, src)

# visualize
sd.viz.viz(surf, dist)


# Calculate distances on native surface and display on fsaverage
fsa4 = nib.freesurfer.read_geometry('/Applications/freesurfer/subjects/fsaverage4/surf/lh.sphere.reg')[0]
native = nib.freesurfer.read_geometry('/Applications/freesurfer/subjects/bert/surf/lh.sphere.reg')[0]
idx_fsa4_to_native = sd.utils.find_node_match(fsa4, native)[0]

surf_fsa4 = nib.freesurfer.read_geometry('/Applications/freesurfer/subjects/fsaverage4/surf/lh.pial')
sd.viz.viz(surf_fsa4, dist[idx_fsa4_to_native])
