import nibabel as nib
import numpy as np
import os
from mayavi import mlab
import surfdist.surfdist as sd

base_dir = '/Applications/freesurfer/subjects/fsaverage4/'

surf = nib.freesurfer.read_geometry(os.path.join(base_dir, 'surf/lh.pial'))
cort = np.sort(nib.freesurfer.read_label(os.path.join(base_dir, 'label/lh.cortex.label')))
annot = nib.freesurfer.read_annot(os.path.join(base_dir, 'label/lh.aparc.a2009s.annot'))

# Calculate geodesic distance from central sulcus:
label_text = 'S_central'
label_ind = annot[2].index(label_text)
labels = np.where(np.in1d(annot[0],label_ind))
labels = np.where(np.in1d(cort, labels))[0]
src = np.array(labels, dtype=np.int32)

dist = sd.dist_calc(surf, cort, src)
data = sd.recort(dist, surf, cort)
sd.viz(surf, data)
