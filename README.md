surfdist
========
[![Build Status](https://travis-ci.com/NeuroanatomyAndConnectivity/surfdist.svg?branch=master)](https://travis-ci.com/NeuroanatomyAndConnectivity/surfdist)

Calculate the exact geodesic distance on a triangular surface mesh using the [gdist package](https://pypi.python.org/pypi/gdist/), which is based on the [c++ library](https://code.google.com/p/geodesic/).

Installation
------------

    pip install surfdist

Example
-------
Freesurfer files:

    import nibabel as nib
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import surfdist as sd
    from surfdist import viz, load, utils, analysis

    # calculate and display distance from central sulcus at each node:
    cmap = 'coolwarm'
    base_dir = '/Applications/freesurfer/subjects/'
    surf = nib.freesurfer.read_geometry(os.path.join(base_dir, 'bert/surf/lh.pial'))
    cort = np.sort(nib.freesurfer.read_label(os.path.join(base_dir, 'bert/label/lh.cortex.label')))
    sulc = nib.freesurfer.read_morph_data(os.path.join(base_dir, 'bert/surf/lh.sulc'))

    # load central sulcus nodes
    src  = sd.load.load_freesurfer_label(os.path.join(base_dir, 'bert/label/lh.aparc.a2009s.annot'), 'S_central', cort)

    # calculate distance
    dist = sd.analysis.dist_calc(surf, cort, src)

    # visualize
    plot_med = sd.viz.viz(surf[0], surf[1], dist, bg_map=sulc, bg_on_stat=True, cmap=cmap)
    plot_lat = sd.viz.viz(surf[0], surf[1], dist, azim=180, bg_map=sulc, bg_on_stat=True, cmap=cmap)

    # Calculate distances on native surface and display on fsaverage
    fsa4 = nib.freesurfer.read_geometry(os.path.join(base_dir,'fsaverage4/surf/lh.sphere.reg'))[0]
    fsa4_sulc=nib.freesurfer.read_morph_data(os.path.join(base_dir, 'fsaverage4/surf/lh.sulc'))
    native = nib.freesurfer.read_geometry(os.path.join(base_dir, 'bert/surf/lh.sphere.reg'))[0]
    idx_fsa4_to_native = sd.utils.find_node_match(fsa4, native)[0]

    surf_fsa4 = nib.freesurfer.read_geometry(os.path.join(base_dir, 'fsaverage4/surf/lh.pial'))
    plot_fsa4_med = sd.viz.viz(surf_fsa4[0], surf_fsa4[1], dist[idx_fsa4_to_native], bg_map=fsa4_sulc, bg_on_stat=True, cmap=cmap)
    plot_fsa4_lat = sd.viz.viz(surf_fsa4[0], surf_fsa4[1], dist[idx_fsa4_to_native], azim=180, bg_map=fsa4_sulc, bg_on_stat=True, cmap=cmap)

    plt.show()

Gifti files:

    import surfdist as sd
    from surfdist import viz, load, utils, analysis

    surf_labels = nib.load("fsLR.32k.L.label.gii")
    # pick only the vertices of the cortex, excluding the medial wall
    cortex = np.where(labels.darrays[0].data != 0)[0]

    surfL = nib.load("sub-1_hemi-L_inflated.32k_fs_LR.surf.gii")
    nodes = surfL.agg_data('NIFTI_INTENT_POINTSET')
    triangles = surfL.agg_data('NIFTI_INTENT_TRIANGLE')
    surf = (nodes, triangles)

    destrieux = nib.load("destrieux-labels_den-32k_hemi-L.label.gii").darrays[0].data

    # pick only the vertices of A1 and angular gyrus.
    a1_vrtx = np.where(destrieux == 32)[0]
    angG_vrtx = np.where(destrieux == 24)[0]

    # calculate distances from A1 to the rest of the vertices
    all_dist = analysis.dist_calc(surf, cortex, a1_vrtx)
    # calculate the shortest distance from A1 to angular gyrus
    dist_min = anaysis.dist_calc(surf, cortex, a1_vrtx, angG_vrtx)

