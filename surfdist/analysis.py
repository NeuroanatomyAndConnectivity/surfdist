import gdist
import numpy as np
from surfdist.utils import surf_keep_cortex, translate_src, recort
from surfdist import load


def dist_calc(surf, cortex, source_nodes, max_distance=None, recortex=True):
    """
    Calculate exact geodesic distance along the cortical surface from a set
    of source nodes to every other node.

    Inputs
    -------
    surf : (vertices, triangles) tuple as returned by
           ``nibabel.freesurfer.read_geometry`` or
           ``(gii.darrays[0].data, gii.darrays[1].data)`` for a gifti.
    cortex : array of cortex vertex indices (medial wall excluded).
    source_nodes : indices of vertices in the source ROI.
    max_distance : float or None, optional
        If provided, geodesic propagation stops at this distance.
        Vertices unreachable within the threshold are returned as a
        large sentinel value by gdist; callers can mask them with
        ``dist > max_distance``. Default ``None`` (no limit).
    recortex : bool, default True
        If True, the returned array spans the full surface (medial-wall
        vertices read 0). If False, the returned array spans only the
        cortex (in cortex-vertex order); useful when ``dist_calc`` is
        called as an inner loop and recortexing is wasted work.

    Returns
    -------
    dist : ndarray
    """
    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)
    translated_source_nodes = translate_src(source_nodes, cortex)
    if max_distance is None:
        data = gdist.compute_gdist(
            cortex_vertices, cortex_triangles,
            source_indices=translated_source_nodes,
        )
    else:
        data = gdist.compute_gdist(
            cortex_vertices, cortex_triangles,
            source_indices=translated_source_nodes,
            max_distance=float(max_distance),
        )
    if recortex:
        return recort(data, surf, cortex)
    return data


def calc_roi_dist(surf, cortex, source_nodes, target_nodes, summary='min'):
    """
    Geodesic distance from source ROI X to target ROI Y, summarized to a scalar.

    Inputs
    -------
    surf : (vertices, triangles) tuple. See dist_calc.
    cortex : array of cortex vertex indices.
    source_nodes : indices of vertices in ROI X (the propagation source).
    target_nodes : indices of vertices in ROI Y (where distances are sampled).
    summary : str, one of 'min', 'mean', 'median', 'max'. How distances
              from each target vertex back to the source set are reduced
              into a single ROI-to-ROI value. Default 'min' (the
              conventional shortest distance between regions).

    Returns
    -------
    roi_dist : float
    """
    dists = dist_calc(surf, cortex, source_nodes)
    dists_to_target = dists[np.asarray(target_nodes).ravel()]
    if summary == 'min':
        return float(np.min(dists_to_target))
    if summary == 'mean':
        return float(np.mean(dists_to_target))
    if summary == 'median':
        return float(np.median(dists_to_target))
    if summary == 'max':
        return float(np.max(dists_to_target))
    raise ValueError(
        f"unknown summary {summary!r}; expected one of "
        "['min', 'mean', 'median', 'max']"
    )


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


def dist_calc_matrix(surf, cortex, labels, exceptions = ['Unknown', 'Medial_wall'], summary = 'min', verbose = True):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    "labels" specifies the freesurfer label file to use. All values will be used other than those
    specified in "exceptions" (default: 'Unknown' and 'Medial_Wall').
    summary defines how the distances are summarized with suppoted values: 'min', 'mean', 'median', 'max'

    returns:
      dist_mat: symmetrical nxn matrix of minimum distance between pairs of labels
      rois: label names in order of n
    """

    cortex_vertices, cortex_triangles = surf_keep_cortex(surf, cortex)

    # remove exceptions from label list:
    label_list = load.get_freesurfer_label(labels, verbose = False)
    # nibabel returns names as bytes; normalize so string exceptions work
    label_list = [n.decode('utf-8') if isinstance(n, (bytes, bytearray)) else n
                  for n in label_list]
    rs = np.where([a not in exceptions for a in label_list])[0]
    rois = [label_list[r] for r in rs]
    if verbose:
        print("# of regions: " + str(len(rois)))

    # calculate distance from each region to all nodes:
    dist_roi = []
    for roi in rois:
        source_nodes = load.load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        dist_roi.append(gdist.compute_gdist(cortex_vertices, cortex_triangles,
                                                source_indices = translated_source_nodes))
        if verbose:
            print(roi)
    dist_roi = np.array(dist_roi)

    # Calculate min distance per region:
    dist_mat = []
    for roi in rois:
        source_nodes = load.load_freesurfer_label(labels, roi)
        translated_source_nodes = translate_src(source_nodes, cortex)
        if summary == 'min':
            dist_mat.append(np.min(dist_roi[:,translated_source_nodes], axis = 1))
        elif summary == 'mean':
            dist_mat.append(np.mean(dist_roi[:,translated_source_nodes], axis = 1))
        elif summary == 'median':
            dist_mat.append(np.median(dist_roi[:,translated_source_nodes], axis = 1))
        elif summary == 'max':
            dist_mat.append(np.max(dist_roi[:,translated_source_nodes], axis = 1))
        else:
            raise(f'undefined summary: {summary}')
    dist_mat = np.array(dist_mat)

    return dist_mat, rois
