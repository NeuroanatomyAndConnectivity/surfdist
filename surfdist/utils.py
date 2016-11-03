import numpy as np


def surf_keep_cortex(surf, cortex):
    """
    Remove medial wall from cortical surface to ensure that shortest paths are only calculated through the cortex.

    Inputs
    -------
    surf : Tuple containing two numpy arrays of shape (n_nodes,3). Each node of the first array specifies the x, y, z
           coordinates one node of the surface mesh. Each node of the second array specifies the indices of the three
           nodes building one triangle of the surface mesh.
           (e.g. the output from nibabel.freesurfer.io.read_geometry)
    cortex : Array with indices of vertices included in within the cortex.
             (e.g. the output from nibabel.freesurfer.io.read_label)
    """

    # split surface into coordinates and faces
    vertices, triangles = surf

    # keep only the vertices within the cortex label
    cortex_vertices = np.array(vertices[cortex], dtype=np.float64)

    # remove faces with nodes not in the newly defined list of vertices
    keep = np.zeros(len(triangles))
    # for or each face/triangle keep only those that contain 3 nodes within the list of cortex nodes
    for i in [0, 1, 2]:
        keep += np.array([item in cortex for item in triangles[:, i]])
    ind = np.where(keep == 3)
    cortex_triangles_old = np.array(triangles[ind], dtype=np.int32)
    cortex_triangles = np.array(triangles[ind], dtype=np.int32)
    # reassign node index before outputting triangles
    for c, i in enumerate(cortex):
        cortex_triangles[np.where(cortex_triangles_old == i)] = c

    return cortex_vertices, cortex_triangles


def translate_src(src, cortex):
    """
    Convert source nodes to new surface (without medial wall).
    """
    src_new = np.array(np.where(np.in1d(cortex, src))[0], dtype=np.int32)

    return src_new


def recort(input_data, surf, cortex):
    """
    Return data values to space of full cortex (including medial wall), with medial wall equal to zero.
    """
    data = np.zeros(len(surf[0]))
    data[cortex] = input_data
    return data


def find_node_match(simple_vertices, complex_vertices):
    '''
    Thanks to juhuntenburg.
    Functions taken from https://github.com/juhuntenburg/brainsurfacescripts

    Finds those points on the complex mesh that correspoind best to the
    simple mesh while forcing a one-to-one mapping.
    '''

    import scipy.spatial

    # make array for writing in final voronoi seed indices
    voronoi_seed_idx = np.zeros((simple_vertices.shape[0],), dtype='int64')-1
    missing = np.where(voronoi_seed_idx == -1)[0].shape[0]
    mapping_single = np.zeros_like(voronoi_seed_idx)

    neighbours = 0
    col = 0

    while missing != 0:

        neighbours += 100
        # find nearest neighbours
        inaccuracy, mapping = scipy.spatial.KDTree(
            complex_vertices).query(simple_vertices, k=neighbours)
        # go through columns of nearest neighbours until unique mapping is
        # achieved, if not before end of neighbours, extend number of
        # neighbours
        while col < neighbours:
            # find all missing voronoi seed indices
            missing_idx = np.where(voronoi_seed_idx == -1)[0]
            missing = missing_idx.shape[0]
            if missing == 0:
                break
            else:
                # for missing entries fill in next neighbour
                mapping_single[missing_idx] = np.copy(
                    mapping[missing_idx, col])
                # find unique values in mapping_single
                unique, double_idx = np.unique(
                    mapping_single, return_inverse=True)
                # empty voronoi seed index
                voronoi_seed_idx = np.zeros(
                    (simple_vertices.shape[0],), dtype='int64')-1
                # fill voronoi seed idx with unique values
                for u in range(unique.shape[0]):
                    # find the indices of this value in mapping
                    entries = np.where(double_idx == u)[0]
                    # set the first entry to the value
                    voronoi_seed_idx[entries[0]] = unique[u]
                # go to next column
                col += 1

    return voronoi_seed_idx, inaccuracy
