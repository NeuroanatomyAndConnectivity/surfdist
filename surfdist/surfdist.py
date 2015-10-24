import numpy as np


def surf_keep_cortex(surf, cortex):
    """
    Remove medial wall from cortical surface to ensure that shortest paths are only calculated through the cortex.
    """
    vertices = np.array(surf[0][cortex], dtype=np.float64)

    keep = np.zeros(len(surf[1]))
    for i in [0, 1, 2]:
        keep += np.array([item in cortex for item in surf[1][:, i]])
    ind = np.where(keep == 3)
    triangles = np.array(surf[1][ind], dtype=np.int32)
    for c, i in enumerate(cortex):
        triangles[np.where(triangles == i)] = c

    return vertices, triangles


def recort(input_data, surf, cortex):
    """
    Return data values to space of full cortex (including medial wall), with medial wall equal to zero.
    """
    data = np.zeros(len(surf[0]))
    data[cortex] = input_data
    return data


def viz(surf, data):
    """
    Visualize results on cortical surace using mayavi.
    """
    vertices = np.array(surf[0], dtype=np.float64)
    triangles = np.array(surf[1], dtype=np.int32)
    from mayavi import mlab
    mlab.triangular_mesh(
        vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles, scalars=data)


def dist_calc(surf, cortex, src):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    """
    import gdist

    vertices, triangles = surf_keep_cortex(surf, cortex)
    dist = gdist.compute_gdist(vertices, triangles, source_indices=src)

    return dist


def zone_calc(surf, cortex, src):
    """
    Calculate closest nodes to each source node using exact geodesic distance along the cortical surface.
    """
    import gdist

    vertices, triangles = surf_keep_cortex(surf, cortex)

    dist_vals = []
    for x in range(len(src)):
        dist_vals.append(
            gdist.compute_gdist(vertices, triangles, source_indices=src[x]))

    zone = argsort(dist_vals)

    return zone
