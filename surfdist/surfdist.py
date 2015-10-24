import numpy as np


def dist_calc(surf, cortex, src):
    """
    Calculate exact geodesic distance along cortical surface from set of source nodes.
    """
    import gdist
    from utils import surf_keep_cortex, recort

    vertices, triangles = surf_keep_cortex(surf, cortex)
    data = gdist.compute_gdist(vertices, triangles, source_indices=src)
    dist = recort(data, surf, cortex)
    del data

    return dist


def zone_calc(surf, cortex, src):
    """
    Calculate closest nodes to each source node using exact geodesic distance along the cortical surface.
    """
    import gdist
    from utils import surf_keep_cortex, recort

    vertices, triangles = surf_keep_cortex(surf, cortex)

    dist_vals = np.zeros((len(src), len(vertices)))
    for x in range(len(src)):
        dist_vals[x, :] = gdist.compute_gdist(vertices, triangles,
                                              source_indices=np.array([src[x]], dtype=np.int32))

    data = np.argsort(dist_vals, axis=0)[-1, :] + 1
    zone = recort(data, surf, cortex)
    del data

    return zone
