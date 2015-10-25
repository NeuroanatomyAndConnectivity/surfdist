import numpy as np

def viz(surf, data):
    """
    Visualize results on cortical surace using mayavi.
    """
    vertices = np.array(surf[0], dtype=np.float64)
    triangles = np.array(surf[1], dtype=np.int32)

    from mayavi import mlab
    mlab.triangular_mesh(
        vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles, scalars=data)
