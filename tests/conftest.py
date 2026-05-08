import numpy as np
import pytest


def _build_grid(N):
    xs, ys = np.meshgrid(np.arange(N), np.arange(N))
    verts = np.column_stack([
        xs.ravel().astype(np.float64),
        ys.ravel().astype(np.float64),
        np.zeros(N * N, dtype=np.float64),
    ])
    tris = []
    for r in range(N - 1):
        for c in range(N - 1):
            i = r * N + c
            tris.append([i, i + 1, i + N])
            tris.append([i + 1, i + N + 1, i + N])
    return verts, np.array(tris, dtype=np.int32)


@pytest.fixture
def grid_surface():
    """A 4x4 planar triangular mesh.

    On a flat polyhedral mesh, the exact geodesic distance from gdist
    equals straight-line Euclidean distance between vertices.
    """
    N = 4
    verts, tris = _build_grid(N)
    return verts, tris, N


@pytest.fixture
def full_cortex(grid_surface):
    _, _, N = grid_surface
    return np.arange(N * N, dtype=np.int32)


@pytest.fixture
def partial_cortex(grid_surface):
    """Cortex excluding the rightmost column (a synthetic 'medial wall')."""
    _, _, N = grid_surface
    medial = np.arange(N - 1, N * N, N)
    return np.array(
        [i for i in range(N * N) if i not in medial], dtype=np.int32
    )
