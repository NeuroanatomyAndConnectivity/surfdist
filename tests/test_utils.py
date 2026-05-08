import numpy as np

from surfdist.utils import (
    surf_keep_cortex,
    triangles_keep_cortex,
    translate_src,
    recort,
    find_node_match,
)


def test_surf_keep_cortex_full(grid_surface, full_cortex):
    verts, tris, N = grid_surface
    cv, ct = surf_keep_cortex((verts, tris), full_cortex)
    assert cv.shape == (N * N, 3)
    assert ct.shape == tris.shape
    np.testing.assert_array_equal(ct, tris)


def test_surf_keep_cortex_partial(grid_surface, partial_cortex):
    verts, tris, _ = grid_surface
    cv, ct = surf_keep_cortex((verts, tris), partial_cortex)
    assert cv.shape == (len(partial_cortex), 3)
    assert ct.shape[0] < tris.shape[0]
    assert ct.max() < len(partial_cortex)
    assert ct.min() >= 0


def test_triangles_keep_cortex_drops_medial(grid_surface, partial_cortex):
    _, tris, _ = grid_surface
    ct = triangles_keep_cortex(tris, partial_cortex)
    medial = np.setdiff1d(np.arange(16), partial_cortex)
    # No reindexed triangle can reference a medial-wall vertex
    assert not np.any(np.isin(partial_cortex[ct], medial))


def test_translate_src_identity(full_cortex):
    src = np.array([0, 5, 10], dtype=np.int32)
    out = translate_src(src, full_cortex)
    np.testing.assert_array_equal(np.sort(out), np.sort(src))


def test_translate_src_with_medial(partial_cortex):
    src = np.array([0, 5], dtype=np.int32)
    out = translate_src(src, partial_cortex)
    expected = np.array([
        np.where(partial_cortex == 0)[0][0],
        np.where(partial_cortex == 5)[0][0],
    ])
    np.testing.assert_array_equal(np.sort(out), np.sort(expected))


def test_recort_zeros_outside_cortex(grid_surface, partial_cortex):
    verts, tris, N = grid_surface
    data = np.ones(len(partial_cortex), dtype=np.float64)
    full = recort(data, (verts, tris), partial_cortex)
    assert full.shape == (N * N,)
    np.testing.assert_array_equal(full[partial_cortex], 1)
    medial = np.setdiff1d(np.arange(N * N), partial_cortex)
    np.testing.assert_array_equal(full[medial], 0)


def test_find_node_match_identity():
    verts = np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [1., 1., 0.],
    ])
    idx, _ = find_node_match(verts, verts)
    np.testing.assert_array_equal(np.sort(idx), np.arange(len(verts)))
