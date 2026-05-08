import nibabel as nib
import numpy as np
import pytest

from surfdist.analysis import dist_calc, dist_calc_matrix


def test_dist_calc_planar_equals_euclidean(grid_surface, full_cortex):
    """On a flat polyhedral mesh, exact geodesic == Euclidean."""
    verts, tris, _ = grid_surface
    src = np.array([0], dtype=np.int32)
    dist = dist_calc((verts, tris), full_cortex, src)
    expected = np.linalg.norm(verts - verts[0], axis=1)
    np.testing.assert_allclose(dist, expected, atol=1e-10)


def test_dist_calc_zero_at_source(grid_surface, full_cortex):
    verts, tris, _ = grid_surface
    src = np.array([7], dtype=np.int32)
    dist = dist_calc((verts, tris), full_cortex, src)
    assert dist[7] == 0.0
    assert (dist >= 0).all()


def test_dist_calc_recorts_to_full_surface(grid_surface, partial_cortex):
    """Returned array spans the full surface; medial-wall vertices read 0."""
    verts, tris, N = grid_surface
    # Use a source that is in the partial cortex
    src = np.array([0], dtype=np.int32)
    dist = dist_calc((verts, tris), partial_cortex, src)
    assert dist.shape == (N * N,)
    medial = np.setdiff1d(np.arange(N * N), partial_cortex)
    np.testing.assert_array_equal(dist[medial], 0)


def _write_quadrant_annot(path, N=4):
    """Split the 4x4 grid into 4 row-quadrants A/B/C/D."""
    n = N * N
    labels = np.repeat(np.arange(N), N).astype(np.int32)
    ctab = np.array([
        [10, 10, 10, 0, 0],
        [20, 20, 20, 0, 1],
        [30, 30, 30, 0, 2],
        [40, 40, 40, 0, 3],
    ], dtype=np.int32)
    names = [b'A', b'B', b'C', b'D']
    nib.freesurfer.io.write_annot(str(path), labels, ctab, names)
    return labels


@pytest.mark.parametrize('summary', ['min', 'mean', 'median', 'max'])
def test_dist_calc_matrix_summary_kwarg(
    grid_surface, full_cortex, tmp_path, summary
):
    """Exercises the summary= kwarg added in PR #19 (kaurao)."""
    verts, tris, _ = grid_surface
    annot = tmp_path / 'grid.annot'
    _write_quadrant_annot(annot)

    mat, rois = dist_calc_matrix(
        (verts, tris), full_cortex, str(annot),
        exceptions=[], summary=summary, verbose=False,
    )
    assert mat.shape == (4, 4)
    assert (mat >= 0).all()
    # Diagonal is always 0 regardless of summary: every region's source set
    # contains its own vertices, so distance from the set to itself is 0.
    np.testing.assert_array_equal(np.diag(mat), 0)


def test_dist_calc_matrix_min_is_symmetric_for_disjoint_rows(
    grid_surface, full_cortex, tmp_path
):
    verts, tris, _ = grid_surface
    annot = tmp_path / 'grid.annot'
    _write_quadrant_annot(annot)
    mat, _ = dist_calc_matrix(
        (verts, tris), full_cortex, str(annot),
        exceptions=[], summary='min', verbose=False,
    )
    np.testing.assert_allclose(mat, mat.T, atol=1e-10)


def test_dist_calc_matrix_excludes_medial_wall(
    grid_surface, full_cortex, tmp_path
):
    verts, tris, _ = grid_surface
    annot = tmp_path / 'grid.annot'
    _write_quadrant_annot(annot)
    # Treat 'A' as a synthetic medial wall and verify it's excluded
    mat, rois = dist_calc_matrix(
        (verts, tris), full_cortex, str(annot),
        exceptions=['A'], summary='min', verbose=False,
    )
    assert 'A' not in rois
    assert mat.shape == (3, 3)
