import nibabel as nib
import numpy as np
import pytest
from nibabel import cifti2
from nibabel.gifti import GiftiDataArray, GiftiImage, GiftiLabel, GiftiLabelTable

from surfdist.load import (
    get_freesurfer_label,
    load_cifti_labels,
    load_freesurfer_label,
    load_gifti_labels,
    load_labels,
    load_surface,
)


def _write_fake_annot(path, n_vertices=16):
    """Two-label annot: id 0 ('Medial_wall') over verts 0..3, id 1 ('Cortex') over 4..15."""
    labels = np.zeros(n_vertices, dtype=np.int32)
    labels[4:] = 1
    ctab = np.array([
        [25, 5, 25, 0, 0],
        [50, 50, 50, 0, 1],
    ], dtype=np.int32)
    names = [b'Medial_wall', b'Cortex']
    nib.freesurfer.io.write_annot(str(path), labels, ctab, names)


def test_get_freesurfer_label(tmp_path):
    annot = tmp_path / 'lh.test.annot'
    _write_fake_annot(annot)
    names = get_freesurfer_label(str(annot), verbose=False)
    decoded = [n.decode('utf-8') if isinstance(n, bytes) else n for n in names]
    assert 'Cortex' in decoded
    assert 'Medial_wall' in decoded


def test_load_freesurfer_label_finds_cortex(tmp_path):
    annot = tmp_path / 'lh.test.annot'
    _write_fake_annot(annot)
    nodes = load_freesurfer_label(str(annot), 'Cortex')
    expected = np.arange(4, 16)
    assert nodes.dtype == np.int32
    np.testing.assert_array_equal(np.sort(nodes.ravel()), expected)


def test_load_freesurfer_label_accepts_bytes(tmp_path):
    annot = tmp_path / 'lh.test.annot'
    _write_fake_annot(annot)
    # The implementation has a special-case for bytes label_name (added by 615bc01)
    nodes = load_freesurfer_label(str(annot), b'Cortex')
    assert nodes.size == 12


def _write_fake_gifti_labels(path, label_data, label_name_by_key):
    """Build a minimal gifti label file with the given per-vertex labels."""
    darr = GiftiDataArray(
        np.asarray(label_data, dtype=np.int32),
        intent='NIFTI_INTENT_LABEL',
        datatype='NIFTI_TYPE_INT32',
    )
    lt = GiftiLabelTable()
    for key, name in label_name_by_key.items():
        gl = GiftiLabel(key=key)
        gl.label = name
        lt.labels.append(gl)
    img = GiftiImage(labeltable=lt, darrays=[darr])
    img.to_filename(str(path))


def test_load_gifti_labels_returns_name_to_indices(tmp_path):
    path = tmp_path / 'lh.test.label.gii'
    data = np.array([0, 0, 1, 1, 2, 2, 2, 2], dtype=np.int32)
    _write_fake_gifti_labels(path, data, {0: 'medial', 1: 'A', 2: 'B'})

    result = load_gifti_labels(str(path))
    assert set(result.keys()) == {'medial', 'A', 'B'}
    np.testing.assert_array_equal(result['medial'], [0, 1])
    np.testing.assert_array_equal(result['A'], [2, 3])
    np.testing.assert_array_equal(result['B'], [4, 5, 6, 7])


def test_load_gifti_labels_handles_missing_label(tmp_path):
    """A label declared in the table but absent from data returns an empty array."""
    path = tmp_path / 'lh.empty.label.gii'
    data = np.array([0, 0, 1, 1], dtype=np.int32)
    _write_fake_gifti_labels(path, data, {0: 'medial', 1: 'A', 2: 'B'})
    result = load_gifti_labels(str(path))
    assert result['B'].size == 0


def _write_fake_dlabel_cifti(
    path,
    L_mask, L_data,
    R_mask, R_data,
    label_dict={0: ('???', (0, 0, 0, 0)),
                1: ('A', (1, 0, 0, 1)),
                2: ('B', (0, 1, 0, 1))},
):
    """Build a minimal dlabel.nii covering both hemispheres with given masks."""
    bm_l = cifti2.BrainModelAxis.from_mask(L_mask, name='CortexLeft')
    bm_r = cifti2.BrainModelAxis.from_mask(R_mask, name='CortexRight')
    bm = bm_l + bm_r
    la = cifti2.LabelAxis(['parcels'], [label_dict])
    data = np.concatenate([L_data, R_data])[None, :].astype(np.float32)
    hdr = cifti2.Cifti2Header.from_axes((la, bm))
    img = cifti2.Cifti2Image(data, hdr)
    img.to_filename(str(path))


def test_load_cifti_labels_left_hemi(tmp_path):
    """L hemi: 8 verts, cifti covers verts 1-6 (medial wall = {0,7})."""
    path = tmp_path / 'test.dlabel.nii'
    L_mask = np.array([0, 1, 1, 1, 1, 1, 1, 0])
    L_data = np.array([1, 1, 2, 2, 1, 2])
    R_mask = np.array([0, 1, 1, 1, 1, 1, 1, 0])
    R_data = np.array([2, 2, 1, 1, 2, 1])
    _write_fake_dlabel_cifti(path, L_mask, L_data, R_mask, R_data)

    result = load_cifti_labels(str(path), 'L')
    np.testing.assert_array_equal(np.sort(result['A']), [1, 2, 5])
    np.testing.assert_array_equal(np.sort(result['B']), [3, 4, 6])
    np.testing.assert_array_equal(np.sort(result['???']), [0, 7])


def test_load_cifti_labels_right_hemi(tmp_path):
    """R hemi reads the right-hemi grayordinates, not the left."""
    path = tmp_path / 'test.dlabel.nii'
    L_mask = np.array([0, 1, 1, 1, 1, 1, 1, 0])
    L_data = np.array([1, 1, 2, 2, 1, 2])
    R_mask = np.array([0, 1, 1, 1, 1, 1, 1, 0])
    R_data = np.array([2, 2, 1, 1, 2, 1])
    _write_fake_dlabel_cifti(path, L_mask, L_data, R_mask, R_data)

    result = load_cifti_labels(str(path), 'R')
    np.testing.assert_array_equal(np.sort(result['A']), [3, 4, 6])
    np.testing.assert_array_equal(np.sort(result['B']), [1, 2, 5])
    np.testing.assert_array_equal(np.sort(result['???']), [0, 7])


def test_load_cifti_labels_asymmetric_hemispheres(tmp_path):
    """Asymmetric mesh sizes: should not cross-contaminate (regression for the
    Lnverts/Rnverts offset bug in PR #20's draft)."""
    path = tmp_path / 'asym.dlabel.nii'
    L_mask = np.ones(6)
    L_data = np.array([1, 1, 2, 2, 1, 2])
    R_mask = np.ones(10)
    R_data = np.array([2, 2, 1, 1, 2, 1, 2, 1, 1, 2])
    _write_fake_dlabel_cifti(path, L_mask, L_data, R_mask, R_data)

    L = load_cifti_labels(str(path), 'L')
    R = load_cifti_labels(str(path), 'R')
    for v in L.values():
        if v.size:
            assert (v < 6).all()
    for v in R.values():
        if v.size:
            assert (v < 10).all()
    assert sum(arr.size for k, arr in L.items() if k != '???') == 6
    assert sum(arr.size for k, arr in R.items() if k != '???') == 10


def test_load_cifti_labels_invalid_hemi_raises(tmp_path):
    path = tmp_path / 'test.dlabel.nii'
    _write_fake_dlabel_cifti(
        path,
        np.ones(4), np.array([1, 1, 2, 2]),
        np.ones(4), np.array([1, 2, 1, 2]),
    )
    with pytest.raises(ValueError, match="hemi must be"):
        load_cifti_labels(str(path), 'X')


# ---- load_surface dispatcher --------------------------------------------------

def test_load_surface_passes_tuple_through():
    verts = np.zeros((5, 3), dtype=np.float64)
    tris = np.zeros((3, 3), dtype=np.int32)
    out = load_surface((verts, tris))
    assert out[0] is verts and out[1] is tris


def test_load_surface_from_freesurfer_geometry(tmp_path):
    verts = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    tris = np.array([[0, 1, 2]], dtype=np.int32)
    p = tmp_path / 'lh.test'
    nib.freesurfer.io.write_geometry(str(p), verts, tris)
    out_v, out_t = load_surface(str(p))
    np.testing.assert_allclose(out_v, verts)
    np.testing.assert_array_equal(out_t, tris)


def test_load_surface_from_gifti(tmp_path):
    verts = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]], dtype=np.float32)
    tris = np.array([[0, 1, 2]], dtype=np.int32)
    darr_v = GiftiDataArray(verts, intent='NIFTI_INTENT_POINTSET',
                             datatype='NIFTI_TYPE_FLOAT32')
    darr_t = GiftiDataArray(tris, intent='NIFTI_INTENT_TRIANGLE',
                             datatype='NIFTI_TYPE_INT32')
    img = GiftiImage(darrays=[darr_v, darr_t])
    p = tmp_path / 'lh.test.surf.gii'
    img.to_filename(str(p))
    out_v, out_t = load_surface(str(p))
    np.testing.assert_allclose(out_v, verts)
    np.testing.assert_array_equal(out_t, tris)


def test_load_surface_invalid_type_raises():
    with pytest.raises(TypeError):
        load_surface(42)


# ---- load_labels dispatcher ---------------------------------------------------

def test_load_labels_dispatches_freesurfer_annot(tmp_path):
    annot = tmp_path / 'lh.test.annot'
    _write_fake_annot(annot)
    nodes, medial = load_labels(str(annot), exceptions=['Medial_wall'])
    assert 'Medial_wall' not in nodes
    assert 'Cortex' in nodes
    assert nodes['Cortex'].size == 12
    np.testing.assert_array_equal(np.sort(medial), np.arange(4))


def test_load_labels_dispatches_gifti(tmp_path):
    p = tmp_path / 'lh.test.label.gii'
    data = np.array([0, 0, 1, 1, 2, 2], dtype=np.int32)
    _write_fake_gifti_labels(p, data, {0: '???', 1: 'A', 2: 'B'})
    nodes, medial = load_labels(str(p))
    assert set(nodes.keys()) == {'A', 'B'}  # '???' moved into medial
    np.testing.assert_array_equal(np.sort(medial), [0, 1])


def test_load_labels_dispatches_cifti(tmp_path):
    p = tmp_path / 'lh.test.dlabel.nii'
    L_mask = np.array([0, 1, 1, 1, 1, 1, 1, 0])
    L_data = np.array([1, 1, 2, 2, 1, 2])
    R_mask = np.array([0, 1, 1, 1, 1, 1, 1, 0])
    R_data = np.array([2, 2, 1, 1, 2, 1])
    _write_fake_dlabel_cifti(p, L_mask, L_data, R_mask, R_data)

    nodes, medial = load_labels(str(p), hemi='L')
    assert set(nodes.keys()) == {'A', 'B'}
    np.testing.assert_array_equal(np.sort(medial), [0, 7])


def test_load_labels_cifti_requires_hemi(tmp_path):
    p = tmp_path / 'lh.test.dlabel.nii'
    _write_fake_dlabel_cifti(
        p,
        np.ones(4), np.array([1, 1, 2, 2]),
        np.ones(4), np.array([1, 2, 1, 2]),
    )
    with pytest.raises(ValueError, match="hemi must be"):
        load_labels(str(p))


def test_load_labels_unrecognized_extension_raises(tmp_path):
    p = tmp_path / 'unknown.txt'
    p.write_text("not a label file")
    with pytest.raises(ValueError, match="unrecognized label file extension"):
        load_labels(str(p))


def test_load_labels_drops_extra_exceptions(tmp_path):
    p = tmp_path / 'lh.test.label.gii'
    data = np.array([0, 0, 1, 1, 2, 2, 3, 3], dtype=np.int32)
    _write_fake_gifti_labels(
        p, data, {0: '???', 1: 'A', 2: 'B', 3: 'C'},
    )
    nodes, _medial = load_labels(str(p), exceptions=['B'])
    assert set(nodes.keys()) == {'A', 'C'}
