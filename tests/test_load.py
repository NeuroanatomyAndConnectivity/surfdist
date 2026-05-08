import nibabel as nib
import numpy as np

from surfdist.load import get_freesurfer_label, load_freesurfer_label


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
