import nibabel as nib
import numpy as np


def load_freesurfer_label(annot_input, label_name, cortex=None):
    """
    Get source node list for a specified freesurfer label.

    Inputs
    -------
    annot_input : freesurfer annotation label file
    label_name : freesurfer label name
    cortex : not used
    """

    if cortex is not None:
        print("Warning: cortex is not used to load the freesurfer label")

    labels, color_table, names = nib.freesurfer.read_annot(annot_input)
    names = [i.decode('utf-8') for i in names]
    if isinstance(label_name, bytes):
        label_name = label_name.decode('utf-8')
    label_value = names.index(label_name)
    label_nodes = np.array(np.where(np.isin(labels, label_value)), dtype=np.int32)

    return label_nodes


def get_freesurfer_label(annot_input, verbose = True):
    """
    Print freesurfer label names.
    """
    labels, color_table, names = nib.freesurfer.read_annot(annot_input)
    if verbose:
        print(names)
    return names


def load_gifti_labels(gifti_label):
    """
    Get a mapping of label name -> vertex indices from a GIFTI label file.

    Inputs
    -------
    gifti_label : str or path-like
        Path to a hemisphere-specific GIFTI label file (e.g. ``*.label.gii``).

    Returns
    -------
    label_nodes : dict
        Dictionary mapping each label name to a 1-D numpy array of vertex
        indices belonging to that label.
    """
    gii = nib.load(gifti_label)
    data = gii.darrays[0].data
    key_to_name = gii.labeltable.get_labels_as_dict()
    return {name: np.where(data == key)[0] for key, name in key_to_name.items()}
