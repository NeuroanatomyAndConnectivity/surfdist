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


_CIFTI_HEMI_STRUCTURE = {
    'L': 'CIFTI_STRUCTURE_CORTEX_LEFT',
    'R': 'CIFTI_STRUCTURE_CORTEX_RIGHT',
}


def load_cifti_labels(cifti_label, hemi, medial_wall_key='???'):
    """
    Get a mapping of label name -> mesh vertex indices for one hemisphere
    of a CIFTI dlabel file.

    Returned indices are in the original surface mesh's vertex space, so
    they can be passed directly to ``analysis.dist_calc`` together with a
    surface loaded from the matching ``*.surf.gii`` file.

    Vertices that the CIFTI does not cover for the requested hemisphere
    (typically the medial wall, which HCP-style CIFTIs exclude from the
    grayordinate set) are returned under ``medial_wall_key``.

    Inputs
    -------
    cifti_label : str or path-like
        Path to a CIFTI dlabel file (e.g. ``*.dlabel.nii``).
    hemi : {'L', 'R'}
        Which hemisphere's labels to extract.
    medial_wall_key : str
        Dict key under which to return the mesh vertices that the CIFTI
        does not cover for ``hemi``. Defaults to ``'???'`` to match the
        HCP convention.

    Returns
    -------
    label_nodes : dict
        Dictionary mapping each label name to a 1-D numpy array of vertex
        indices in the original surface mesh.
    """
    if hemi not in _CIFTI_HEMI_STRUCTURE:
        raise ValueError(
            f"hemi must be 'L' or 'R', got {hemi!r}"
        )
    structure = _CIFTI_HEMI_STRUCTURE[hemi]

    img = nib.load(cifti_label)
    label_axis = img.header.get_axis(0)
    bm_axis = img.header.get_axis(1)

    # Locate the requested hemisphere's grayordinate slice and mesh mapping
    for name, slc, struct in bm_axis.iter_structures():
        if name == structure:
            data = img.get_fdata().squeeze()[slc].astype(np.int64)
            mesh_vertex_indices = np.asarray(struct.vertex)
            n_mesh = struct.nvertices[name]
            break
    else:
        raise ValueError(
            f"hemisphere {hemi!r} ({structure}) not found in {cifti_label}"
        )

    label_dict = label_axis.label[0]  # {key: (name, rgba)}
    out = {}
    for key in np.unique(data):
        key_int = int(key)
        if key_int in label_dict:
            label_name = label_dict[key_int][0]
        else:
            label_name = str(key_int)
        out[label_name] = mesh_vertex_indices[data == key]

    # Mesh vertices the CIFTI does not represent for this hemisphere
    medial_mask = np.ones(n_mesh, dtype=bool)
    medial_mask[mesh_vertex_indices] = False
    medial_indices = np.where(medial_mask)[0]
    if medial_indices.size > 0:
        out[medial_wall_key] = medial_indices
    elif medial_wall_key not in out:
        out[medial_wall_key] = medial_indices  # empty array

    return out


def load_surface(surface):
    """
    Load a cortical surface from a path, or pass a pre-built tuple through.

    Inputs
    -------
    surface : tuple, str, or path-like
        Either a ``(vertices, triangles)`` tuple already in surfdist's
        expected format, a freesurfer geometry filename (e.g. ``lh.pial``),
        or a path to a GIFTI surface file (``*.surf.gii``).

    Returns
    -------
    surf : (vertices, triangles) tuple
    """
    if isinstance(surface, tuple) and len(surface) == 2:
        return surface
    if not (isinstance(surface, (str, bytes)) or hasattr(surface, '__fspath__')):
        raise TypeError(
            f"surface must be a (vertices, triangles) tuple or a path, "
            f"got {type(surface).__name__}"
        )
    path = str(surface)
    if path.endswith('.surf.gii'):
        gii = nib.load(path)
        return (gii.darrays[0].data, gii.darrays[1].data)
    return nib.freesurfer.read_geometry(path)


def _load_freesurfer_annot_dict(annot_path):
    """Read a freesurfer annot into ``{label_name: vertex_indices}``."""
    label_data, _ctab, names = nib.freesurfer.read_annot(annot_path)
    decoded = [n.decode('utf-8') if isinstance(n, (bytes, bytearray)) else n
               for n in names]
    return {
        name: np.where(label_data == idx)[0]
        for idx, name in enumerate(decoded)
    }


def load_labels(labels, hemi=None, exceptions=()):
    """
    Load parcel labels from a freesurfer ``.annot``, GIFTI ``.label.gii``,
    or CIFTI ``.dlabel.nii`` file, dispatching by extension.

    Inputs
    -------
    labels : str or path-like
        Path to the label file.
    hemi : {'L', 'R'} or None
        Required for CIFTI dlabel files; ignored otherwise.
    exceptions : iterable of str
        Label names to drop from the returned dict. For freesurfer annots
        only, the first entry doubles as the medial-wall label (since
        annots have no built-in '???' key).

    Returns
    -------
    label_nodes : dict
        ``{label_name: vertex_indices}`` for non-medial-wall labels.
    medial_wall : ndarray
        Vertex indices of the medial wall (empty if the file does not
        encode one).
    """
    path = str(labels)
    exceptions = list(exceptions)

    if path.endswith('.label.gii'):
        nodes = load_gifti_labels(path)
    elif path.endswith('.dlabel.nii'):
        if hemi is None:
            raise ValueError(
                "hemi must be 'L' or 'R' for a CIFTI dlabel file"
            )
        nodes = load_cifti_labels(path, hemi)
    elif path.endswith('.annot'):
        nodes = _load_freesurfer_annot_dict(path)
    else:
        raise ValueError(
            f"unrecognized label file extension: {path!r} "
            "(expected .label.gii, .dlabel.nii, or .annot)"
        )

    medial_key = '???'
    if medial_key in nodes:
        medial_wall = nodes.pop(medial_key)
    elif exceptions and exceptions[0] in nodes:
        medial_wall = nodes.pop(exceptions[0])
        exceptions = exceptions[1:]
    else:
        medial_wall = np.array([], dtype=np.int64)

    for ex in exceptions:
        nodes.pop(ex, None)

    return nodes, medial_wall
