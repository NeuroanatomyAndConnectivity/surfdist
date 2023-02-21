import nibabel as nib
import numpy as np
import sys
from surfdist.utils import intSettoList

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
    label_value = names.index(label_name)
    label_nodes = np.array(np.where(np.in1d(labels, label_value)), dtype=np.int32)

    return label_nodes

def get_freesurfer_label(annot_input, verbose = True):
    """
    Print freesurfer label names.
    """
    labels, color_table, names = nib.freesurfer.read_annot(annot_input)
    if verbose:
        print(names)
    return names

def load_gifti_labels(giftiLabel):
    """Load in the values from a gifti label. Input is heimsphere specific.
    This means you load in a left or right gifti file separately"""
    gii = nib.load(giftiLabel)
    gii_data = gii.darrays[0].data 
    
    labels=gii.labeltable.get_labels_as_dict()
    label_vals={y: x for x, y in labels.items()}
    labels=label_vals
    del label_vals
    labelNodes={}
    for k in labels:
        labelNodes[k]=np.where(labels[k]==gii_data)[0]
    return labelNodes

def load_cifti_labels(ciftiLabel):
    """Load in the values from a cifti label. Input is expected to contain both left and right hemispheres. 
    Current iteration has been tested on Yeo and Schaefer cifti files in fs_lr 32k space
     Working on HCP flag which takes HCP generated cifti files. until then use HCP giftis """
    cifti = nib.load( ciftiLabel)
    cifti_data = cifti.get_fdata(dtype=np.float32)
    cifti_hdr = cifti.header
    nifti_hdr = cifti.nifti_header
    axes = [cifti_hdr.get_axis(i) for i in range(cifti.ndim)]
    Lnverts=axes[1].nvertices['CIFTI_STRUCTURE_CORTEX_LEFT']
    Rnverts=axes[1].nvertices['CIFTI_STRUCTURE_CORTEX_RIGHT']
    
    vertexValues=cifti.get_fdata().squeeze()
    if len(vertexValues)==(Lnverts+Rnverts):
        LvertVals=vertexValues[0:Lnverts]
        RvertVals=vertexValues[Rnverts:]
    else:
        print('labels are not split by hemisphere. function will use')
        print('hcp indices') #### these indices are specified for cortex in the HCP output ciftis
        lcort=slice(0,29696)
        rcort=slice(29696, 59412)
        LvertVals=vertexValues[lcort]
        RvertVals=vertexValues[rcort]

    #### create a dictionary of labels to vertices 
    LvertIdx=intSettoList(LvertVals)
    LlabelNodes={}
    for value in LvertIdx:
        LlabelNodes[axes[0].label[0][value][0]]=np.where(LvertVals==value)[0]
    
    #### do the right now
    RvertIdx=intSettoList(RvertVals)
    RlabelNodes={}
    for value in RvertIdx:
        RlabelNodes[axes[0].label[0][value][0]]=np.where(RvertVals==value)[0]
    return LlabelNodes,RlabelNodes