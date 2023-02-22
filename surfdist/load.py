import nibabel as nib
import numpy as np
import sys
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

def load_FS_annot(annot_input,exceptions=[]):
    label_list = get_freesurfer_label(annot_input,verbose=False)
    tmp=[i.decode('utf-8') for i in label_list]
    label_list=tmp
    del tmp
    if len(exceptions)>0:
        label_list=list(set(label_list)-set(exceptions))
    labelNodes=[load_freesurfer_label(annot_input, roi) for roi in label_list]
    labelNodes=dict(zip(label_list,labelNodes))
    return labelNodes

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

def intSettoList(data):
    #### helper function for loading ciftis### 
    data=list(set(data))
    return [int(x) for x in data]

def load_cifti_labels(ciftiLabel,hemi):
    """Load in the values from a cifti label. Input is expected to contain both left and right hemispheres. 
    Current iteration has been tested on Yeo and Schaefer cifti files in fs_lr 32k space
     Working on HCP flag which takes HCP generated cifti files. until then use HCP giftis """
    cifti = nib.load( ciftiLabel)
    cifti_data = cifti.get_fdata(dtype=np.float32).squeeze()
    cifti_hdr = cifti.header
    nifti_hdr = cifti.nifti_header
    axes = [cifti_hdr.get_axis(i) for i in range(cifti.ndim)]
    Lnverts=axes[1].nvertices['CIFTI_STRUCTURE_CORTEX_LEFT']
    Rnverts=axes[1].nvertices['CIFTI_STRUCTURE_CORTEX_RIGHT']
    
    vertexValues=cifti.get_fdata().squeeze()
    if len(vertexValues)==(Lnverts+Rnverts):
        print('this is NOT a HCP cifti. Labels are split by dimensions of total surface')
        LvertVals=vertexValues[0:Lnverts]
        RvertVals=vertexValues[Rnverts:]

        rois={}
        for i in axes[0].label[0]:
            rois[i]=np.where(cifti_data==i)[0]
        MedialWall=rois.pop(0)
        LMedialWall=MedialWall[MedialWall<Lnverts]
        RMedialWall=MedialWall[MedialWall>Lnverts]-Rnverts
        split=int(len(rois)/2)+1
        
        Lrois=range(1,split)
        Rrois=range(split,len(rois)+1)

        LlabelNodes={}
        for i in Lrois:
            LlabelNodes[axes[0].label[0][i][0]]=rois[i]

        
        RlabelNodes={}
        for i in Rrois:
            RlabelNodes[axes[0].label[0][i][0]]=rois[i]-Rnverts
        if hemi=='L':
            LlabelNodes['???']=LMedialWall
            return LlabelNodes
        elif hemi=='R':
            RlabelNodes['???']=RMedialWall
            return RlabelNodes
    else:
        print('labels are not split by hemisphere. function will use')
        print('hcp indices') #### these indices are specified for cortex in the HCP output ciftis
        lcort=slice(0,29696)
        rcort=slice(29696, 59412)
        LvertVals=vertexValues[lcort]
        RvertVals=vertexValues[rcort]
        #### create a dictionary of labels to vertices 
        if hemi=='L':
            LvertIdx=intSettoList(LvertVals)
            LlabelNodes={}
            for value in LvertIdx:
                LlabelNodes[axes[0].label[0][value][0]]=np.where(LvertVals==value)[0]
            return LlabelNodes
        elif hemi=='R':
            RvertIdx=intSettoList(RvertVals)
            RlabelNodes={}
            for value in RvertIdx:
                RlabelNodes[axes[0].label[0][value][0]]=np.where(RvertVals==value)[0]
        return RlabelNodes