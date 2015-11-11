from nipype.interfaces.io import FreeSurferSource, DataSink
from nipype.interfaces.utility import IdentityInterface
from nipype import Workflow, Node, MapNode, JoinNode, Function
import nibabel as nib
import numpy as np
import os
import surfdist as sd
import csv

atlas = 'aparc.a2009s.annot'
labs = 'cortex.label'
subjects_dir = '/home/mfalkiewicz/experiments/SDT1/Surfaces'
subject_list = ['AD0303','AH2805']
sources = ['S_central']
target = ['fsaverage4','fsaverage5']
hemi = ['lh','rh']


def trimming(itemz, phrase):
  item = [x for x in itemz if phrase in x][0]
  return item


def calc_surfdist(surface, labels, annot, reg, origin, target):
  import nibabel as nib
  import numpy as np
  import os
  from surfdist import load, utils, surfdist
  import csv

  """ inputs:
  surface - surface file (e.g. lh.pial, with full path)
  labels - label file (e.g. lh.cortex.label, with full path)
  annot - annot file (e.g. lh.aparc.a2009s.annot, with full path)
  reg - registration file (lh.sphere.reg)
  origin - the label from which we calculate distances
  target - target surface (e.g. fsaverage4)
  """

  # Load stuff
  surf = nib.freesurfer.read_geometry(surface)
  cort = np.sort(nib.freesurfer.read_label(labels))
  src  = load.load_freesurfer_label(annot, origin, cort)

  # Calculate distances
  dist = surfdist.dist_calc(surf, cort, src)

  # Project distances to target
  trg = nib.freesurfer.read_geometry(target)[0]
  native = nib.freesurfer.read_geometry(reg)[0]
  idx_trg_to_native = utils.find_node_match(trg, native)[0]

  # Get indices in trg space 
  distt = dist[idx_trg_to_native]
  
  # Write to file and return file handle
  filename = os.path.join(os.getcwd(),'distances.csv')
  distt.tofile(filename,sep=",")

  return filename

def stack_files(files,source,target):
  """
  This function takes a list of files as input and vstacks them
  """
  import csv
  import os
  import numpy as np

  fname = "sdist_%s_%s.csv" % (source, target)
  filename = os.path.join(os.getcwd(),fname)

  alldist = []

  for dfile in files:
    alldist.append(np.genfromtxt(dfile, delimiter=','))

  alldist = np.array(alldist)
  alldist.tofile(filename,",")

  return filename


def create_surfdist_workflow(subjects_dir,
                             subject_list,
                             sources,
                             target,
                             hemi,
                             atlas,
                             labs):

  sd = Workflow(name='SurfDist')
    
  # Run a separate tree for each template, hemisphere and source structure
  infosource = Node(IdentityInterface(fields=['template','hemi','source']), name="infosource")
  infosource.iterables = [('template', target),('hemi', hemi),('source',sources)]

  # Get template files
  fsst = Node(FreeSurferSource(),name='FS_Source_template')
  fsst.inputs.subjects_dir = subjects_dir

  sd.connect(infosource,'template',fsst,'subject_id')
  sd.connect(infosource,'hemi',fsst,'hemi')

  # Get subjects
  fss = MapNode(FreeSurferSource(),iterfield='subject_id',name='FS_Source')
  fss.inputs.subject_id = subject_list
  fss.inputs.subjects_dir = subjects_dir

  sd.connect(infosource,'hemi',fss,'hemi')

  # Trim labels
  tlab = MapNode(Function(input_names=['itemz','phrase'],
                        output_names=['item'], function=trimming),
                        iterfield='itemz', name='tlab')
  tlab.inputs.phrase = labs
  sd.connect(fss,'label',tlab,'itemz')

  # Trim annotations
  tannot = MapNode(Function(input_names=['itemz','phrase'],
                        output_names=['item'], function=trimming),
                        iterfield = ['itemz'], name='tannot')
  tannot.inputs.phrase = atlas
  sd.connect(fss,'annot',tannot,'itemz')

  # Calculate distances for each hemi
  sdist = MapNode(Function(input_names=['surface','labels','annot','reg','origin','target'],
                        output_names=['distances'], function=calc_surfdist), 
                        iterfield = ['surface','labels','annot','reg'], name='sdist')
  sd.connect(infosource,'source',sdist,'origin')
  sd.connect(fss,'pial',sdist,'surface')
  sd.connect(tlab,'item',sdist,'labels')
  sd.connect(tannot,'item',sdist,'annot')
  sd.connect(fss,'sphere_reg',sdist,'reg')
  sd.connect(fsst,'sphere_reg',sdist,'target')
  
  # Gather data for each hemi from all subjects
  bucket = Node(Function(input_names=['files','source','target'],output_names=['group_dist'], 
                         function=stack_files), name='bucket')
  sd.connect(infosource,'source',bucket,'source')
  sd.connect(infosource,'template',bucket,'target')
  sd.connect(sdist,'distances',bucket,'files')

  # Sink the data
  datasink = Node(DataSink(), name='sinker')
  datasink.inputs.base_directory = os.getcwd()
  sd.connect(infosource,'hemi',datasink,'container')
  sd.connect(bucket,'group_dist',datasink,'group_distances')

  return sd

sdist = create_surfdist_workflow(subjects_dir, subject_list, sources, target, hemi, atlas, labs)
sdist.base_dir = os.getcwd()
#sdist.run()
sdist.run(plugin='CondorDAGMan')
