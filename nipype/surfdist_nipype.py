#!/usr/bin/env python

from nipype.interfaces.io import FreeSurferSource, DataSink
from nipype.interfaces.utility import IdentityInterface
from nipype import Workflow, Node, MapNode, JoinNode, Function
import nibabel as nib
import numpy as np
import os
import surfdist as sd
import csv

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

def stack_files(files, hemi, source, target):
  """
  This function takes a list of files as input and vstacks them
  """
  import csv
  import os
  import numpy as np

  fname = "sdist_%s_%s_%s.csv" % (hemi, source, target)
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
                             labs,
                             name):

  sd = Workflow(name=name)
    
  # Run a separate tree for each template, hemisphere and source structure
  infosource = Node(IdentityInterface(fields=['template','hemi','source']), name="infosource")
  infosource.iterables = [('template', target),('hemi', hemi),('source',sources)]

  # Get template files
  fsst = Node(FreeSurferSource(),name='FS_Source_template')
  fsst.inputs.subjects_dir = subjects_dir

  sd.connect(infosource,'template',fsst,'subject_id')
  sd.connect(infosource,'hemi',fsst,'hemi')

  # Get subjects
  fss = Node(FreeSurferSource(),iterfield='subject_id',name='FS_Source')
  fss.iterables = ('subject_id', subject_list)
  fss.inputs.subjects_dir = subjects_dir

  sd.connect(infosource,'hemi',fss,'hemi')

  # Trim labels
  tlab = Node(Function(input_names=['itemz','phrase'],
                        output_names=['item'], function=trimming),
                        name='tlab')
  tlab.inputs.phrase = labs
  sd.connect(fss,'label',tlab,'itemz')

  # Trim annotations
  tannot = Node(Function(input_names=['itemz','phrase'],
                        output_names=['item'], function=trimming),
                        name='tannot')
  tannot.inputs.phrase = atlas
  sd.connect(fss,'annot',tannot,'itemz')

  # Calculate distances for each hemi
  sdist = Node(Function(input_names=['surface','labels','annot','reg','origin','target'],
                        output_names=['distances'], function=calc_surfdist), 
                        name='sdist')
  sd.connect(infosource,'source',sdist,'origin')
  sd.connect(fss,'pial',sdist,'surface')
  sd.connect(tlab,'item',sdist,'labels')
  sd.connect(tannot,'item',sdist,'annot')
  sd.connect(fss,'sphere_reg',sdist,'reg')
  sd.connect(fsst,'sphere_reg',sdist,'target')
  
  # Gather data for each hemi from all subjects
  bucket = JoinNode(Function(input_names=['files','hemi','source','target'],output_names=['group_dist'], 
                         function=stack_files), joinsource = 'fss', joinfield = 'files', name='bucket')
  sd.connect(infosource,'source',bucket,'source')
  sd.connect(infosource,'template',bucket,'target')
  sd.connect(infosource,'hemi',bucket,'hemi')
  sd.connect(sdist,'distances',bucket,'files')

  # Sink the data
  datasink = Node(DataSink(parametrization=False), name='sinker')
  datasink.inputs.base_directory = os.getcwd()
  sd.connect(infosource,'hemi',datasink,'container')
  sd.connect(bucket,'group_dist',datasink,'group_distances')

  return sd



def create_workflow(args, name=None):

    with open(args.subject_file) as f:
      subject_list = f.read().splitlines()
    
    if name is None:
        name = 'surfdist'

    kwargs = dict(subjects_dir = args.subjects_dir,
                  subject_list = subject_list,
                  sources = args.sources,
                  target = args.target_surfs,
                  hemi = args.hemi,
                  atlas = args.annot,
                  labs = args.labels,
                  name=name)
    wf = create_surfdist_workflow(**kwargs)

    return wf


if __name__ == "__main__":
    from argparse import ArgumentParser, RawTextHelpFormatter
    import os
    defstr = ' (default %(default)s)'
    parser = ArgumentParser(description='''This script generates and runs a nipype pipeline for calculating distances from source label(s) 
on a Freesurfer surface. After calculating the distances in native space it transforms 
the distances into selected target space and creates a CSV file containing data for all 
subjects. This table can be used for permutation testing in PALM.''',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-s", "--subject_ids", dest="subject_file",
                        help="Subject list file", required=True)
    parser.add_argument("-sd", "--subjects_dir", dest="subjects_dir",
                        help="FreeSurfer subject directory", required=True)
    parser.add_argument("-t", "--target_surfaces", dest="target_surfs", nargs="+",
                        default=['fsaverage5'],
                        help="FreeSurfer target surfaces" + defstr)
    parser.add_argument("-a", "--annot", dest="annot",
                        default='aparc.a2009s', 
                        help="Annotation for source label(s)" + defstr)
    parser.add_argument("-l", "--label", dest="labels",
                        default='cortex',
                        help="Label(s)" + defstr)
    parser.add_argument("-src", "--source", dest="sources", nargs = "+",
                        default=['S_central'], 
                        help="Label(s) to calculate distances from" + defstr)
    parser.add_argument("-hemi", "--hemi", dest="hemi", nargs = "+",
                        default=['lh','rh'],
                        help="Hemisphere(s) for distance calculation" + defstr)
    parser.add_argument("-o", "--output_dir", dest="sink",
                        help="Output directory base")
    parser.add_argument("-w", "--work_dir", dest="work_dir",
                        help="Output directory base")
    parser.add_argument("-p", "--plugin", dest="plugin",
                        default='Linear',
                        help="Plugin to use")
    parser.add_argument("--plugin_args", dest="plugin_args",
                        help="Plugin arguments")
    args = parser.parse_args()


wf = create_workflow(args)
if args.work_dir:
  work_dir = os.path.abspath(args.work_dir)
else:
  work_dir = os.getcwd()

wf.base_dir = work_dir

if args.plugin_args:
  wf.run(args.plugin, plugin_args=eval(args.plugin_args))
else:
  wf.run(args.plugin)

#wf.write_graph(dotfilename='func_preproc.dot', graph2use='exec', format='pdf', simple_form=False)

