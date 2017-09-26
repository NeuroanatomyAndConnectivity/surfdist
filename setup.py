from distutils.core import setup
setup(
  name = 'surfdist',
  packages = ['surfdist'],
  install_requires = ['numpy',
                      'gdist',
                      'nibabel',
                      'scipy'],
  version = '0.13.14',
  description = 'For calculating exact geodesic distances on cortical surface meshes',
  author = 'Daniel Margulies',
  author_email = 'margulies@cbs.mpg.de',
  url = 'https://github.com/NeuroanatomyAndConnectivity/surfdist',
  download_url = 'https://github.com/NeuroanatomyAndConnectivity/surfdist/tarball/0.13.14',
  keywords = ['geodesic', 'distance', 'brain', 'cortex'],
  classifiers = [],
)
