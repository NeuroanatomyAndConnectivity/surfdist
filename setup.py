from distutils.core import setup
setup(
  name = 'surfdist',
  packages = ['surfdist'],
  install_requires = ['numpy',
                      'gdist',
                      'nibabel',
                      'scipy',
                      'numba'],
  version = '0.14.0',
  description = 'For calculating exact geodesic distances on cortical surface meshes',
  author = 'Daniel Margulies',
  author_email = 'daniel.margulies@gmail.com',
  url = 'https://github.com/NeuroanatomyAndConnectivity/surfdist',
  download_url = 'https://github.com/NeuroanatomyAndConnectivity/surfdist/tarball/0.14.0',
  keywords = ['geodesic', 'distance', 'brain', 'cortex'],
  classifiers = [],
)
