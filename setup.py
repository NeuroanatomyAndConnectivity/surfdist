import setuptools 

with open("README.md", "r") as fh:
    long_description = fh.read()
    
try:
    # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:
    # for pip <= 9.0.3
    from pip.req import parse_requirements

def load_requirements(fname):
    reqs = parse_requirements(fname, session="test")
    return [str(ir.req) for ir in reqs]

setuptools.setup(
  name = 'surfdist',
  packages = ['surfdist'],
  version = '0.15.4',
  description = 'For calculating exact geodesic distances on cortical surface meshes',
  long_description=long_description,
  long_description_content_type="text/markdown",
  author = 'Daniel Margulies',
  author_email = 'daniel.margulies@gmail.com',
  url = 'https://github.com/NeuroanatomyAndConnectivity/surfdist',
  keywords = ['geodesic', 'distance', 'brain', 'cortex'],
  license='LICENSE.txt',
  python_requires='>=3.6',  
  install_requires=load_requirements("requirements.txt")
)
