import setuptools 

with open("README.md", "r") as fh:
    long_description = fh.read()

from pip._internal.req import parse_requirements
def load_requirements(fname):
    install_reqs = parse_requirements(fname, session=False)
    try:
        requirements = [str(ir.req) for ir in install_reqs]
    except:
        requirements = [str(ir.requirement) for ir in install_reqs]    
    return requirements 

setuptools.setup(
  name = 'surfdist',
  packages = ['surfdist'],
  version = '0.15.5',
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
