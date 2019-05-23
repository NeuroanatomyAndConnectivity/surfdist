from setuptools import setup

try:
    # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:
    # for pip <= 9.0.3
    from pip.req import parse_requirements

def load_requirements(fname):
    reqs = parse_requirements(fname, session="test")
    return [str(ir.req) for ir in reqs]

setup(
  name = 'surfdist',
  packages = ['surfdist'],
  version = '0.14.0',
  description = 'For calculating exact geodesic distances on cortical surface meshes',
  author = 'Daniel Margulies',
  author_email = 'daniel.margulies@gmail.com',
  url = 'https://github.com/NeuroanatomyAndConnectivity/surfdist',
  keywords = ['geodesic', 'distance', 'brain', 'cortex'],
  license='LICENSE.txt',
  long_description=open('README.md').read(),
  install_requires=load_requirements("requirements.txt")
)
