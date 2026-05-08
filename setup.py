import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


def load_requirements(fname):
    with open(fname) as fh:
        return [
            line.strip()
            for line in fh
            if line.strip() and not line.lstrip().startswith("#")
        ]


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
  python_requires='>=3.8',
  install_requires=load_requirements("requirements.txt")
)
