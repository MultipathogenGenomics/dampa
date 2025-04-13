from setuptools import setup,find_packages
from dampa import __version__

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='dampa',
      version=__version__,
      description='Diversity Aware Metagenomic Panel Assignment',
      long_description=readme(),
      classifiers=[
          'License :: OSI Approved :: GPLv3',
          'Programming Language :: Python :: 3.10',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research',
      ],
      keywords='targetted metagenomics, probe design, metagenomics, bioinformatics',
      url='https://github.com/MultipathogenGenomics/dampa',
      author='Michael Payne',
      author_email='michael.payne@sydney.edu.au',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['dampa=dampa.dampa:main'],
      },
      zip_safe=False)