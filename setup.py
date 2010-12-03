from setuptools import setup, find_packages
import os

from bioscripts.blasaid import __version__ as VERSION

setup (
  name='bioscripts-blasaid',
  version=VERSION,
  description="Automatic blasting and summarizing of NGS data",
  long_description=open("README.txt").read() + "\n" +
	 open(os.path.join("docs", "HISTORY.txt")).read(),
  # Get more strings from
  # http://pypi.python.org/pypi?%3Aaction=list_classifiers
  classifiers=[
	 'Programming Language :: Python',
	 'Development Status :: 4 - Beta',
	 'Environment :: Console',
	 'Intended Audience :: Science/Research',
	 'Topic :: Scientific/Engineering :: Bio-Informatics',
  ],
  keywords='bioinformatics sequencing blast',
  author='Paul Agapow',
  author_email='paul-michael.agapow@hpa.org.uk',
  url='http://www.agapow.net/software/bioscripts.blasaid',
  license='GPL',
  packages=find_packages(exclude=['ez_setup']),
  namespace_packages=['bioscripts'],
  include_package_data=True,
  zip_safe=False,
  install_requires=[
	 'setuptools',
	 'biopython >= 1.49',
	 'PyYaml',
  ],
  entry_points={
	 'console_scripts': [
		'convbioseq = bioscripts.convert.convbioseq:main',
		'convalign = bioscripts.convert.convalign:main',
	 ],
  },

)
