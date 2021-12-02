from setuptools import setup
from distutils.core import setup

setup(name='obgraph',
      version='0.0.5',
      description='obgraph',
      url='http://github.com/ivargr/obgraph',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=["obgraph"],
      zip_safe=False,
      install_requires=['numpy==1.20.3', 'tqdm', 'pathos', 'graph_kmer_index'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points={
            'console_scripts': ['obgraph=obgraph.command_line_interface:main']
      },
)

"""

rm -rf dist
python3 setup.py sdist
twine upload --skip-existing dist/*

"""