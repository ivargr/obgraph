from setuptools import setup
from distutils.core import setup
from Cython.Build import cythonize

setup(name='obgraph',
      version='0.0.9',
      description='obgraph',
      url='http://github.com/ivargr/obgraph',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=["obgraph"],
      zip_safe=False,
      install_requires=['numpy==1.20.3', 'tqdm', 'pathos', 'graph_kmer_index>=0.0.15', 'shared_memory_wrapper>=0.0.7'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points={
            'console_scripts': ['obgraph=obgraph.command_line_interface:main']
      },
      ext_modules = cythonize("obgraph/cython_traversing.pyx"),
)

"""

rm -rf dist
python3 setup.py sdist
twine upload --skip-existing dist/*

"""