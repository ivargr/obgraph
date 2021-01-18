from setuptools import setup
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("obgraph.cython_traversing",
              ["obgraph/cython_traversing.pyx"],
              libraries=["m"],
              extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
              extra_link_args=['-fopenmp']
              )]

setup(name='obgraph',
      version='0.0.1',
      description='obgraph',
      url='http://github.com/ivargr/obgraph',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=["obgraph"],
      zip_safe=False,
      install_requires=['numpy', 'tqdm'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points={
            'console_scripts': ['obgraph=obgraph.command_line_interface:main']
      },
      cmdclass = {"build_ext": build_ext},
      ext_modules = ext_modules
)
