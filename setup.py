from setuptools import setup
from distutils.core import setup

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
)
