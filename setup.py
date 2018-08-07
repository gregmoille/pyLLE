from setuptools import setup
from setuptools.command.install import install
from setuptools.command.build_py import build_py
import subprocess as sub
import os

class MyInstall(install):
    def run(self):
        install.run(self)
        print('-'*10)
        print("TEST")
        print('-'*10)
        print('-'*10)
        print('-'*10)
        print('-'*10)
        print('-'*10)
        # sub.call(['julia', 'InstallPkg.jl'])
        

class MyBuild(build_py):
    def run(self):
        build_py.run(self)
        for ii in range(20):
            print('-'*10)
        # sub.call(['julia', 'InstallPkg.jl'])
        
        

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(name='pyLLE',
      version='0.0.11',
      description='LLE Solver',
      url='https://github.com/gregmoille/pyLLE',
      author='Greg Moille',
      author_email='gregory.moille@nist.gov',
      license='MIT',
      long_description=long_description,
      packages=['pyLLE'],
      install_requires=[
          'scipy',
          'numpy',
          'matplotlib',
          'h5py',
          'prettytable',
          'matplotlib',
      ],
      package_data={'': ['*.jl']},
      include_package_data=True,
      zip_safe=False,
      # cmdclass={'install': MyInstall},
      cmdclass={'build_py': MyBuild},
      classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),)


