from setuptools import setup
from setuptools.command.install import install
import subprocess as sub
import os

class MyInstall(install):
    def run(self):
        install.run(self)
        sub.call(['julia', 'InstallPkg.jl'])
        

setup(name='pyLLE',
      version='0.1',
      description='LLE Solver',
      url='',
      author='Greg Moille',
      author_email='gregory.moille@nist.gov',
      license='MIT',
      packages=['pyLLE'],
      install_requires=[
          'scipy',
          'numpy',
          'matplotlib',
          'h5py',
          'pickle',
          'logging',
          'time',
          'sys',
          'subprocess',
          'os',
          'inspect',
          'copy',
          'shutil',
          'tempfile',
          'prettytable',
          'warnings',
          'matplotlib',
          'logging',
      ],
      package_data={'': ['*.jl']},
      include_package_data=True,
      zip_safe=False,
      cmdclass={'install': MyInstall},)