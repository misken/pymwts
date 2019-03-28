

from setuptools import setup

setup(name='pymwts',
      version='1.0',
      description='Pyomo based mwts',
      author='Mark Isken',
      author_email='isken@oakland.edu',
      url='http://www.hselab.org/machinery',
      packages=['pymwts', 'pymwts.pymwtsio'],
      entry_points = {
        'console_scripts': [
            'mwts = pymwts.mwts:main', 'mwts_phase2 = pymwts.solvemwts_phase2:main']}
      
     )


"""
distutils approach

#from distutils.core import setup
scripts = ['scripts/mwts'],
"""
