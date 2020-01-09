
# TODO: See https://www.jeffknupp.com/blog/2013/08/16/open-sourcing-a-python-project-the-right-way/ for
#       how to create proper setup.py file

from setuptools import setup

setup(name='pymwts',
      version='0.1',
      description='Pyomo based multi-week tour scheduling model',
      author='Mark Isken',
      author_email='isken@oakland.edu',
      url='',
      packages=['pymwts', 'pymwts.pymwtsio'],
      entry_points = {
        'console_scripts': [
            'mwts = pymwts.mwts:main', 'mwts_phase2 = pymwts.solvemwts_phase2:main']},
      install_requires=['pyomo']

      )


"""
distutils approach

#from distutils.core import setup
scripts = ['scripts/mwts'],
"""
