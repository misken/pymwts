"""Setup file for chime scenario runner
"""
__version__ = "0.2.0"
__author__ = "misken"

from setuptools import setup

setup(name='pymwts',
      version=__version__,
      description='Pyomo based multi-week tour scheduling model',
      author=__author__,
      author_email='isken@oakland.edu',
      url='',
      packages=['pymwts', 'pymwts.pymwtsio'],
      entry_points = {
        'console_scripts': [
            'pymwts = pymwts.mwts:main']},
      python_requires='>=3.7',
      install_requires=['pyomo',
                        'pandas']

      )

