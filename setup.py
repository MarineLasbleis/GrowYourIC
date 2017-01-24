from __future__ import absolute_import, print_function
import os
import re
import sys


here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    # intentionally *not* adding an encoding option to open
    return codecs.open(os.path.join(here, *parts), 'r').read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
raise RuntimeError("Unable to find version string.")

version = find_version('GrowYourIC', '__init__.py')

metadata = dict(name="GrowYourIC",
                version=version,
                description="a toolkit to propagate seismic rays through models of Earth's inner core",
                url="https://github.com/MarineLasbleis/GrowYourIC",
                licence='GPL',
                long_description="a toolkit to propagate seismic rays through models of Earth's inner core",
                packages=['GrowYourIC'],
                classifiers=[
                    'Programming Language :: Python :: 3.4'],
                author='Marine Lasbleis',
                author_email='marine.lasleis@elsi.jp',
               )

try:
    from setuptools import setup
    metadata['install_requires'] = ['numpy', 'matplotlib', 'scipy']
except ImportError:
    from distutils.core import setup

setup(**metadata)


