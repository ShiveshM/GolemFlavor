#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
from setuptools import setup, find_packages

NAME = 'GolemFlavor'
DESCRIPTION = 'GolemFlavor: A Python package for Astrophysical Flavor analysis with GolemFit'
MAINTAINER = 'Shivesh Mandalia'
MAINTAINER_EMAIL = 's.p.mandalia@qmul.ac.uk'
URL = 'https://github.com/ShiveshM/GolemFlavor'
LICENSE = 'MIT'

here = os.path.abspath(os.path.dirname(__file__))

def read(path, encoding='utf-8'):
    with io.open(path, encoding=encoding) as f:
        content = f.read()
    return content

def get_install_requirements(path):
    content = read(path)
    requirements = [req for req in content.split("\n")
                    if req != '' and not req.startswith('#')]
    return requirements

LONG_DESCRIPTION = read(os.path.join(here,'README.md'))

# Want to read in package version number from __version__.py
about = {}
with io.open(os.path.join(here, 'golemflavor', '__version__.py'), encoding='utf-8') as f:
    exec(f.read(), about)
    VERSION = about['__version__']

INSTALL_REQUIRES = get_install_requirements(os.path.join(here, 'requirements.txt'))

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url=URL,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    license=LICENSE,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    packages=find_packages(),
    install_requires=INSTALL_REQUIRES,
    setup_requires=['setuptools>=38.6.0'],
    scripts=[
        'scripts/contour.py',
        'scripts/fr.py',
        'scripts/mc_texture.py',
        'scripts/mc_unitary.py',
        'scripts/mc_x.py',
        'scripts/plot_sens.py',
        'scripts/sens.py'
    ]
)
