#!/usr/bin/env python3

from setuptools import setup

# setup is the gateway to the package build process.
# The only required components for a package are
# the name, author and contact, and description.
setup(
    name='gluex_moments_analysis',
    version='0.2.03',
    author='Richard Jones',
    author_email='richard.t.jones@uconn.edu',
    description='toolkit for implementation of support-vector moments analysis',
    license='GPL',
    url='https://github.com/rjones30/gluex_moments_analysis',
    packages=['gluex_moments_analysis'],
    package_dir={"gluex_moments_analysis": "gluex_moments_analysis"},
    package_data={"gluex_moments_analysis": ["src/*.C", "src/*.h", "src/*.c", "src/Makefile"]},
    install_requires=[
            'numpy',
            'awkward',
            'uproot',
            'h5py',
            'scipy',
    ]
)
