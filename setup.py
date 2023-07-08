#!/usr/bin/env python3

from setuptools import setup

# setup is the gateway to the package build process.
# The only required components for a package are
# the name, author and contact, and description.
setup(
    name='gluex_moments_analysis',
    version='0.1.2',
    author='Richard Jones',
    author_email='richard.t.jones@uconn.edu',
    description='toolkit for implementation of support-vector moments analysis',
    license='GPL',
    url='https://github.com/rjones30/gluex_moments_analysis',
    packages=['gluex_moments_analysis'],
    package_dir={"src": "src"},
    package_data={"src": ["*.C", "*.h", "*.c", "Makefile"]},
    install_requires=[
            'numpy',
            'awkward',
            'uproot',
            'h5py',
            'scipy',
    ]
)
