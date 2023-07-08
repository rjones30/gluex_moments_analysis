#!/usr/bin/env python3

from setuptools import setup

# setup is the gateway to the package build process.
# The only required components for a package are
# the name, author and contact, and description.
setup(
    name='gluex_moments_analysis',
    version='0.1.0',
    author='Richard Jones',
    author_email='richard.t.jones@uconn.edu',
    description='toolkit for implementation of support-vector moments analysis',
    license='GPL',
    url='https://github.com/rjones30/gluex_moments_analysis',
    packages=['gluex_moments_analysis'],
    install_requires=[
            'numpy',
            'math',
            'ROOT',
            'uproot',
            'h5py',
            'threading',
            'scipy',
            'time',
            'math'
    ]
)
