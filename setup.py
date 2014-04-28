#!/usr/bin/env python

# NOTE: multiprocessing import required for issues with nose tests.
# See: http://bugs.python.org/issue15881#msg170215
import multiprocessing

from setuptools import setup

setup(
    name='reference_genome_maker',
    version='0.1',
    author='Church Lab',
    author_email='gleb@mit.edu',
    maintainer='Gleb Kuznetsov',
    maintainer_email='gleb@mit.edu',
    url='http://churchlab.github.io/millstone/',
    package_dir={'': 'src'},
    packages=['reference_genome_maker'],
    install_requires=[
        'biopython >= 1.6.1',
        'PyVCF >= 0.6.7'
    ],
    test_suite='nose.collector',
)
