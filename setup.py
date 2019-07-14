#!/usr/bin/env python3

import setuptools


with open('README.md', 'r') as f:
    long_description = f.read()

setuptools.setup(
        name='stir',
        scripts=['bin/stir'],
        version='0.1-alpha',
        author='Lorenzo Gaifas',
        author_email='brisvag@gmail.com',
        description='Wrapper for a collection of tools to visualize Martini coarse-grained trajectories in PyMOL',
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/brisvag/stir',
        packages=setuptools.find_packages(),
        classifiers=[
            'Programming Language :: Python :: 3',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
            'Operating System :: OS Independent',
        ],
        install_requires=[
            'pymol',
            'garnish',
            'psutil',
        ],
        python_requires='>=3.5'
        )
