#!/usr/bin/env python3

import setuptools


with open("README.md", 'r') as f:
    long_description = f.read()

setuptools.setup(
        name="martinitoolspymol",
        version="",
        author="Lorenzo Gaifas",
        author_email="brisvag@gmail.com",
        description="Collection of tools to visualize Martini coarse-grained trajectories in PyMOL",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/martinitoolspymol/mtools",
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
            "Operating System :: OS Independent",
            ],
        install_requires=[
            'networkx',
            'numpy'
            ],
        python_requires='>=3.5'
        )
