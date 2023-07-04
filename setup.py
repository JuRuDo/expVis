#!/bin/env python

#######################################################################
# Copyright (C) 2022 Julian Dosch
#
# This file is part of SpICE.
#
#  SpICE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SpICE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with expNet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="SpICE",
    version="0.2",
    python_requires='>=3.9.0',
    description="A dashboard to look at transcript expression weighted with FAS score",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Julian Dosch",
    author_email="Dosch@bio.uni-frankfurt.de",
    url="",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'dash',
        'dash_cytoscape',
        'dash_bootstrap_components',
        'colour',
        'numpy',
        'pandas',
        'pyyaml',
        'dash_daq',
        'scikit-learn',
        'scipy',
        'kaleido'
    ],
    entry_points={
        'console_scripts': [
            'spice.dash = expVis.ExpVis:main',
            'spice.pca = expVis.create_PCA_data:get_options'
        ]
    },
    license="GPL-3.0",
    classifiers=[
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
