# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 12:58:48 2018

@author: babin
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='recan',
      version='0.1',
      author='Yuriy Babin',
      author_email='babin.yurii@gmail.com',
      description='recombination analysis tool',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/babinyurii/RECAN',
      packages=setuptools.find_packages(),
      classifiers=("Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",),)
    