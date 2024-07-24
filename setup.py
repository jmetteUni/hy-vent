#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 17:21:20 2024

@author: jmette
"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='hy_vent',
      version='0.1',
      description='Python package to analyse CTD Data of hydrothermal venting ',
#      long_description=readme(),
      url='https://github.com/jmetteUni/hy-vent',
      author='Jonathan Mette',
      author_email='j.mette@posteo.de',
      license='MIT',
      packages=['hy_vent'],
      install_requires=[
          'os',
          'pandas',
          'datetime'],
      zip_safe=False)


