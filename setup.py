#!/usr/bin/env python

from setuptools import setup

setup(name='gmo_labeler',
      version='0.1.53',
      description='GMO labeler using hyperspectral image data and machine/deep learning',
      author='Michael Nagle',
      author_email='michael.nagle@oregonstate.edu',
      packages=['gmo_labeler'],
      install_requires=['pandas', 'numpy', 'PIL'])
