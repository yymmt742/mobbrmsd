#!/usr/bin/env python
# -*- coding: utf-8 -*-
from skbuild import setup

setup(
    name='mobbrmsd',
    version="0.0.1",
    author="yymmt742",
    author_email="yymmt@kuchem.kyoto-u.ac.jp",
    packages=["mobbrmsd"],
    cmake_install_dir="mobbrmsd/bin",
    url="https://github.com/yymmt742/mobbrmsd",
    download_url="https://github.com/yymmt742/mobbrmsd",
    description="molecular orientation corrected rmsd with branch-and-bound",
    long_description="molecular orientation corrected rmsd with branch-and-bound",
    classifiers=[
                'Programming Language :: Fortran',
                'Programming Language :: Python'
                ],
    license="MIT",
    python_requires=">=3.8"
    )
