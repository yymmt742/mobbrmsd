#!/bin/bash

python sample.py

python sample_stats.py sample_03_*npy
python sample_stats.py sample_04_*npy
python sample_stats.py sample_05_*npy
python sample_stats.py sample_06_*npy
python sample_stats.py sample_07_*npy
python sample_stats.py sample_08_*npy
python sample_stats.py sample_09_*npy
python sample_stats.py sample_10_*npy

python stats_plot.py sample_*.gz
python stats_fit.py sample_*.gz > stats_fit.log
python vs_lprmsd.py > vs_lprmsd.log
