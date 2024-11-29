#!/bin/bash

datdir=./

python timing_vpaa.py input/vpaa1_02.json 100000
python timing_vpaa.py input/vpaa1_03.json 100000
python timing_vpaa.py input/vpaa1_04.json 100000
python timing_vpaa.py input/vpaa1_05.json 100000
python timing_vpaa.py input/vpaa1_06.json 100000
python timing_vpaa.py input/vpaa1_07.json 100000
python timing_vpaa.py input/vpaa1_08.json 100000
python timing_vpaa.py input/vpaa1_09.json 100000
python timing_vpaa.py input/vpaa1_10.json 100000

for i in `seq -w 2 10`; do
  for j in `seq 2`; do
    python timing_vpaa.py input/vpaa${j}_${i}.json >> timing.log
    \mv timing_data_vpaa${j}_${i}.pkl.gz ${datdir}
  done
done
python timing_dummy.py input/vpaa1_0[2-9].json input/vpaa1_10.json >> timing_dummy.log
\mv timing_dummy_data.pkl.gz ${datdir}
python vpaa_plot.py ${datdir}/timing_data_vpaa*.pkl.gz ${datdir}/timing_dummy_data.pkl.gz
