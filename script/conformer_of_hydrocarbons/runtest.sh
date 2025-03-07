#!/bin/bash

sys=("ethane" "propane" "isobutane" "neopentane" "tetramethylbutane" "hexamethylpentane" "tetramethylhexane" "tetramethylpentane" "octamethylhexane")
n_sample=("10" "10" "10" "10" "10" "10" "10" "10" "10")
run_getbestrms=("1" "1" "1" "1" "1" "0" "0" "0" "0")

for i in "${!sys[@]}"; do
  echo ${sys[$i]}
  python gen_conformers.py ${sys[$i]} ${n_sample[$i]} | obabel -i gzmat -o pdb > ${sys[$i]}_ref.pdb
  python gen_conformers.py ${sys[$i]} ${n_sample[$i]} | obabel -i gzmat -o pdb > ${sys[$i]}_trg.pdb
  python conformer_of_hydrogens.py input/${sys[$i]}.json ${run_getbestrms[$i]}
done
