#!/bin/bash

sys=("water" "ammonia" "methane")

for i in "${!sys[@]}"; do
  python molecular_aggregates.py input/${sys[$i]}.json
done
