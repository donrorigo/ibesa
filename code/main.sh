#!/bin/bash

mkdir "IBEA Results"
for i in {1..5}
do
  mkdir IBEA\ Results/Rep$i
  while read field1 field2 field3; do
    ./ibesa 100 100 10 $field1 $field2 $field3
    ./IBESA\ convert\ to\ HYPERVOLUME $field1.txt
    ./normalize RESULTADOS\:\ $field1 $field1 0.40
    ./hyp_ind hyp_ind_param.txt NORMALIZADO\:\ $field1 null.txt output.txt
    echo "$field1: $(./convert2decimal output.txt)" >> "hypervolumes.txt"
    echo "$field1: $(./convert2decimal output.txt)" >> $field1"_utility.txt"
    rm RESULTADOS\:\ $field1 output.txt
    mv $field1.txt NORMALIZADO\:\ $field1 IBEA\ Results/Rep$i
  done < proteins.txt
  mv "hypervolumes.txt" IBEA\ Results/Rep$i
done