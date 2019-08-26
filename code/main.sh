#!/bin/bash

while read field1 field2 field3; do
  ./ibesa 100 100 10 $field1 $field2 $field3 
  ./normalize RESULTADOS\:\ $field1 $field1
  ./hyp_ind hyp_ind_param.txt NORMALIZADO\:\ $field1 afjj.txt "$field1_output.txt"
  echo $(./convert2decimal $field1_output.txt) >> "hypervolume: $field1.txt"
done < proteins.txt