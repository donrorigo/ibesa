#!/bin/bash

function rep()
{
while read field1 field2 field3; do
  ./ibesa 100 100 10 $field1 $field2 $field3
  ./IBESA\ convert\ to\ HYPERVOLUME $field1.txt
  ./normalize RESULTADOS\:\ $field1 $field1
  ./hyp_ind hyp_ind_param.txt NORMALIZADO\:\ $field1 null.txt output.txt
  echo "$field1:$(./convert2decimal output.txt)" >> "hypervolume_rep$1.txt"
  rm $field1.txt RESULTADOS\:\ $field1 NORMALIZADO\:\ $field1 output.txt
done < proteins.txt
}

for i in {1..5}
do
rep i
done