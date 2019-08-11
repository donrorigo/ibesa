#!/bin/bash

while read field1 field2 field3; do
  ./ibesa 100 100 80 $field1 $field2 $field3
done < proteins.txt