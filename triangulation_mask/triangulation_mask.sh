#!/bin/bash
#
g++ -c -g triangulation_mask.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_mask.cpp"
  exit
fi
rm compiler.txt
#
mv triangulation_mask.o ~/libcpp/$ARCH/
#
echo "Library installed as ~/libcpp/$ARCH/triangulation_mask.o"
