#!/bin/bash
#
g++ -c -g p15_mask.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling p15_mask.cpp"
  exit
fi
rm compiler.txt
#
g++ ~/libcpp/$ARCH/triangulation_mask.o p15_mask.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_mask.o + p15_mask.o"
  exit
fi
rm p15_mask.o
#
chmod ugo+x a.out
mv a.out triangulation_mask
./triangulation_mask p15 > p15_mask_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangulation_mask."
  exit
fi
rm triangulation_mask
#
echo "Program output written to p15_mask_output.txt"
