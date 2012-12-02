#!/bin/bash
#
g++ -c -g -I/$HOME/include normal_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ normal_prb.o /$HOME/libcpp/$ARCH/normal.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading normal_prb.o."
  exit
fi
#
rm normal_prb.o
#
mv a.out normal_prb
./normal_prb > normal_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running normal_prb."
  exit
fi
rm normal_prb
#
echo "Program output written to normal_prb_output.txt"
