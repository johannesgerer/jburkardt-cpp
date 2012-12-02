#!/bin/bash
#
g++ -c -g -I/$HOME/include quality_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quality_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ quality_prb.o /$HOME/libcpp/$ARCH/quality.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quality_prb.o."
  exit
fi
#
rm quality_prb.o
#
mv a.out quality_prb
./quality_prb > quality_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quality_prb."
  exit
fi
rm quality_prb
#
echo "Program output written to quality_prb_output.txt"
