#!/bin/bash
#
g++ -c -I/$HOME/include exactness_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling exactness_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ exactness_prb.o /$HOME/libcpp/$ARCH/exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading exactness_prb.o."
  exit
fi
#
rm exactness_prb.o
#
mv a.out exactness_prb
./exactness_prb > exactness_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running exactness_prb."
  exit
fi
rm exactness_prb
#
echo "Program output written to exactness_prb_output.txt"
