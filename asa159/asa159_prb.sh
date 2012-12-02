#!/bin/bash
#
g++ -c -g -I/$HOME/include asa159_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa159_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa159_prb.o /$HOME/libcpp/$ARCH/asa159.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa159_prb.o."
  exit
fi
#
rm asa159_prb.o
#
mv a.out asa159_prb
./asa159_prb > asa159_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa159_prb."
  exit
fi
rm asa159_prb
#
echo "Program output written to asa159_prb_output.txt"
