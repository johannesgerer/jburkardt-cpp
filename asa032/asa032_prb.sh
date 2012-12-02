#!/bin/bash
#
g++ -c -g -I/$HOME/include asa032_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa032_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa032_prb.o /$HOME/libcpp/$ARCH/asa032.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa032_prb.o."
  exit
fi
#
rm asa032_prb.o
#
mv a.out asa032_prb
./asa032_prb > asa032_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa032_prb."
  exit
fi
rm asa032_prb
#
echo "Program output written to asa032_prb_output.txt"
