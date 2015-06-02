#!/bin/bash
#
g++ -c -I/$HOME/include asa053_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa053_prb.cpp"
  exit
fi
#
g++ asa053_prb.o /$HOME/libcpp/$ARCH/asa053.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa053_prb.o."
  exit
fi
#
rm asa053_prb.o
#
mv a.out asa053_prb
./asa053_prb > asa053_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa053_prb."
  exit
fi
rm asa053_prb
#
echo "Program output written to asa053_prb_output.txt"
