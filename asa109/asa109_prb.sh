#!/bin/bash
#
g++ -c -I/$HOME/include asa109_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa109_prb.cpp"
  exit
fi
#
g++ asa109_prb.o /$HOME/libcpp/$ARCH/asa109.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa109_prb.o."
  exit
fi
#
rm asa109_prb.o
#
mv a.out asa109_prb
./asa109_prb > asa109_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa109_prb."
  exit
fi
rm asa109_prb
#
echo "Program output written to asa109_prb_output.txt"
