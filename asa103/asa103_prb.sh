#!/bin/bash
#
g++ -c -I/$HOME/include asa103_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa103_prb.cpp"
  exit
fi
#
g++ asa103_prb.o /$HOME/libcpp/$ARCH/asa103.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa103_prb.o."
  exit
fi
#
rm asa103_prb.o
#
mv a.out asa103_prb
./asa103_prb > asa103_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa103_prb."
  exit
fi
rm asa103_prb
#
echo "Program output written to asa103_prb_output.txt"
