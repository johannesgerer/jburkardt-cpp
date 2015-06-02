#!/bin/bash
#
g++ -c -I/$HOME/include asa147_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa147_prb.cpp"
  exit
fi
#
g++ asa147_prb.o /$HOME/libcpp/$ARCH/asa147.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa147_prb.o."
  exit
fi
#
rm asa147_prb.o
#
mv a.out asa147_prb
./asa147_prb > asa147_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa147_prb."
  exit
fi
rm asa147_prb
#
echo "Program output written to asa147_prb_output.txt"
