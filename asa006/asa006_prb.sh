#!/bin/bash
#
g++ -c -I/$HOME/include asa006_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa006_prb.cpp"
  exit
fi
#
g++ asa006_prb.o /$HOME/libcpp/$ARCH/asa006.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa006_prb.o."
  exit
fi
#
rm asa006_prb.o
#
mv a.out asa006_prb
./asa006_prb > asa006_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa006_prb."
  exit
fi
rm asa006_prb
#
echo "Program output written to asa006_prb_output.txt"
