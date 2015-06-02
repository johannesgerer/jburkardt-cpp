#!/bin/bash
#
g++ -c -I/$HOME/include asa144_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa144_prb.cpp"
  exit
fi
#
g++ asa144_prb.o /$HOME/libcpp/$ARCH/asa144.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa144_prb.o."
  exit
fi
#
rm asa144_prb.o
#
mv a.out asa144_prb
./asa144_prb > asa144_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa144_prb."
  exit
fi
rm asa144_prb
#
echo "Program output written to asa144_prb_output.txt"
