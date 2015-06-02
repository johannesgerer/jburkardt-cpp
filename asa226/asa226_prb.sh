#!/bin/bash
#
g++ -c -I/$HOME/include asa226_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa226_prb.cpp"
  exit
fi
#
g++ asa226_prb.o /$HOME/libcpp/$ARCH/asa226.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa226_prb.o."
  exit
fi
#
rm asa226_prb.o
#
mv a.out asa226_prb
./asa226_prb > asa226_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa226_prb."
  exit
fi
rm asa226_prb
#
echo "Program output written to asa226_prb_output.txt"
