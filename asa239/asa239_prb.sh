#!/bin/bash
#
g++ -c -I/$HOME/include asa239_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa239_prb.cpp"
  exit
fi
#
g++ asa239_prb.o /$HOME/libcpp/$ARCH/asa239.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa239_prb.o."
  exit
fi
#
rm asa239_prb.o
#
mv a.out asa239_prb
./asa239_prb > asa239_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa239_prb."
  exit
fi
rm asa239_prb
#
echo "Program output written to asa239_prb_output.txt"
