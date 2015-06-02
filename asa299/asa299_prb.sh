#!/bin/bash
#
g++ -c -I/$HOME/include asa299_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa299_prb.cpp"
  exit
fi
#
g++ asa299_prb.o /$HOME/libcpp/$ARCH/asa299.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa299_prb.o."
  exit
fi
#
rm asa299_prb.o
#
mv a.out asa299_prb
./asa299_prb > asa299_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa299_prb."
  exit
fi
rm asa299_prb
#
echo "Program output written to asa299_prb_output.txt"
