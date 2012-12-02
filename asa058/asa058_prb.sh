#!/bin/bash
#
g++ -c -g -I/$HOME/include asa058_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa058_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ asa058_prb.o /$HOME/libcpp/$ARCH/asa058.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa058_prb.o."
  exit
fi
#
rm asa058_prb.o
#
mv a.out asa058_prb
./asa058_prb > asa058_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa058_prb."
  exit
fi
rm asa058_prb
#
echo "Program output written to asa058_prb_output.txt"
