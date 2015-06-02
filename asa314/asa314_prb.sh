#!/bin/bash
#
g++ -c -I/$HOME/include asa314_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa314_prb.cpp"
  exit
fi
#
g++ asa314_prb.o /$HOME/libcpp/$ARCH/asa314.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa314_prb.o."
  exit
fi
#
rm asa314_prb.o
#
mv a.out asa314_prb
./asa314_prb > asa314_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa314_prb."
  exit
fi
rm asa314_prb
#
echo "Program output written to asa314_prb_output.txt"
