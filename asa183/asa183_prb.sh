#!/bin/bash
#
g++ -c -I/$HOME/include asa183_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa183_prb.cpp"
  exit
fi
#
g++ asa183_prb.o /$HOME/libcpp/$ARCH/asa183.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa183_prb.o."
  exit
fi
#
rm asa183_prb.o
#
mv a.out asa183_prb
./asa183_prb > asa183_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa183_prb."
  exit
fi
rm asa183_prb
#
echo "Program output written to asa183_prb_output.txt"
