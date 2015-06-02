#!/bin/bash
#
g++ -c -I/$HOME/include asa007_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa007_prb.cpp"
  exit
fi
#
g++ asa007_prb.o /$HOME/libcpp/$ARCH/asa007.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa007_prb.o."
  exit
fi
#
rm asa007_prb.o
#
mv a.out asa007_prb
./asa007_prb > asa007_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa007_prb."
  exit
fi
rm asa007_prb
#
echo "Program output written to asa007_prb_output.txt"
