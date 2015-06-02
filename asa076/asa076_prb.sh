#!/bin/bash
#
g++ -c -I/$HOME/include asa076_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling asa076_prb.cpp"
  exit
fi
#
g++ asa076_prb.o /$HOME/libcpp/$ARCH/asa076.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa076_prb.o."
  exit
fi
#
rm asa076_prb.o
#
mv a.out asa076_prb
./asa076_prb > asa076_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa076_prb."
  exit
fi
rm asa076_prb
#
echo "Program output written to asa076_prb_output.txt"
