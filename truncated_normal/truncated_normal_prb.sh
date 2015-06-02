#!/bin/bash
#
g++ -c -I/$HOME/include truncated_normal_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling truncated_normal_prb.cpp"
  exit
fi
#
g++ truncated_normal_prb.o /$HOME/libcpp/$ARCH/truncated_normal.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading truncated_normal_prb.o."
  exit
fi
#
rm truncated_normal_prb.o
#
mv a.out truncated_normal_prb
./truncated_normal_prb > truncated_normal_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running truncated_normal_prb."
  exit
fi
rm truncated_normal_prb
#
echo "Program output written to truncated_normal_prb_output.txt"
