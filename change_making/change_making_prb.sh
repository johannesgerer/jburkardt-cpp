#!/bin/bash
#
g++ -c -I/$HOME/include change_making_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling change_making_prb.cpp"
  exit
fi
#
g++ -o change_making_prb change_making_prb.o /$HOME/libcpp/$ARCH/change_making.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading change_making_prb.o."
  exit
fi
#
rm change_making_prb.o
#
./change_making_prb > change_making_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running change_making_prb."
  exit
fi
rm change_making_prb
#
echo "Program output written to change_making_prb_output.txt"
