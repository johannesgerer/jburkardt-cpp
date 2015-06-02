#!/bin/bash
#
g++ -c -I/$HOME/include cg_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_prb.cpp"
  exit
fi
#
g++ cg_prb.o /$HOME/libcpp/$ARCH/cg.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cg_prb.o"
  exit
fi
#
rm cg_prb.o
#
mv a.out cg_prb
./cg_prb > cg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cg_prb."
  exit
fi
rm cg_prb
#
echo "Program output written to cg_prb_output.txt"
