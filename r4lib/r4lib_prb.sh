#!/bin/bash
#
g++ -c -I/$HOME/include r4lib_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling r4lib_prb.cpp"
  exit
fi
#
g++ r4lib_prb.o /$HOME/libcpp/$ARCH/r4lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r4lib_prb.o."
  exit
fi
#
rm r4lib_prb.o
#
mv a.out r4lib_prb
./r4lib_prb > r4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running r4lib_prb."
  exit
fi
rm r4lib_prb
#
echo "Program output written to r4lib_prb_output.txt"
