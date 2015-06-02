#!/bin/bash
#
g++ -c -I/$HOME/include wathen_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling wathen_prb.cpp"
  exit
fi
#
g++ wathen_prb.o /$HOME/libcpp/$ARCH/wathen.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wathen_prb.o."
  exit
fi
#
rm wathen_prb.o
#
mv a.out wathen_prb
./wathen_prb > wathen_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wathen_prb."
  exit
fi
rm wathen_prb
#
echo "Program output written to wathen_prb_output.txt"
