#!/bin/bash
#
g++ -c -I/$HOME/include hermite_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_prb.cpp"
  exit
fi
#
g++ hermite_prb.o /$HOME/libcpp/$ARCH/hermite.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_prb.o."
  exit
fi
#
rm hermite_prb.o
#
mv a.out hermite_prb
./hermite_prb > hermite_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hermite_prb."
  exit
fi
rm hermite_prb
#
echo "Program output written to hermite_prb_output.txt"
