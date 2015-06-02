#!/bin/bash
#
g++ -c -I/$HOME/include hermite_cubic_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_cubic_prb.cpp"
  exit
fi
#
g++ hermite_cubic_prb.o /$HOME/libcpp/$ARCH/hermite_cubic.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_cubic_prb.o."
  exit
fi
#
rm hermite_cubic_prb.o
#
mv a.out hermite_cubic_prb
./hermite_cubic_prb > hermite_cubic_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hermite_cubic_prb."
  exit
fi
rm hermite_cubic_prb
#
echo "Program output written to hermite_cubic_prb_output.txt"
