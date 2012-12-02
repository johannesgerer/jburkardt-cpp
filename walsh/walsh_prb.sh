#!/bin/bash
#
g++ -c -g -I/$HOME/include walsh_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling walsh_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ walsh_prb.o /$HOME/libcpp/$ARCH/walsh.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading walsh_prb.o."
  exit
fi
#
rm walsh_prb.o
#
mv a.out walsh_prb
./walsh_prb > walsh_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running walsh_prb."
  exit
fi
rm walsh_prb
#
echo "Program output written to walsh_prb_output.txt"
