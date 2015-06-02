#!/bin/bash
#
g++ -c -g -I/$HOME/include ou_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ou_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ou_prb.o /$HOME/libcpp/$ARCH/ou.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ou_prb.o."
  exit
fi
#
rm ou_prb.o
#
mv a.out ou_prb
./ou_prb > ou_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ou_prb."
  exit
fi
rm ou_prb
#
echo "Program output written to ou_prb_output.txt"
