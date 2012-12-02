#!/bin/bash
#
g++ -c -g -I/$HOME/include keast_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling keast_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ keast_prb.o /$HOME/libcpp/$ARCH/keast.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading keast_prb.o."
  exit
fi
#
rm keast_prb.o
#
mv a.out keast_prb
./keast_prb > keast_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running keast_prb."
  exit
fi
rm keast_prb
#
echo "Program output written to keast_prb_output.txt"
