#!/bin/bash
#
g++ -c -g -I/$HOME/include sftpack_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sftpack_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sftpack_prb.o /$HOME/libcpp/$ARCH/sftpack.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sftpack_prb.o."
  exit
fi
#
rm sftpack_prb.o
#
mv a.out sftpack_prb
./sftpack_prb > sftpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sftpack_prb."
  exit
fi
rm sftpack_prb
#
echo "Program output written to sftpack_prb_output.txt"
