#!/bin/bash
#
g++ -c -g -I/$HOME/include divdif_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling divdif_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ divdif_prb.o /$HOME/libcpp/$ARCH/divdif.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading divdif_prb.o."
  exit
fi
#
rm divdif_prb.o
#
mv a.out divdif_prb
./divdif_prb > divdif_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running divdif_prb."
  exit
fi
rm divdif_prb
#
echo "Program output written to divdif_prb_output.txt"
