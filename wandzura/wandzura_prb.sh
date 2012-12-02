#!/bin/bash
#
g++ -c -g -I/$HOME/include wandzura_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wandzura_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ wandzura_prb.o /$HOME/libcpp/$ARCH/wandzura.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wandzura_prb.o."
  exit
fi
#
rm wandzura_prb.o
#
mv a.out wandzura_prb
./wandzura_prb > wandzura_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wandzura_prb."
  exit
fi
rm wandzura_prb
#
echo "Program output written to wandzura_prb_output.txt"
