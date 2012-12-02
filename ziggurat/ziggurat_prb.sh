#!/bin/bash
#
g++ -c -g -I/$HOME/include ziggurat_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ ziggurat_prb.o /$HOME/libcpp/$ARCH/ziggurat.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ziggurat_prb.o."
  exit
fi
#
rm ziggurat_prb.o
#
mv a.out ziggurat_prb
./ziggurat_prb > ziggurat_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ziggurat_prb."
  exit
fi
rm ziggurat_prb
#
echo "Program output written to ziggurat_prb_output.txt"
