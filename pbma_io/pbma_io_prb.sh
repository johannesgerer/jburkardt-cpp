#!/bin/bash
#
g++ -c -g -I/$HOME/include pbma_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbma_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pbma_io_prb.o /$HOME/libcpp/$ARCH/pbma_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbma_io_prb.o."
  exit
fi
#
rm pbma_io_prb.o
#
mv a.out pbma_io_prb
./pbma_io_prb > pbma_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pbma_io_prb."
  exit
fi
rm pbma_io_prb
#
echo "Program output written to pbma_io_prb_output.txt"
