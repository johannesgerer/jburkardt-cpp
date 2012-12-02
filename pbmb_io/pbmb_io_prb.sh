#!/bin/bash
#
g++ -c -g -I/$HOME/include pbmb_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbmb_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pbmb_io_prb.o /$HOME/libcpp/$ARCH/pbmb_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbmb_io_prb.o."
  exit
fi
#
rm pbmb_io_prb.o
#
mv a.out pbmb_io_prb
./pbmb_io_prb > pbmb_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pbmb_io_prb."
  exit
fi
rm pbmb_io_prb
#
echo "Program output written to pbmb_io_prb_output.txt"
