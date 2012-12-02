#!/bin/bash
#
g++ -c -g -I/$HOME/include bmp_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bmp_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ bmp_io_prb.o /$HOME/libcpp/$ARCH/bmp_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bmp_io_prb.o."
  exit
fi
#
rm bmp_io_prb.o
#
mv a.out bmp_io_prb
./bmp_io_prb > bmp_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bmp_io_prb."
  exit
fi
rm bmp_io_prb
#
echo "Program output written to bmp_io_prb_output.txt"
