#!/bin/bash
#
g++ -c -g triangle_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ triangle_io_prb.o /$HOME/libcpp/$ARCH/triangle_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_io_prb.o."
  exit
fi
#
rm triangle_io_prb.o
#
mv a.out triangle_io_prb
./triangle_io_prb > triangle_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle_io_prb."
  exit
fi
rm triangle_io_prb
#
echo "Program output written to triangle_io_prb_output.txt"
