#!/bin/bash
#
g++ -c -g -I/$HOME/include xy_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xy_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ xy_io_prb.o /$HOME/libcpp/$ARCH/xy_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xy_io_prb.o."
  exit
fi
#
rm xy_io_prb.o
#
mv a.out xy_io_prb
./xy_io_prb > xy_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xy_io_prb."
  exit
fi
rm xy_io_prb
#
echo "Program output written to xy_io_prb_output.txt"
