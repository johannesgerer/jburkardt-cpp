#!/bin/bash
#
g++ -c -g -I/$HOME/include fem_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ fem_io_prb.o /$HOME/libcpp/$ARCH/fem_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_io_prb.o."
  exit
fi
#
rm fem_io_prb.o
#
mv a.out fem_io_prb
./fem_io_prb > fem_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem_io_prb."
  exit
fi
rm fem_io_prb
#
echo "Program output written to fem_io_prb_output.txt"
