#!/bin/bash
#
g++ -c -g -I/$HOME/include grf_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grf_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ grf_io_prb.o /$HOME/libcpp/$ARCH/grf_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grf_io_prb.o."
  exit
fi
#
rm grf_io_prb.o
#
mv a.out grf_io_prb
./grf_io_prb > grf_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running grf_io_prb."
  exit
fi
rm grf_io_prb
#
echo "Program output written to grf_io_prb_output.txt"
