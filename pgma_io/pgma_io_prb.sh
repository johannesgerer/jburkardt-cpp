#!/bin/bash
#
g++ -c -g -I/$HOME/include pgma_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgma_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pgma_io_prb.o /$HOME/libcpp/$ARCH/pgma_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pgma_io_prb.o."
  exit
fi
#
rm pgma_io_prb.o
#
mv a.out pgma_io_prb
./pgma_io_prb > pgma_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pgma_io_prb."
  exit
fi
rm pgma_io_prb
#
echo "Program output written to pgma_io_prb_output.txt"
