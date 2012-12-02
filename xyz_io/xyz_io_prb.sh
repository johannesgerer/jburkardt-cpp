#!/bin/bash
#
g++ -c -g -I/$HOME/include xyz_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xyz_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ xyz_io_prb.o /$HOME/libcpp/$ARCH/xyz_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xyz_io_prb.o."
  exit
fi
#
rm xyz_io_prb.o
#
mv a.out xyz_io_prb
./xyz_io_prb > xyz_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running xyz_io_prb."
  exit
fi
rm xyz_io_prb
#
echo "Program output written to xyz_io_prb_output.txt"
