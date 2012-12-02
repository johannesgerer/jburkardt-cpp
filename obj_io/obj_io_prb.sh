#!/bin/bash
#
g++ -c -g -I/$HOME/include obj_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling obj_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ obj_io_prb.o /$HOME/libcpp/$ARCH/obj_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading obj_io_prb.o."
  exit
fi
#
rm obj_io_prb.o
#
mv a.out obj_io_prb
./obj_io_prb > obj_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running obj_io_prb."
  exit
fi
rm obj_io_prb
#
echo "Program output written to obj_io_prb_output.txt"
