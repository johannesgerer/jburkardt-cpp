#!/bin/bash
#
g++ -c -g -I/$HOME/include ppmb_io_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ppmb_io_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ ppmb_io_prb.o /$HOME/libcpp/$ARCH/ppmb_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ppmb_io_prb.o."
  exit
fi
#
rm ppmb_io_prb.o
#
mv a.out ppmb_io_prb
./ppmb_io_prb > ppmb_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ppmb_io_prb."
  exit
fi
rm ppmb_io_prb
#
echo "Program output written to ppmb_io_prb_output.txt"
