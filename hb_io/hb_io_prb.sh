#!/bin/bash
#
g++ -c -I/$HOME/include hb_io_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hb_io_prb.cpp"
  exit
fi
#
g++ hb_io_prb.o /$HOME/libcpp/$ARCH/hb_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hb_io_prb.o."
  exit
fi
#
rm hb_io_prb.o
#
mv a.out hb_io_prb
./hb_io_prb > hb_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hb_io_prb."
  exit
fi
rm hb_io_prb
#
echo "Program output written to hb_io_prb_output.txt"
