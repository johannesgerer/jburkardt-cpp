#!/bin/bash
#
g++ -c -I/$HOME/include cc_io_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_io_prb.cpp"
  exit
fi
#
g++ cc_io_prb.o /$HOME/libcpp/$ARCH/cc_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cc_io_prb.o"
  exit
fi
#
rm cc_io_prb.o
#
mv a.out cc_io_prb
./cc_io_prb > cc_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cc_io_prb."
  exit
fi
rm cc_io_prb
#
echo "Program output written to cc_io_prb_output.txt"
