#!/bin/bash
#
g++ -c -g -I/$HOME/include timestamp_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timestamp_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ timestamp_prb.o /$HOME/libcpp/$ARCH/timestamp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading timestamp_prb.o."
  exit
fi
#
rm timestamp_prb.o
#
mv a.out timestamp_prb
./timestamp_prb > timestamp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running timestamp_prb."
  exit
fi
rm timestamp_prb
#
echo "Program output written to timestamp_prb_output.txt"
