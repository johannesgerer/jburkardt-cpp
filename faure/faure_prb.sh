#!/bin/bash
#
g++ -c -g -I/$HOME/include faure_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling faure_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ faure_prb.o /$HOME/libcpp/$ARCH/faure.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading faure_prb.o."
  exit
fi
#
rm faure_prb.o
#
mv a.out faure_prb
./faure_prb > faure_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running faure_prb."
  exit
fi
rm faure_prb
#
echo "Program output written to faure_prb_output.txt"
