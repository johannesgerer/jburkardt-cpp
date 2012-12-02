#!/bin/bash
#
g++ -c -g -I/$HOME/include i8lib_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i8lib_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ i8lib_prb.o /$HOME/libcpp/$ARCH/i8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading i8lib_prb.o."
  exit
fi
#
rm i8lib_prb.o
#
mv a.out i8lib_prb
./i8lib_prb > i8lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running i8lib_prb."
  exit
fi
rm i8lib_prb
#
echo "Program output written to i8lib_prb_output.txt"
