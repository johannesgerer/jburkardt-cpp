#!/bin/bash
#
g++ -c -g -I/$HOME/include i4lib_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i4lib_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ i4lib_prb.o /$HOME/libcpp/$ARCH/i4lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading i4lib_prb.o."
  exit
fi
#
rm i4lib_prb.o
#
mv a.out i4lib_prb
./i4lib_prb > i4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running i4lib_prb."
  exit
fi
rm i4lib_prb
#
echo "Program output written to i4lib_prb_output.txt"
