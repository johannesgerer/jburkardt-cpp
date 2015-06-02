#!/bin/bash
#
g++ -c -I/$HOME/include c8lib_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling c8lib_prb.cpp"
  exit
fi
#
g++ c8lib_prb.o /$HOME/libcpp/$ARCH/c8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c8lib_prb.o."
  exit
fi
#
rm c8lib_prb.o
#
mv a.out c8lib_prb
./c8lib_prb > c8lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c8lib_prb."
  exit
fi
rm c8lib_prb
#
echo "Program output written to c8lib_prb_output.txt"
