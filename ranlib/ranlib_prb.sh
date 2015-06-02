#!/bin/bash
#
g++ -c -I/$HOME/include ranlib_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling ranlib_prb.cpp"
  exit
fi
#
g++ ranlib_prb.o /$HOME/libcpp/$ARCH/ranlib.o /$HOME/libcpp/$ARCH/rnglib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ranlib_prb.o."
  exit
fi
#
rm ranlib_prb.o
#
mv a.out ranlib_prb
./ranlib_prb > ranlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ranlib_prb."
  exit
fi
rm ranlib_prb
#
echo "Program output written to ranlib_prb_output.txt"
