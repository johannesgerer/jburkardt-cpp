#!/bin/bash
#
g++ -c -g -I/$HOME/include bernstein_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bernstein_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ bernstein_prb.o /$HOME/libcpp/$ARCH/bernstein.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bernstein_prb.o."
  exit
fi
#
rm bernstein_prb.o
#
mv a.out bernstein_prb
./bernstein_prb > bernstein_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bernstein_prb."
  exit
fi
rm bernstein_prb
#
echo "Program output written to bernstein_prb_output.txt"
