#!/bin/bash
#
g++ -c -g -I/$HOME/include minpack_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling minpack_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ minpack_prb.o /$HOME/libcpp/$ARCH/minpack.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading minpack_prb.o."
  exit
fi
#
rm minpack_prb.o
#
mv a.out minpack_prb
./minpack_prb > minpack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running minpack_prb."
  exit
fi
rm minpack_prb
#
echo "Program output written to minpack_prb_output.txt"
