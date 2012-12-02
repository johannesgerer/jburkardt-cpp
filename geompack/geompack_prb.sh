#!/bin/bash
#
g++ -c -g -I/$HOME/include geompack_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geompack_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ geompack_prb.o /$HOME/libcpp/$ARCH/geompack.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading geompack_prb.o."
  exit
fi
#
rm geompack_prb.o
#
mv a.out geompack_prb
./geompack_prb > geompack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running geompack_prb."
  exit
fi
rm geompack_prb
#
echo "Program output written to geompack_prb_output.txt"
