#!/bin/bash
#
g++ -c -g -I/$HOME/include qwv_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ qwv_prb.o /$HOME/libcpp/$ARCH/qwv.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qwv_prb.o."
  exit
fi
#
rm qwv_prb.o
#
mv a.out qwv_prb
./qwv_prb > qwv_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qwv_prb."
  exit
fi
rm qwv_prb
#
echo "Program output written to qwv_prb_output.txt"
