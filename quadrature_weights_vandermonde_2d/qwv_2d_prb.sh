#!/bin/bash
#
g++ -c -I/$HOME/include qwv_2d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv_2d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ qwv_2d_prb.o /$HOME/libcpp/$ARCH/qwv_2d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qwv_2d_prb.o."
  exit
fi
#
rm qwv_2d_prb.o
#
mv a.out qwv_2d_prb
./qwv_2d_prb > qwv_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qwv_2d_prb."
  exit
fi
rm qwv_2d_prb
#
echo "Program output written to qwv_2d_prb_output.txt"
