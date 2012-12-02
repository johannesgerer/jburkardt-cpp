#!/bin/bash
#
g++ -c -g -I/$HOME/include tanh_quad_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tanh_quad_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ tanh_quad_prb.o /$HOME/libcpp/$ARCH/tanh_quad.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tanh_quad_prb.o."
  exit
fi
#
rm tanh_quad_prb.o
#
mv a.out tanh_quad_prb
./tanh_quad_prb > tanh_quad_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tanh_quad_prb."
  exit
fi
rm tanh_quad_prb
#
echo "Program output written to tanh_quad_prb_output.txt"
