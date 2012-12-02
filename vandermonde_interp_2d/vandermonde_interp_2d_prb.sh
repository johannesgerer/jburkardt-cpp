#!/bin/bash
#
g++ -c -g -I/$HOME/include vandermonde_interp_2d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_interp_2d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ vandermonde_interp_2d_prb.o /$HOME/libcpp/$ARCH/vandermonde_interp_2d.o \
                                /$HOME/libcpp/$ARCH/qr_solve.o \
                                /$HOME/libcpp/$ARCH/test_interp_2d.o \
                                /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vandermonde_interp_2d_prb.o"
  exit
fi
#
rm vandermonde_interp_2d_prb.o
#
mv a.out vandermonde_interp_2d_prb
./vandermonde_interp_2d_prb > vandermonde_interp_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vandermonde_interp_2d_prb."
  exit
fi
rm vandermonde_interp_2d_prb
#
echo "Program output written to vandermonde_interp_2d_prb_output.txt"
