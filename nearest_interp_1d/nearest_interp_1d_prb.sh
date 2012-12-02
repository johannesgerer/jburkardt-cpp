#!/bin/bash
#
g++ -c -g -I/$HOME/include nearest_interp_1d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nearest_interp_1d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ nearest_interp_1d_prb.o /$HOME/libcpp/$ARCH/nearest_interp_1d.o \
                            /$HOME/libcpp/$ARCH/test_interp.o \
                            /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nearest_interp_1d_prb.o"
  exit
fi
#
rm nearest_interp_1d_prb.o
#
mv a.out nearest_interp_1d_prb
./nearest_interp_1d_prb > nearest_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nearest_interp_1d_prb."
  exit
fi
rm nearest_interp_1d_prb
#
echo "Program output written to nearest_interp_1d_prb_output.txt"
