#!/bin/bash
#
g++ -c -g -I/$HOME/include pwl_interp_2d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pwl_interp_2d_prb.o /$HOME/libcpp/$ARCH/pwl_interp_2d.o \
                        /$HOME/libcpp/$ARCH/test_interp_2d.o \
                        /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pwl_interp_2d_prb.o."
  exit
fi
#
rm pwl_interp_2d_prb.o
#
mv a.out pwl_interp_2d_prb
./pwl_interp_2d_prb > pwl_interp_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pwl_interp_2d_prb."
  exit
fi
rm pwl_interp_2d_prb
#
echo "Program output written to pwl_interp_2d_prb_output.txt"
