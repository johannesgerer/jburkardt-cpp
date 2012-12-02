#!/bin/bash
#
g++ -c -g -I/$HOME/include pwl_approx_1d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_approx_1d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pwl_approx_1d_prb.o /$HOME/libcpp/$ARCH/pwl_approx_1d.o \
                        /$HOME/libcpp/$ARCH/test_interp_1d.o \
                        /$HOME/libcpp/$ARCH/qr_solve.o \
                        /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pwl_approx_1d_prb.o."
  exit
fi
#
rm pwl_approx_1d_prb.o
#
mv a.out pwl_approx_1d_prb
./pwl_approx_1d_prb > pwl_approx_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pwl_approx_1d_prb."
  exit
fi
rm pwl_approx_1d_prb
#
echo "Program output written to pwl_approx_1d_prb_output.txt"
