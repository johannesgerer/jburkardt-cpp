#!/bin/bash
#
g++ -c -g -I/$HOME/include chebyshev_interp_1d_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_interp_1d_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ chebyshev_interp_1d_prb.o /$HOME/libcpp/$ARCH/chebyshev_interp_1d.o \
                              /$HOME/libcpp/$ARCH/test_interp.o \
                              /$HOME/libcpp/$ARCH/qr_solve.o \
                              /$HOME/libcpp/$ARCH/r8lib.o -lm 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_interp_1d_prb.o."
  exit
fi
#
rm chebyshev_interp_1d_prb.o
#
mv a.out chebyshev_interp_1d_prb
./chebyshev_interp_1d_prb > chebyshev_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_interp_1d_prb."
  exit
fi
rm chebyshev_interp_1d_prb
#
echo "Program output written to chebyshev_interp_1d_prb_output.txt"
