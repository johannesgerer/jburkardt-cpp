#!/bin/bash
#
g++ -c -I/$HOME/include qls_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qls_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ qls_prb.o /$HOME/libcpp/$ARCH/qls.o /$HOME/libcpp/$ARCH/qr_solve.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qls_prb.o."
  exit
fi
#
rm qls_prb.o
#
mv a.out qls_prb
./qls_prb > qls_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qls_prb."
  exit
fi
rm qls_prb
#
echo "Program output written to qls_prb_output.txt"
