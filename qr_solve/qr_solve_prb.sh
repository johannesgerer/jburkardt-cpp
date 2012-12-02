#!/bin/bash
#
g++ -c -g -I/$HOME/include qr_solve_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qr_solve_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ qr_solve_prb.o /$HOME/libcpp/$ARCH/qr_solve.o /$HOME/libcpp/$ARCH/test_ls.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qr_solve_prb.o."
  exit
fi
#
rm qr_solve_prb.o
#
mv a.out qr_solve_prb
./qr_solve_prb > qr_solve_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qr_solve_prb."
  exit
fi
rm qr_solve_prb
#
echo "Program output written to qr_solve_prb_output.txt"
