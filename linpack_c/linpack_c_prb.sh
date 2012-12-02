#!/bin/bash
#
g++ -c -g -I/$HOME/include linpack_c_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling linpack_c_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ linpack_c_prb.o /$HOME/libcpp/$ARCH/linpack_c.o /$HOME/libcpp/$ARCH/blas1_c.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linpack_c_prb.o."
  exit
fi
#
rm linpack_c_prb.o
#
mv a.out linpack_c_prb
./linpack_c_prb > linpack_c_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linpack_c_prb."
  exit
fi
rm linpack_c_prb
#
echo "Program output written to linpack_c_prb_output.txt"
