#!/bin/bash
#
g++ -c -g -I/$HOME/include test_eigen_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_eigen_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ test_eigen_prb.o /$HOME/libcpp/$ARCH/test_eigen.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_eigen_prb.o."
  exit
fi
#
rm test_eigen_prb.o
#
mv a.out test_eigen_prb
./test_eigen_prb > test_eigen_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_eigen_prb."
  exit
fi
rm test_eigen_prb
#
echo "Program output written to test_eigen_prb_output.txt"
