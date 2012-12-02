#!/bin/bash
#
g++ -c -g -I/$HOME/include test_matrix_exponential_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_matrix_exponential_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ test_matrix_exponential_prb.o /$HOME/libcpp/$ARCH/test_matrix_exponential.o \
                                  /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_matrix_exponential_prb.o."
  exit
fi
#
rm test_matrix_exponential_prb.o
#
mv a.out test_matrix_exponential_prb
./test_matrix_exponential_prb > test_matrix_exponential_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_matrix_exponential_prb."
  exit
fi
rm test_matrix_exponential_prb
#
echo "Program output written to test_matrix_exponential_prb_output.txt"
