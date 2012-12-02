#!/bin/bash
#
g++ -c -g -I/$HOME/include test_optimization_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_optimization_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ test_optimization_prb.o /$HOME/libcpp/$ARCH/test_optimization.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_optimization_prb.o."
  exit
fi
#
rm test_optimization_prb.o
#
mv a.out test_optimization_prb
./test_optimization_prb > test_optimization_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_optimization_prb."
  exit
fi
rm test_optimization_prb
#
echo "Program output written to test_optimization_prb_output.txt"
