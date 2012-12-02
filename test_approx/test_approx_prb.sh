#!/bin/bash
#
g++ -c -g -I/$HOME/include test_approx_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_approx_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ test_approx_prb.o /$HOME/libcpp/$ARCH/test_approx.o /$HOME/libcpp/$ARCH/spline.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_approx_prb.o."
  exit
fi
#
rm test_approx_prb.o
#
mv a.out test_approx_prb
./test_approx_prb > test_approx_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_approx_prb."
  exit
fi
rm  test_approx_prb
#
echo "Program output written to test_approx_prb_output.txt"
