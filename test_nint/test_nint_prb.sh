#!/bin/bash
#
g++ -c -g -I/$HOME/include test_nint_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nint_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ test_nint_prb.o /$HOME/libcpp/$ARCH/test_nint.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nint_prb.o."
  exit
fi
#
rm test_nint_prb.o
#
mv a.out test_nint_prb
./test_nint_prb > test_nint_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_nint_prb."
  exit
fi
rm test_nint_prb
#
echo "Program output written to test_nint_prb_output.txt"
