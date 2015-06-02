#!/bin/bash
#
g++ -c -I/$HOME/include test_values_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling test_values_prb.cpp"
  exit
fi
#
g++ test_values_prb.o /$HOME/libcpp/$ARCH/test_values.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_values_prb.o."
  exit
fi
#
rm test_values_prb.o
#
mv a.out test_values_prb
./test_values_prb > test_values_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_values_prb."
  exit
fi
rm test_values_prb
#
echo "Program output written to test_values_prb_output.txt"
