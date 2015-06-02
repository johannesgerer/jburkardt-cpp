#!/bin/bash
#
g++ -c -I/$HOME/include hermite_test_int_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_test_int_prb.cpp"
  exit
fi
#
g++ hermite_test_int_prb.o /$HOME/libcpp/$ARCH/hermite_test_int.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_test_int_prb.o."
  exit
fi
#
rm hermite_test_int_prb.o
#
mv a.out hermite_test_int_prb
./hermite_test_int_prb > hermite_test_int_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hermite_test_int_prb."
  exit
fi
rm hermite_test_int_prb
#
echo "Program output written to hermite_test_int_prb_output.txt"
