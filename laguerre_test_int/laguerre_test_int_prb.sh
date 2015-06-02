#!/bin/bash
#
g++ -c -I/$HOME/include laguerre_test_int_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_test_int_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ laguerre_test_int_prb.o /$HOME/libcpp/$ARCH/laguerre_test_int.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_test_int_prb.o."
  exit
fi
#
rm laguerre_test_int_prb.o
#
mv a.out laguerre_test_int_prb
./laguerre_test_int_prb > laguerre_test_int_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running laguerre_test_int_prb."
  exit
fi
rm laguerre_test_int_prb
#
echo "Program output written to laguerre_test_int_prb_output.txt"
