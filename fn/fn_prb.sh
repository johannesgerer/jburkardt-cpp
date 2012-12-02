#!/bin/bash
#
g++ -c -g -I/$HOME/include fn_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fn_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ fn_prb.o /$HOME/libcpp/$ARCH/fn.o /$HOME/libcpp/$ARCH/test_values.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fn_prb.o."
  exit
fi
#
rm fn_prb.o
#
mv a.out fn_prb
./fn_prb > fn_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fn_prb."
  exit
fi
rm fn_prb
#
echo "Program output written to fn_prb_output.txt"
