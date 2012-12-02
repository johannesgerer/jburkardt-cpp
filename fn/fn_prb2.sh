#!/bin/bash
#
g++ -c -g -I/$HOME/include fn_prb2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fn_prb2.cpp"
  exit
fi
rm compiler.txt
#
g++ fn_prb2.o /$HOME/libcpp/$ARCH/fn.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fn_prb2.o."
  exit
fi
#
rm fn_prb2.o
#
mv a.out fn_prb2
./fn_prb2 > fn_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fn_prb2."
  exit
fi
rm fn_prb2
#
echo "Program output written to fn_prb2_output.txt"
