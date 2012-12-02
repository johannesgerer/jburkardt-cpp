#!/bin/bash
#
g++ -c -g -I/$HOME/include test_tri_int_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_tri_int_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ test_tri_int_prb.o /$HOME/libcpp/$ARCH/test_tri_int.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_tri_int_prb.o."
  exit
fi
#
rm test_tri_int_prb.o
#
mv a.out test_tri_int_prb
./test_tri_int_prb > test_tri_int_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_tri_int_prb."
  exit
fi
rm test_tri_int_prb
#
echo "Program output written to test_tri_int_prb_output.txt"
