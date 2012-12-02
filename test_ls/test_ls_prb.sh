#!/bin/bash
#
g++ -c -g -I/$HOME/include test_ls_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_ls_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ test_ls_prb.o /$HOME/libcpp/$ARCH/test_ls.o /$HOME/libcpp/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_ls_prb.o"
  exit
fi
#
rm test_ls_prb.o
#
mv a.out test_ls_prb
./test_ls_prb > test_ls_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_ls_prb."
  exit
fi
rm test_ls_prb
#
echo "Program output written to test_ls_prb_output.txt"
