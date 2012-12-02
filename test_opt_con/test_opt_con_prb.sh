#!/bin/bash
#
g++ -c -g -I/$HOME/include test_opt_con_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_opt_con_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ test_opt_con_prb.o /$HOME/libcpp/$ARCH/test_opt_con.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_opt_con_prb.o."
  exit
fi
#
rm test_opt_con_prb.o
#
mv a.out test_opt_con_prb
./test_opt_con_prb > test_opt_con_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_opt_con_prb."
  exit
fi
rm test_opt_con_prb
#
echo "Program output written to test_opt_con_prb_output.txt"
