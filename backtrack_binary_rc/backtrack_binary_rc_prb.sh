#!/bin/bash
#
g++ -c -g -I/$HOME/include backtrack_binary_rc_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling backtrack_binary_rc_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ backtrack_binary_rc_prb.o /$HOME/libcpp/$ARCH/backtrack_binary_rc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading backtrack_binary_rc_prb.o."
  exit
fi
#
rm backtrack_binary_rc_prb.o
#
mv a.out backtrack_binary_rc_prb
./backtrack_binary_rc_prb > backtrack_binary_rc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running backtrack_binary_rc_prb."
  exit
fi
rm backtrack_binary_rc_prb
#
echo "Program output written to backtrack_binary_rc_prb_output.txt"
