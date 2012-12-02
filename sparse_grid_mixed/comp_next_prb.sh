#!/bin/bash
#
g++ -c -g -I/$HOME/include comp_next_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling comp_next_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ comp_next_prb.o /$HOME/libcpp/$ARCH/sparse_grid_mixed.o /$HOME/libcpp/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading comp_next_prb.o."
  exit
fi
#
rm comp_next_prb.o
#
mv a.out comp_next_prb
./comp_next_prb > comp_next_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running comp_next_prb."
  exit
fi
rm comp_next_prb
#
echo "Program output written to comp_next_prb_output.txt"
