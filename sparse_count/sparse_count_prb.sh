#!/bin/bash
#
g++ -c -g -I/$HOME/include sparse_count_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_count_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sparse_count_prb.o /$HOME/libcpp/$ARCH/sparse_count.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_count_prb.o."
  exit
fi
#
rm sparse_count_prb.o
#
mv a.out sparse_count_prb
./sparse_count_prb > sparse_count_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparse_count_prb."
  exit
fi
rm sparse_count_prb
#
echo "Program output written to sparse_count_prb_output.txt"
