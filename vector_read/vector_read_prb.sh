#!/bin/bash
#
g++ -c -g -I/$HOME/include vector_read_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vector_read_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ vector_read_prb.o /$HOME/libcpp/$ARCH/vector_read.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vector_read_prb.o."
  exit
fi
#
rm vector_read_prb.o
#
mv a.out vector_read_prb
./vector_read_prb < vector_read_prb_input.txt > vector_read_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vector_read_prb."
  exit
fi
rm vector_read_prb
#
echo "Program output written to vector_read_prb_output.txt"
