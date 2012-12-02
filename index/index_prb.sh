#!/bin/bash
#
g++ -c -g -I/$HOME/include index_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling index_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ index_prb.o /$HOME/libcpp/$ARCH/index.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading index_prb.o."
  exit
fi
#
rm index_prb.o
#
mv a.out index_prb
./index_prb > index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running index_prb."
  exit
fi
rm index_prb
#
echo "Program output written to index_prb_output.txt"
