#!/bin/bash
#
g++ -c -g -I/$HOME/include treepack_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling treepack_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ treepack_prb.o /$HOME/libcpp/$ARCH/treepack.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading treepack_prb.o."
  exit
fi
#
rm treepack_prb.o
#
mv a.out treepack_prb
./treepack_prb > treepack_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running treepack_prb."
  exit
fi
rm treepack_prb
#
echo "Program output written to treepack_prb_output.txt"
