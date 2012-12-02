#!/bin/bash
#
g++ -c -g -I/$HOME/include kmeans_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kmeans_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ kmeans_prb.o /$HOME/libcpp/$ARCH/kmeans.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading kmeans_prb.o."
  exit
fi
#
rm kmeans_prb.o
#
mv a.out kmeans_prb
./kmeans_prb > kmeans_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running kmeans_prb."
  exit
fi
rm kmeans_prb
#
echo "Program output written to kmeans_prb_output.txt"
