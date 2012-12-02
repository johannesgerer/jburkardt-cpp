#!/bin/bash
#
g++ -c -g -I/$HOME/include floyd_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling floyd_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ floyd_prb.o /$HOME/libcpp/$ARCH/floyd.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading floyd_prb.o."
  exit
fi
#
rm floyd_prb.o
#
mv a.out floyd_prb
./floyd_prb > floyd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running floyd_prb."
  exit
fi
rm floyd_prb
#
echo "Program output written to floyd_prb_output.txt"
