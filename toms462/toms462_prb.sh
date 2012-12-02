#!/bin/bash
#
g++ -c -g -I/$HOME/include toms462_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms462_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ toms462_prb.o /$HOME/libcpp/$ARCH/toms462.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms462_prb.o."
  exit
fi
#
rm toms462_prb.o
#
mv a.out toms462_prb
./toms462_prb > toms462_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms462_prb."
  exit
fi
rm toms462_prb
#
echo "Program output written to toms462_prb_output.txt"
