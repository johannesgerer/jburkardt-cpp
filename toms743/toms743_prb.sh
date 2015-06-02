#!/bin/bash
#
g++ -c -I/$HOME/include toms743_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling toms743_prb.cpp"
  exit
fi
#
g++ toms743_prb.o /$HOME/libcpp/$ARCH/toms743.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms743_prb.o."
  exit
fi
#
rm toms743_prb.o
#
mv a.out toms743_prb
./toms743_prb > toms743_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms743_prb."
  exit
fi
rm toms743_prb
#
echo "Program output written to toms743_prb_output.txt"
