#!/bin/bash
#
g++ -c -g -I/$HOME/include toms655_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms655_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ toms655_prb.o /$HOME/libcpp/$ARCH/toms655.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms655_prb.o."
  exit
fi
#
rm toms655_prb.o
#
mv a.out toms655_prb
./toms655_prb > toms655_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms655_prb."
  exit
fi
rm toms655_prb
#
echo "Program output written to toms655_prb_output.txt"
