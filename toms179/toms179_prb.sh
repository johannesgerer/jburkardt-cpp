#!/bin/bash
#
g++ -c -g -I/$HOME/include toms179_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms179_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ toms179_prb.o /$HOME/libcpp/$ARCH/toms179.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms179_prb.o."
  exit
fi
#
rm toms179_prb.o
#
mv a.out toms179_prb
./toms179_prb > toms179_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms179_prb."
  exit
fi
rm toms179_prb
#
echo "Program output written to toms179_prb_output.txt"
