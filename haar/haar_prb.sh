#!/bin/bash
#
g++ -c -g -I/$HOME/include haar_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling haar_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ haar_prb.o /$HOME/libcpp/$ARCH/haar.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading haar_prb.o."
  exit
fi
#
rm haar_prb.o
#
mv a.out haar_prb
./haar_prb > haar_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running haar_prb."
  exit
fi
rm haar_prb
#
echo "Program output written to haar_prb_output.txt"
