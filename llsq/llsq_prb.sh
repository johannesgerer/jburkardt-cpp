#!/bin/bash
#
g++ -c -g -I/$HOME/include llsq_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling llsq_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ llsq_prb.o /$HOME/libcpp/$ARCH/llsq.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading llsq_prb.o."
  exit
fi
#
rm llsq_prb.o
#
mv a.out llsq_prb
./llsq_prb > llsq_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running llsq_prb."
  exit
fi
rm llsq_prb
#
echo "Program output written to llsq_prb_output.txt"
