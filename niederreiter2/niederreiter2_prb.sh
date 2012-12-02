#!/bin/bash
#
g++ -c -g -I/$HOME/include niederreiter2_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter2_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ niederreiter2_prb.o /$HOME/libcpp/$ARCH/niederreiter2.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading niederreiter2_prb.o."
  exit
fi
#
rm niederreiter2_prb.o
#
mv a.out niederreiter2_prb
./niederreiter2_prb > niederreiter2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running niederreiter2_prb."
  exit
fi
rm niederreiter2_prb
#
echo "Program output written to niederreiter2_prb_output.txt"
