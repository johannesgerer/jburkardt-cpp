#!/bin/bash
#
g++ -c -g -I/$HOME/include niederreiter_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling niederreiter_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ niederreiter_prb.o /$HOME/libcpp/$ARCH/niederreiter.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading niederreiter_prb.o."
  exit
fi
#
rm niederreiter_prb.o
#
mv a.out niederreiter_prb
./niederreiter_prb > niederreiter_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running niederreiter_prb."
  exit
fi
rm niederreiter_prb
#
echo "Program output written to niederreiter_prb_output.txt"
