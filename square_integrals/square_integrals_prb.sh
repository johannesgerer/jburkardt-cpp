#!/bin/bash
#
g++ -c -g -I/$HOME/include square_integrals_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square_integrals_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ square_integrals_prb.o /$HOME/libcpp/$ARCH/square_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_integrals_prb.o."
  exit
fi
#
rm square_integrals_prb.o
#
mv a.out square_integrals_prb
./square_integrals_prb > square_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_integrals_prb."
  exit
fi
rm square_integrals_prb
#
echo "Program output written to square_integrals_prb_output.txt"
