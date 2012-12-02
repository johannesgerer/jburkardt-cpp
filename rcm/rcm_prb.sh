#!/bin/bash
#
g++ -c -g -I/$HOME/include rcm_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rcm_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ rcm_prb.o /$HOME/libcpp/$ARCH/rcm.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rcm_prb.o."
  exit
fi
#
rm rcm_prb.o
#
mv a.out rcm_prb
./rcm_prb > rcm_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rcm_prb."
  exit
fi
rm rcm_prb
#
echo "Program output written to rcm_prb_output.txt"
