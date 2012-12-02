#!/bin/bash
#
g++ -c -g sde_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sde_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sde_prb.o /$HOME/libcpp/$ARCH/sde.o /$HOME/libcpp/$ARCH/qr_solve.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sde_prb.o."
  exit
fi
#
rm sde_prb.o
#
mv a.out sde_prb
./sde_prb > sde_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sde_prb."
  exit
fi
rm sde_prb
#
echo "Program output written to sde_prb_output.txt"
