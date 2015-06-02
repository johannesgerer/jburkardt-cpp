#!/bin/bash
#
g++ -c -I/$HOME/include naca_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling naca_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ naca_prb.o /$HOME/libcpp/$ARCH/naca.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading naca_prb.o."
  exit
fi
#
rm naca_prb.o
#
mv a.out naca_prb
./naca_prb > naca_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running naca_prb."
  exit
fi
rm naca_prb
#
echo "Program output written to naca_prb_output.txt"
