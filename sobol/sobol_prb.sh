#!/bin/bash
#
g++ -c -g -I/$HOME/include sobol_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sobol_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ sobol_prb.o /$HOME/libcpp/$ARCH/sobol.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sobol_prb.o."
  exit
fi
#
rm sobol_prb.o
#
mv a.out sobol_prb
./sobol_prb > sobol_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sobol_prb."
  exit
fi
rm sobol_prb
#
echo "Program output written to sobol_prb_output.txt"
