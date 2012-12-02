#!/bin/bash
#
g++ -c -g -I/$HOME/include ode_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ode_prb.cpp."
  exit
fi
rm compiler.txt
#
g++ ode_prb.o /$HOME/libcpp/$ARCH/ode.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ode_prb.o."
  exit
fi
#
rm ode_prb.o
#
mv a.out ode_prb
./ode_prb > ode_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ode_prb."
  exit
fi
rm ode_prb
#
echo "Program output written to ode_prb_output.txt"
