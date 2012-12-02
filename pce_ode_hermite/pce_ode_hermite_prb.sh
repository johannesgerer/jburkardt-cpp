#!/bin/bash
#
g++ -c -g -I/$HOME/include pce_ode_hermite_prb.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pce_ode_hermite_prb.cpp"
  exit
fi
rm compiler.txt
#
g++ pce_ode_hermite_prb.o /$HOME/libcpp/$ARCH/pce_ode_hermite.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pce_ode_hermite_prb.o"
  exit
fi
#
rm pce_ode_hermite_prb.o
#
mv a.out pce_ode_hermite_prb
./pce_ode_hermite_prb > pce_ode_hermite_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pce_ode_hermite_prb."
  exit
fi
rm pce_ode_hermite_prb
#
echo "Program output written to pce_ode_hermite_prb_output.txt"
