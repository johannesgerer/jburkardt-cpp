#!/bin/bash
#
g++ -c -g ising_2d_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ising_2d_simulation.cpp"
  exit
fi
rm compiler.txt
#
g++ ising_2d_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ising_2d_simulation.o."
  exit
fi
#
rm ising_2d_simulation.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/ising_2d_simulation
#
echo "Executable installed as ~/bincpp/$ARCH/ising_2d_simulation"
