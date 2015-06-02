#!/bin/bash
#
g++ -c lorenz_ode.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lorenz_ode.cpp"
  exit
fi
rm compiler.txt
#
g++ lorenz_ode.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lorenz_ode.o"
  exit
fi
rm lorenz_ode.o
#
mv a.out ~/bincpp/$ARCH/lorenz_ode
#
echo "Executable installed as ~/bincpp/$ARCH/lorenz_ode"
