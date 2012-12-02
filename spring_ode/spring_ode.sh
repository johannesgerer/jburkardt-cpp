#!/bin/bash
#
g++ -c spring_ode.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode.cpp"
  exit
fi
rm compiler.txt
#
g++ spring_ode.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode.o"
  exit
fi
rm spring_ode.o
#
mv a.out ~/bincpp/$ARCH/spring_ode
#
echo "Executable installed as ~/bincpp/$ARCH/spring_ode"
