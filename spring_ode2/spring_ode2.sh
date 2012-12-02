#!/bin/bash
#
g++ -c spring_ode2.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spring_ode2.cpp"
  exit
fi
rm compiler.txt
#
g++ spring_ode2.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spring_ode2.o"
  exit
fi
rm spring_ode2.o
#
mv a.out ~/bincpp/$ARCH/spring_ode2
#
echo "Executable installed as ~/bincpp/$ARCH/spring_ode2"
