#!/bin/bash
#
g++ -c string_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling string_simulation.cpp"
  exit
fi
rm compiler.txt
#
g++ string_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking string_simulation.o"
  exit
fi
#
rm string_simulation.o
mv a.out ~/bincpp/$ARCH/string_simulation
#
echo "Executable installed as ~/bincpp/$ARCH/string_simulation"
