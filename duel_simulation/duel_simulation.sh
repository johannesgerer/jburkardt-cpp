#!/bin/bash
#
g++ -c -g duel_simulation.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling duel_simulation.cpp"
  exit
fi
rm compiler.txt
#
g++ duel_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading duel_simulation.o"
  exit
fi
rm duel_simulation.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/duel_simulation
#
echo "Executable installed as ~/bincpp/$ARCH/duel_simulation"
