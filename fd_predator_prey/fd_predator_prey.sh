#!/bin/bash
#
g++ -c -g fd_predator_prey.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd_predator_prey.cpp"
  exit
fi
rm compiler.txt
#
g++ fd_predator_prey.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd_predator_prey.o"
  exit
fi
rm fd_predator_prey.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/fd_predator_prey
#
echo "Executable installed as ~/bincpp/$ARCH/fd_predator_prey"
