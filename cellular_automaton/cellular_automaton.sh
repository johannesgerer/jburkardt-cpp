#!/bin/bash
#
g++ -c -g cellular_automaton.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cellular_automaton.cpp"
  exit
fi
rm compiler.txt
#
g++ cellular_automaton.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cellular_automaton.o"
  exit
fi
rm cellular_automaton.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/cellular_automaton
#
echo "Program installed as ~/bincpp/$ARCH/cellular_automaton"
