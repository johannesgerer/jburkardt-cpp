#!/bin/bash
#
g++ -c -g hermite_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ hermite_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_rule.o"
  exit
fi
rm hermite_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/hermite_rule
#
echo "Executable installed as ~/bincpp/$ARCH/hermite_rule"
