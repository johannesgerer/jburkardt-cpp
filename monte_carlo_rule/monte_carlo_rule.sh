#!/bin/bash
#
g++ -c -g monte_carlo_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling monte_carlo_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ monte_carlo_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading monte_carlo_rule.o"
  exit
fi
rm monte_carlo_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/monte_carlo_rule
#
echo "Executable installed as ~/bincpp/$ARCH/monte_carlo_rule"
