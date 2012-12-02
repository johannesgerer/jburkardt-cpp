#!/bin/bash
#
g++ -c -g patterson_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling patterson_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ patterson_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading patterson_rule.o"
  exit
fi
rm patterson_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/patterson_rule
#
echo "Executable installed as ~/bincpp/$ARCH/patterson_rule"
