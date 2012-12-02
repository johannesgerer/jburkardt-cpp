#!/bin/bash
#
g++ -c -g gegenbauer_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gegenbauer_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ gegenbauer_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gegenbauer_rule.o"
  exit
fi
rm gegenbauer_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/gegenbauer_rule
#
echo "Executable installed as ~/bincpp/$ARCH/gegenbauer_rule"
