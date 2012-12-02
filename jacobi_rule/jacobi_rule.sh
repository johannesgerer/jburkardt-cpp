#!/bin/bash
#
g++ -c -g jacobi_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ jacobi_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jacobi_rule.o"
  exit
fi
rm jacobi_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/jacobi_rule
#
echo "Executable installed as ~/bincpp/$ARCH/jacobi_rule"
