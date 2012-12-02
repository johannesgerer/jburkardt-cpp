#!/bin/bash
#
g++ -c -g tanh_sinh_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tanh_sinh_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ tanh_sinh_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tanh_sinh_rule.o"
  exit
fi
rm tanh_sinh_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/tanh_sinh_rule
#
echo "Executable installed as ~/bincpp/$ARCH/tanh_sinh_rule"
