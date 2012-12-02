#!/bin/bash
#
g++ -c -g -I$HOME/include product_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling product_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ product_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading product_rule.o"
  exit
fi
rm product_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/product_rule
#
echo "Executable installed as ~/bincpp/$ARCH/product_rule"
