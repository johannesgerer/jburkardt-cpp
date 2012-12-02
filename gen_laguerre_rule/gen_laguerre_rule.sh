#!/bin/bash
#
g++ -c -g gen_laguerre_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gen_laguerre_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ gen_laguerre_rule.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gen_laguerre_rule.o"
  exit
fi
rm gen_laguerre_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/gen_laguerre_rule
#
echo "Executable installed as ~/bincpp/$ARCH/gen_laguerre_rule"
