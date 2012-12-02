#!/bin/bash
#
g++ -c -g chebyshev1_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev1_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ chebyshev1_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev1_rule.o"
  exit
fi
rm chebyshev1_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/chebyshev1_rule
#
echo "Executable installed as ~/bincpp/$ARCH/chebyshev1_rule"
