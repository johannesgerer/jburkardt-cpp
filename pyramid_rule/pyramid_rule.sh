#!/bin/bash
#
g++ -c -g pyramid_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_rule.cpp"
  exit
fi
rm compiler.txt
#
g++ pyramid_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_rule.o"
  exit
fi
rm pyramid_rule.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/pyramid_rule
#
echo "Executable installed as ~/bincpp/$ARCH/pyramid_rule"
