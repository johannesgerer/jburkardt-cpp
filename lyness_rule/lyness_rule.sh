#!/bin/bash
#
cp lyness_rule.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include lyness_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lyness_rule.cpp"
  exit
fi
rm compiler.txt
#
mv lyness_rule.o ~/libcpp/$ARCH/lyness_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/lyness_rule.o"
