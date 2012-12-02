#!/bin/bash
#
cp gm_rule.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include gm_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gm_rule.cpp"
  exit
fi
rm compiler.txt
#
mv gm_rule.o ~/libcpp/$ARCH/gm_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/gm_rule.o"
