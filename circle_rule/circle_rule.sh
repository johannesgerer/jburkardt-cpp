#!/bin/bash
#
cp circle_rule.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include circle_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_rule.cpp"
  exit
fi
rm compiler.txt
#
mv circle_rule.o ~/libcpp/$ARCH/circle_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/circle_rule.o"
