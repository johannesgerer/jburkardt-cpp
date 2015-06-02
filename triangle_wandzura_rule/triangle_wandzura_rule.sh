#!/bin/bash
#
cp triangle_wandzura_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle_wandzura_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_wandzura_rule.cpp"
  exit
fi
#
mv triangle_wandzura_rule.o ~/libcpp/$ARCH/triangle_wandzura_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_wandzura_rule.o"
