#!/bin/bash
#
cp triangle_dunavant_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle_dunavant_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_dunavant_rule.cpp"
  exit
fi
#
mv triangle_dunavant_rule.o ~/libcpp/$ARCH/triangle_dunavant_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_dunavant_rule.o"
