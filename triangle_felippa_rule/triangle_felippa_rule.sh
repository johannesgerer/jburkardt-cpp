#!/bin/bash
#
cp triangle_felippa_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle_felippa_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_felippa_rule.cpp"
  exit
fi
#
mv triangle_felippa_rule.o ~/libcpp/$ARCH/triangle_felippa_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_felippa_rule.o"
