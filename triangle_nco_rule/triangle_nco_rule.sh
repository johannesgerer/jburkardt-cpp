#!/bin/bash
#
cp triangle_nco_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle_nco_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_nco_rule.cpp"
  exit
fi
#
mv triangle_nco_rule.o ~/libcpp/$ARCH/triangle_nco_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_nco_rule.o"
