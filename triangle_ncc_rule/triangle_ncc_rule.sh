#!/bin/bash
#
cp triangle_ncc_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle_ncc_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_ncc_rule.cpp"
  exit
fi
#
mv triangle_ncc_rule.o ~/libcpp/$ARCH/triangle_ncc_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_ncc_rule.o"
