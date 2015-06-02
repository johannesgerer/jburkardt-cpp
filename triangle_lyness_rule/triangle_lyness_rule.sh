#!/bin/bash
#
cp triangle_lyness_rule.hpp /$HOME/include
#
g++ -c -I /$HOME/include triangle_lyness_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_lyness_rule.cpp"
  exit
fi
#
mv triangle_lyness_rule.o ~/libcpp/$ARCH/triangle_lyness_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_lyness_rule.o"
