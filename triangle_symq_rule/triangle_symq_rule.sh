#!/bin/bash
#
cp triangle_symq_rule.hpp /$HOME/include
#
g++ -c -I/$HOME/include triangle_symq_rule.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_symq_rule.cpp"
  exit
fi
#
mv triangle_symq_rule.o ~/libcpp/$ARCH/triangle_symq_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_symq_rule.o"
