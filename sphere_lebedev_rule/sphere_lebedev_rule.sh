#!/bin/bash
#
cp sphere_lebedev_rule.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sphere_lebedev_rule.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_lebedev_rule.cpp"
  exit
fi
rm compiler.txt
#
mv sphere_lebedev_rule.o ~/libcpp/$ARCH/sphere_lebedev_rule.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_lebedev_rule.o"
