#!/bin/bash
#
cp sphere_stereograph.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include sphere_stereograph.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_stereograph.cpp"
  exit
fi
rm compiler.txt
#
mv sphere_stereograph.o ~/libcpp/$ARCH/sphere_stereograph.o
#
echo "Library installed as ~/libcpp/$ARCH/sphere_stereograph.o"
