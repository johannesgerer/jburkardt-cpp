#!/bin/bash
#
cp power_method.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include power_method.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling power_method.cpp"
  exit
fi
rm compiler.txt
#
mv power_method.o ~/libcpp/$ARCH/power_method.o
#
echo "Library installed as ~/libcpp/$ARCH/power_method.o"
