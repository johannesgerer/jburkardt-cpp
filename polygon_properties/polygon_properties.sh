#!/bin/bash
#
cp polygon_properties.hpp /$HOME/include
#
g++ -c -I/$HOME/include polygon_properties.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_properties.cpp"
  exit
fi
rm compiler.txt
#
mv polygon_properties.o ~/libcpp/$ARCH/polygon_properties.o
#
echo "Library installed as ~/libcpp/$ARCH/polygon_properties.o"
