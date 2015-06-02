#!/bin/bash
#
cp hypersphere_properties.hpp /$HOME/include
#
g++ -c -I/$HOME/include hypersphere_properties.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling hypersphere_properties.cpp"
  exit
fi
#
mv hypersphere_properties.o ~/libcpp/$ARCH/hypersphere_properties.o
#
echo "Library installed as ~/libcpp/$ARCH/hypersphere_properties.o"
