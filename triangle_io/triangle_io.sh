#!/bin/bash
#
cp triangle_io.hpp /$HOME/include
#
g++ -c triangle_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_io.cpp"
  exit
fi
#
mv triangle_io.o ~/libcpp/$ARCH/triangle_io.o
#
echo "Library installed as ~/libcpp/$ARCH/triangle_io.o"
