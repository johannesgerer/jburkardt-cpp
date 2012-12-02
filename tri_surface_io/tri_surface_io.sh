#!/bin/bash
#
cp tri_surface_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include tri_surface_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tri_surface_io.cpp"
  exit
fi
rm compiler.txt
#
mv tri_surface_io.o ~/libcpp/$ARCH/tri_surface_io.o
#
echo "Library installed as ~/libcpp/$ARCH/tri_surface_io.o"
