#!/bin/bash
#
cp mesh_io.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include mesh_io.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mesh_io.cpp"
  exit
fi
rm compiler.txt
#
mv mesh_io.o ~/libcpp/$ARCH/mesh_io.o
#
echo "Library installed as ~/libcpp/$ARCH/mesh_io.o"
