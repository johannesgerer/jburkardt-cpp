#!/bin/bash
#
cp quad_mesh.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include quad_mesh.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_mesh.cpp"
  exit
fi
rm compiler.txt
#
mv quad_mesh.o ~/libcpp/$ARCH/quad_mesh.o
#
echo "Library installed as ~/libcpp/$ARCH/quad_mesh.o"
