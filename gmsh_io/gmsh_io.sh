#!/bin/bash
#
cp gmsh_io.hpp /$HOME/include
#
g++ -c -I /$HOME/include gmsh_io.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling gmsh_io.cpp"
  exit
fi
#
mv gmsh_io.o ~/libcpp/$ARCH/gmsh_io.o
#
echo "Library installed as ~/libcpp/$ARCH/gmsh_io.o"
