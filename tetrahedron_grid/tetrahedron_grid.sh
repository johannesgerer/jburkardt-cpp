#!/bin/bash
#
cp tetrahedron_grid.H /$HOME/include
#
g++ -c -g tetrahedron_grid.C >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_grid.C."
  exit
fi
rm compiler.txt
#
mv tetrahedron_grid.o ~/libcpp/$ARCH/tetrahedron_grid.o
#
echo "Library installed as ~/libcpp/$ARCH/tetrahedron_grid.o"
