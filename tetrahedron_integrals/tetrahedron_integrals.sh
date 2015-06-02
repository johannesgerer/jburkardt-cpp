#!/bin/bash
#
cp tetrahedron_integrals.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include tetrahedron_integrals.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_integrals.cpp"
  exit
fi
rm compiler.txt
#
mv tetrahedron_integrals.o ~/libcpp/$ARCH/tetrahedron_integrals.o
#
echo "Library installed as ~/libcpp/$ARCH/tetrahedron_integrals.o"
