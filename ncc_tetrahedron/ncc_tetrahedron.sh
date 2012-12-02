#!/bin/bash
#
cp ncc_tetrahedron.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include ncc_tetrahedron.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ncc_tetrahedron.cpp"
  exit
fi
rm compiler.txt
#
mv ncc_tetrahedron.o ~/libcpp/$ARCH/ncc_tetrahedron.o
#
echo "Library installed as ~/libcpp/$ARCH/ncc_tetrahedron.o"
