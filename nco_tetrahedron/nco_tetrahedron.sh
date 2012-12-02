#!/bin/bash
#
cp nco_tetrahedron.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include nco_tetrahedron.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nco_tetrahedron.cpp"
  exit
fi
rm compiler.txt
#
mv nco_tetrahedron.o ~/libcpp/$ARCH/nco_tetrahedron.o
#
echo "Library installed as ~/libcpp/$ARCH/nco_tetrahedron.o"
