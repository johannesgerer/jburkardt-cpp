#!/bin/bash
#
cp tet_mesh.hpp /$HOME/include
#
g++ -c -g -I /$HOME/include tet_mesh.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tet_mesh.cpp"
  exit
fi
rm compiler.txt
#
mv tet_mesh.o ~/libcpp/$ARCH/tet_mesh.o
#
echo "Library installed as ~/libcpp/$ARCH/tet_mesh.o"
