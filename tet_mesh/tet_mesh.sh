#!/bin/bash
#
cp tet_mesh.hpp /$HOME/include
#
g++ -c -I /$HOME/include tet_mesh.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling tet_mesh.cpp"
  exit
fi
#
mv tet_mesh.o ~/libcpp/$ARCH/tet_mesh.o
#
echo "Library installed as ~/libcpp/$ARCH/tet_mesh.o"
