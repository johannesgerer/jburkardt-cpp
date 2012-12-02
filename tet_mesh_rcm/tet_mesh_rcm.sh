#!/bin/bash
#
g++ -c -g tet_mesh_rcm.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling tet_mesh_rcm.cpp"
  exit
fi
rm compiler.txt
#
g++ tet_mesh_rcm.o -lm
if [ $? -ne 0 ]; then
  echo "Error loading tet_mesh_rcm.o"
  exit
fi
rm tet_mesh_rcm.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/tet_mesh_rcm
#
echo "Executable installed as ~/bincp/$ARCH/tet_mesh_rcm"
