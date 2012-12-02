#!/bin/bash
#
g++ -c -g tet_mesh_l2q.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling tet_mesh_l2q.cpp"
  exit
fi
rm compiler.txt
#
g++ tet_mesh_l2q.o -lm
if [ $? -ne 0 ]; then
  echo "Error loading tet_mesh_l2q.o"
  exit
fi
rm tet_mesh_l2q.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/tet_mesh_l2q
#
echo "Executable installed as ~/bincpp/$ARCH/tet_mesh_l2q"
