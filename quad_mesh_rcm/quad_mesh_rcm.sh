#!/bin/bash
#
g++ -c -g quad_mesh_rcm.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_mesh_rcm.cpp"
  exit
fi
rm compiler.txt
#
g++ quad_mesh_rcm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quad_mesh_rcm.o."
  exit
fi
#
rm quad_mesh_rcm.o
#
chmod ugo+x a.out
mv a.out ~/bincpp/$ARCH/quad_mesh_rcm
#
echo "Executable installed as ~/bincpp/$ARCH/quad_mesh_rcm"
