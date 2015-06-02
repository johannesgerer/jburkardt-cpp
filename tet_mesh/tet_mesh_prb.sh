#!/bin/bash
#
g++ -c -I/$HOME/include tet_mesh_prb.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling tet_mesh_prb.cpp"
  exit
fi
#
g++ tet_mesh_prb.o /$HOME/libcpp/$ARCH/tet_mesh.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tet_mesh_prb.o."
  exit
fi
#
rm tet_mesh_prb.o
#
mv a.out tet_mesh_prb
./tet_mesh_prb > tet_mesh_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running tet_mesh_prb."
  exit
fi
rm tet_mesh_prb
#
echo "Program output written to tet_mesh_prb_output.txt"
