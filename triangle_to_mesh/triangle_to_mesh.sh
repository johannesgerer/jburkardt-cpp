#!/bin/bash
#
g++ -c -g triangle_to_mesh.cpp >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_to_mesh.cpp"
  exit
fi
rm compiler.txt
#
g++ triangle_to_mesh.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_to_mesh.o."
  exit
fi
rm triangle_to_mesh.o
#
mv a.out ~/bincpp/$ARCH/triangle_to_mesh
#
echo "Executable installed as ~/bincpp/$ARCH/triangle_to_mesh"
