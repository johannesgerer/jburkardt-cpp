#!/bin/bash
#
g++ -c gmsh_to_fem.cpp
if [ $? -ne 0 ]; then
  echo "Errors compiling gmsh_to_fem.cpp"
  exit
fi
#
g++ gmsh_to_fem.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gmsh_to_fem.o."
  exit
fi
rm gmsh_to_fem.o
#
mv a.out ~/bincpp/$ARCH/gmsh_to_fem
#
echo "Executable installed as ~/bincpp/$ARCH/gmsh_to_fem"
